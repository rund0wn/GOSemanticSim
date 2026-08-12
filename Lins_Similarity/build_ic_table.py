import argparse
import csv
import json
from pathlib import Path

from goatools.anno.gaf_reader import GafReader
from goatools.obo_parser import GODag
from goatools.semantic import TermCounts, get_info_content

# --- Input / output paths ---
GO_RELEASE = "2025-07-22"                           # fetch with: bash ../go_data/download.sh
GO_OBO_PATH = f"../go_data/{GO_RELEASE}/go.obo"     # GO DAG (must match the release the annotations were computed against)
UNIPROT_GAF = f"../go_data/{GO_RELEASE}/filtered_goa_uniprot_all_noiea.gaf"  # background annotation corpus
NAMESPACES = {"MF": "molecular_function", "BP": "biological_process", "CC": "cellular_component"}

# --- Run as 'python build_ic_table.py <output_path> [--terms <dir> --namespace MF --covered-only]' ---
parser = argparse.ArgumentParser(description="Build a table of GO terms and their information content.")
parser.add_argument("output_path", help="Path to write the IC table TSV")
parser.add_argument("--terms", help="Directory of per-sub-ontology .json term files, e.g. ../terms_deepFRI2")
parser.add_argument("--namespace", choices=list(NAMESPACES), help="Restrict to one branch (default: all three)")
parser.add_argument("--covered-only", action="store_true", help="Omit terms with no annotation support")
parser.add_argument("--go-obo", default=GO_OBO_PATH, help="Path to the GO OBO file")
parser.add_argument("--gaf", default=UNIPROT_GAF, help="Path to the background GOA GAF file")
args = parser.parse_args()

# --- Load GO DAG and background corpus ---
# get_id2gos already drops ND ("no biological data") and NOT-qualified annotations,
# so counts here match those used by calc_MF_Lin_set_sim.py.
godag = GODag(args.go_obo)
gaf = GafReader(args.gaf)
namespaces = [args.namespace] if args.namespace else list(NAMESPACES)

# --- Optional term restriction: {"0": "GO:...", ...} per sub-ontology ---
# Each term's branch comes from the DAG, not the filename, so the files need no
# particular naming and a term filed under the wrong one still lands correctly.
selected = None
if args.terms:
    files = sorted(Path(args.terms).glob("*.json"))
    if not files:
        parser.error(f"no .json term files found in {args.terms}")
    wanted = [go_id for f in files for go_id in json.load(open(f)).values()]
    missing = sorted({t for t in wanted if t not in godag})
    selected = {godag[t].id for t in wanted if t in godag}  # alt_ids -> primary
    print(f"Terms read : {len(wanted):,} from {len(files)} file(s), {len(selected):,} unique")
    if missing:
        print(f"WARNING: {len(missing)} term(s) not in GO {GO_RELEASE}, skipped: "
              + ", ".join(missing[:5]) + (" ..." if len(missing) > 5 else ""))

# --- Compute IC per branch and collect rows ---
# IC is relative to its own branch root and corpus size, so each namespace needs
# its own TermCounts. get_info_content returns 0.0 both for a root and for terms
# absent from the corpus, so n_annotated is what tells those two apart.
rows, summary = [], []
for ns in namespaces:
    termcounts = TermCounts(godag, gaf.get_id2gos(namespace=ns))
    branch = sorted({o.id for o in godag.values()
                     if o.namespace == NAMESPACES[ns] and not o.is_obsolete})
    # Normalise against the whole branch, not the selected subset, so IC_norm stays
    # comparable across runs that pass different --terms lists.
    max_ic = max((get_info_content(t, termcounts) for t in branch), default=0.0)
    go_ids = [t for t in branch if t in selected] if selected is not None else branch

    uncovered = 0
    for go_id in go_ids:
        count = termcounts.get_count(go_id)
        if not count:
            uncovered += 1
            if args.covered_only:
                continue
        ic = get_info_content(go_id, termcounts)
        rows.append({
            "GO_ID": go_id,
            "name": godag[go_id].name,
            "namespace": ns,
            "depth": godag[go_id].depth,
            "n_annotated": count,
            "IC": f"{ic:.6f}" if count else "",
            "IC_norm": f"{ic / max_ic:.6f}" if count and max_ic else "",
        })
    summary.append((ns, termcounts.get_total_count(NAMESPACES[ns]), len(go_ids), uncovered, max_ic))

# --- Write the table ---
with open(args.output_path, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"\n{'branch':<8}{'corpus':>10}{'terms':>10}{'no support':>12}{'max IC':>9}")
for ns, total, n_terms, uncovered, max_ic in summary:
    print(f"{ns:<8}{total:>10,}{n_terms:>10,}{uncovered:>12,}{max_ic:>9.4f}")
print(f"\n{len(rows):,} rows written to {args.output_path}"
      + (" (unsupported terms omitted)" if args.covered_only else " (unsupported terms have blank IC)"))
