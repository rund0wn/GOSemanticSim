import argparse
import csv

from goatools.anno.gaf_reader import GafReader
from goatools.obo_parser import GODag
from goatools.semantic import TermCounts, get_info_content

# --- Input / output paths ---
GO_RELEASE = "2025-07-22"                           # fetch with: bash ../go_data/download.sh
GO_OBO_PATH = f"../go_data/{GO_RELEASE}/go.obo"     # GO DAG (must match the release the annotations were computed against)
UNIPROT_GAF = f"../go_data/{GO_RELEASE}/filtered_goa_uniprot_all_noiea.gaf"  # background annotation corpus
NAMESPACES = {"MF": "molecular_function", "BP": "biological_process", "CC": "cellular_component"}

# --- Run as 'python build_ic_table.py <output_path> [--namespace MF --terms <file> --covered-only]' ---
parser = argparse.ArgumentParser(description="Build a table of GO terms and their information content.")
parser.add_argument("output_path", help="Path to write the IC table TSV")
parser.add_argument("--namespace", default="MF", choices=sorted(NAMESPACES), help="GO branch to tabulate")
parser.add_argument("--go-obo", default=GO_OBO_PATH, help="Path to the GO OBO file")
parser.add_argument("--gaf", default=UNIPROT_GAF, help="Path to the background GOA GAF file")
parser.add_argument("--terms", help="Optional file of GO IDs (one per line) to restrict the table to")
parser.add_argument("--covered-only", action="store_true", help="Omit terms with no annotation support")
args = parser.parse_args()

# --- Load GO DAG and background corpus ---
# get_id2gos already drops ND ("no biological data") and NOT-qualified annotations,
# so counts here match those used by calc_MF_Lin_set_sim.py.
godag = GODag(args.go_obo)
gene2gos = GafReader(args.gaf).get_id2gos(namespace=args.namespace)
termcounts = TermCounts(godag, gene2gos)
total = termcounts.get_total_count(NAMESPACES[args.namespace])

# --- Select terms: the whole branch, or a caller-supplied list ---
branch = sorted({o.id for o in godag.values()
                 if o.namespace == NAMESPACES[args.namespace] and not o.is_obsolete})
if args.terms:
    with open(args.terms) as fh:
        wanted = [ln.split()[0] for ln in fh if ln.strip() and not ln.startswith("#")]
    missing = [t for t in wanted if t not in godag]
    # Normalise alt_ids to primary IDs, preserving the caller's order
    seen, go_ids = set(), []
    for t in wanted:
        if t in godag and godag[t].id not in seen:
            seen.add(godag[t].id)
            go_ids.append(godag[t].id)
    if missing:
        print(f"WARNING: {len(missing)} term(s) not in this GO release, skipped: {', '.join(missing[:5])}"
              + (" ..." if len(missing) > 5 else ""))
else:
    go_ids = branch

# --- Compute IC and write the table ---
# get_info_content returns 0.0 both for the root and for terms absent from the corpus,
# so n_annotated is what distinguishes "maximally general" from "no data".
# Normalise against the whole branch, not the selected subset, so IC_norm stays
# comparable across runs that pass different --terms lists.
max_ic = max((get_info_content(t, termcounts) for t in branch), default=0.0)
rows, uncovered = [], 0
for go_id in go_ids:
    term = godag[go_id]
    count = termcounts.get_count(go_id)
    if not count:
        uncovered += 1
        if args.covered_only:
            continue
    ic = get_info_content(go_id, termcounts)
    rows.append({
        "GO_ID": go_id,
        "name": term.name,
        "namespace": args.namespace,
        "depth": term.depth,
        "n_annotated": count,
        "IC": f"{ic:.6f}" if count else "",
        "IC_norm": f"{ic / max_ic:.6f}" if count and max_ic else "",
    })

with open(args.output_path, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"\nCorpus     : {total:,} annotated gene products ({args.namespace})")
print(f"Terms      : {len(rows):,} written, {uncovered:,} with no annotation support"
      + (" (omitted)" if args.covered_only else " (IC left blank)"))
print(f"IC range   : 0 - {max_ic:.4f}")
print(f"Done. Table written to {args.output_path}")
