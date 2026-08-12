import argparse
from functools import lru_cache

from goatools.anno.gaf_reader import GafReader
from goatools.obo_parser import GODag
from goatools.semantic import TermCounts, lin_sim

from utils import Ontology

# --- Input / output paths ---
GO_RELEASE = "2025-07-22"                           # fetch with: bash ../go_data/download.sh
GO_OBO_PATH = f"../go_data/{GO_RELEASE}/go.obo"     # GO DAG (must match the release the annotations/IC were computed against)
UNIPROT_GAF = f"../go_data/{GO_RELEASE}/filtered_goa_uniprot_all_noiea.gaf"  # background annotation corpus used for IC

# --- Run as 'python calc_MF_Lin_set_sim.py <input_path> <output_path> [Optionaal: --go-obo <go_obo_path> --gaf <gaf_path>'] ---
parser = argparse.ArgumentParser(description="Compute a BMA Lin-similarity matrix for a set of proteins.")
parser.add_argument("input_path", help="TSV: Protein<TAB>GO1<TAB>GO2... (optional header row starting with 'Protein')")
parser.add_argument("output_path", help="Path to write the labelled similarity matrix TSV")
parser.add_argument("--go-obo", default=GO_OBO_PATH, help="Path to the GO OBO file")
parser.add_argument("--gaf", default=UNIPROT_GAF, help="Path to the background GOA GAF file (MF namespace)")
args = parser.parse_args()

# --- Load GO DAG, background annotation corpus, and term counts (for IC) ---
godag = GODag(args.go_obo)
gene2gos = GafReader(args.gaf).get_id2gos(namespace="MF")
termcounts = TermCounts(godag, gene2gos)

# --- Load ontology for most-specific-term filtering ---
ont = Ontology(args.go_obo, with_rels=True)


# --- Parse input: Protein<TAB>GO1<TAB>GO2... ---
def read_protein_annotations(path):
    annotations = {}
    with open(path) as f:
        for line in f:
            parts = [p.strip() for p in line.rstrip("\n").split("\t")]
            if not parts[0] or parts[0].startswith("Protein"):
                continue
            protein = parts[0]
            go_ids = list(dict.fromkeys(p for p in parts[1:] if p))
            valid_go_ids = [g for g in go_ids if g in godag]
            for g in set(go_ids) - set(valid_go_ids):
                print(f"Warning: {protein} annotation {g} not found in GO DAG, skipping")
            annotations[protein] = valid_go_ids
    return annotations


# --- Keep only the most specific (leaf) terms per protein's annotation set ---
def filter_most_specific(go_ids):
    go_set = set(go_ids)
    for go_id in list(go_set):
        ancestors = ont.get_ancestors(go_id)
        ancestors.discard(go_id)
        go_set -= ancestors
    return go_set


# --- Cached pairwise Lin similarity between GO terms ---
@lru_cache(maxsize=None)
def term_sim(go_a, go_b):
    if go_a not in godag or go_b not in godag:
        return 0.0
    sim = lin_sim(go_a, go_b, godag, termcounts)
    return sim if sim is not None else 0.0


# --- BMA (Best Match Average) groupwise similarity ---
def bma_sim(annots_a, annots_b):
    if not annots_a or not annots_b:
        return 0.0
    total = sum(max(term_sim(a, b) for b in annots_b) for a in annots_a)
    total += sum(max(term_sim(a, b) for a in annots_a) for b in annots_b)
    return total / (len(annots_a) + len(annots_b))


# --- Build the N x N similarity matrix (symmetric, upper triangle only) ---
annotations = read_protein_annotations(args.input_path)
for protein, go_ids in annotations.items():
    annotations[protein] = filter_most_specific(go_ids)

proteins = list(annotations.keys())
n = len(proteins)
matrix = [[0.0] * n for _ in range(n)]

for i in range(n):
    for j in range(i, n):
        sim = bma_sim(annotations[proteins[i]], annotations[proteins[j]])
        matrix[i][j] = sim
        matrix[j][i] = sim

# --- Write labelled matrix ---
with open(args.output_path, "w") as out:
    out.write("\t" + "\t".join(proteins) + "\n")
    for i, protein in enumerate(proteins):
        row = [protein] + [f"{matrix[i][j]:.4f}" for j in range(n)]
        out.write("\t".join(row) + "\n")

print(f"Done. Matrix written to {args.output_path}")
