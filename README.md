# SemanticSimilarity

GO-based semantic similarity between proteins, using Lin's measure with
Best-Match Average (BMA) aggregation over annotation sets.

Two independent implementations live here:

| Directory | Stack | What it does |
|---|---|---|
| `Lins_Similarity/` | Python + [goatools](https://github.com/tanghaibao/goatools) | Computes a labelled BMA Lin-similarity matrix from a TSV of protein → GO term annotations |
| `groovy_semantics/` | Groovy + [SML / slib](https://www.semantic-measures-library.org/) | Information content and pairwise/set similarity over the GO DAG |

## Data files (not in this repo)

`go_data/` is gitignored — it holds ~96 MB of public reference data that you
need to download yourself. Recreate the layout below:

```
go_data/
├── go.obo                                   # current GO release
├── 2019-06-09/
│   └── go.obo                               # pinned 2019-06-09 GO release
└── filtered_goa_uniprot_all_noiea.gaf       # UniProt GOA, IEA evidence removed
```

- **GO ontology** — https://geneontology.org/docs/download-ontology/
  (archived releases: http://release.geneontology.org/)
- **UniProt GOA** — https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/
  The `.gaf` here has been filtered to drop `IEA` (electronically inferred)
  annotations, so it is an experimental-evidence background corpus.

The GO release used to build the annotation corpus and the release used to
compute information content must match, which is why `2019-06-09/go.obo` is
pinned separately from the current `go.obo`.

## Python usage

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

```bash
cd Lins_Similarity
python calc_MF_Lin_set_sim.py test_in.tsv test_out.tsv
```

Input is a TSV of `Protein<TAB>GO1<TAB>GO2...` (an optional header row may
start with `Protein`). Output is a labelled square similarity matrix.

Override the defaults with `--go-obo <path>` and `--gaf <path>`.

`GUIDE_semantic_similarity.ipynb` walks through the method; `utils.py` holds a
standalone `Ontology` class (OBO parsing, ancestor traversal, IC calculation)
that does not depend on goatools.

## Groovy usage

The scripts use `@Grab` to pull `slib-sml` and `gpars` on first run, so Groovy
and a JVM are the only prerequisites:

```bash
groovy groovy_semantics/Sim.groovy
```

- `IC.groovy`, `ICVectorSim.groovy` — information content, IC-vector similarity
- `Sim.groovy`, `Set_Sim.groovy`, `Set_Sim_Lin.groovy` — pairwise and set-level similarity
