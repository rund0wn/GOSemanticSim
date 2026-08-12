# SemanticSimilarity

GO-based semantic similarity between proteins, using Lin's measure with
Best-Match Average (BMA) aggregation over annotation sets.

Two independent implementations live here:

| Directory | Stack | What it does |
|---|---|---|
| `Lins_Similarity/` | Python + [goatools](https://github.com/tanghaibao/goatools) | Computes a labelled BMA Lin-similarity matrix from a TSV of protein → GO term annotations |
| `groovy_semantics/` | Groovy + [SML / slib](https://www.semantic-measures-library.org/) | Information content and pairwise/set similarity over the GO DAG |

## Data files

Download GO reference data:

```bash
bash go_data/download.sh            # DEFAULT (r2025-07-22)
bash go_data/download.sh {YYYY-MM-DD} # To download any other release
```

Expected output:

```
go_data/
├── download.sh
└── 2025-07-22/
    ├── go.obo                             # GO DAG
    └── filtered_goa_uniprot_all_noiea.gaf # UniProt GOA, IEA evidence removed
```

Both files come from the same release directory
http://release.geneontology.org/, and that is not incidental: information
content is computed from the UniProt database annotations `.gaf` against the GO DAG `go.obo`.

**Note: The `noiea` gaf does not include `IEA` (electronically inferred) annotations.**

To use a release other than the default, either edit `GO_RELEASE` in
`calc_MF_Lin_set_sim.py` or pass `--go-obo` / `--gaf` explicitly.

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
