# SemanticSimilarity

GO-based semantic similarity between proteins, using Lin's measure with
Best-Match Average (BMA) aggregation over annotation sets.

Two independent implementations:

| Directory | Stack | What it does |
|---|---|---|
| `Lins_Similarity/` | Python + [goatools](https://github.com/tanghaibao/goatools) | Computes a pairwise functional similarity matrix (Lin's similarity + BMA) from protein → GO term annotations |
| `groovy_semantics/` | Groovy + [SML / slib](https://www.semantic-measures-library.org/) | Information content and pairwise/set similarity over the GO DAG |

---

## Data files (Do this first)

Download GO reference data:

```bash
bash go_data/download.sh            # DEFAULT (r2025-07-22)
```

Expected output:

```
go_data/
├── download.sh
└── 2025-07-22/
    ├── go.obo                             # GO DAG
    └── filtered_goa_uniprot_all_noiea.gaf # UniProt GOA, IEA evidence removed
```

### More information
#### Why these specific data files?
- Both files come from the same release directory http://release.geneontology.org/, and that is not incidental: information content is computed from the UniProt database annotations `.gaf` against the GO DAG `go.obo`. This ensures that the annotations and the DAG will include the exact same terms, and thus all terms will have a valid IC value.
- The `noiea` gaf does not include `IEA` (electronically inferred) annotations, as there are not "high quality".

#### What if I want to use another predictor/GO version?
If you use a tool/model to generate the annotations, use the *same* release used by the tool/model, which you generally find in their Github repo or publication. If you use database annotations, try to closely match the date of the database release. Databases like UniProt align their releases with GO releases. After identifying the correct version:

1. Download data from the required release:

```bash
bash go_data/download.sh {YYYY-MM-DD} # Date of other GO release
```

2. Either edit `GO_RELEASE` in `calc_MF_Lin_set_sim.py` OR pass `--go-obo` / `--gaf` explicitly.


---

## Python usage

1. Create virtual environment

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

2. Generate functional similarity matrix (Lin's + BMA)

```bash
cd Lins_Similarity
python calc_MF_Lin_set_sim.py {INPUT}.tsv {OUTPUT}.tsv
```

**Note:** `GUIDE_semantic_similarity.ipynb` walks through the method.


### Input
`.tsv` of `Protein<TAB>GO1<TAB>GO2...` (an optional header row may
start with `Protein`).

### Output
`.tsv` of pairwise similarity scores between proteins.

3. Generate a GO term → information content table

```bash
cd Lins_Similarity
python build_ic_table.py {OUTPUT}.tsv                      # whole MF branch
python build_ic_table.py {OUTPUT}.tsv --terms my_terms.txt # only these GO IDs
python build_ic_table.py {OUTPUT}.tsv --namespace BP --covered-only
```

Columns: `GO_ID`, `name`, `namespace`, `depth`, `n_annotated`, `IC`, `IC_norm`.

`IC` is `-log(n_annotated / corpus_size)` over annotations propagated to
ancestors — the same counts `calc_MF_Lin_set_sim.py` uses, so the two agree.
`IC_norm` is always scaled by the branch-wide maximum, so it stays comparable
across runs that pass different `--terms` lists.

**Reading the table:** `IC = 0` means *maximally general* (the branch root), but
`goatools` also returns `0.0` for terms absent from the corpus. The two are
distinguished by `n_annotated` — unsupported terms get a blank `IC`, never a
zero. Roughly 57% of MF terms have no annotation support in the `noiea` corpus;
`--covered-only` drops them.

`ND` ("no biological data") and `NOT`-qualified annotations are excluded — this
is `goatools`' default in `get_id2gos`, not something either script adds.

---

## Groovy usage

The scripts use `@Grab` to pull `slib-sml` and `gpars` on first run, so Groovy
and a JVM are the only prerequisites:

```bash
groovy groovy_semantics/Sim.groovy
```

- `IC.groovy`, `ICVectorSim.groovy` — information content, IC-vector similarity
- `Sim.groovy`, `Set_Sim.groovy`, `Set_Sim_Lin.groovy` — pairwise and set-level similarity
