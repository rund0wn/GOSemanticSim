# SemanticSimilarity

GO-based semantic similarity between proteins, using Lin's measure with
Best-Match Average (BMA) aggregation over annotation sets.

Two independent implementations:

| Directory | Stack | What it does |
|---|---|---|
| `Lins_Similarity/` | Python + [goatools](https://github.com/tanghaibao/goatools) | Computes a pairwise functional similarity matrix (Lin's similarity + BMA) from protein → GO term annotations |
| `groovy_semantics/` | Groovy + [SML / slib](https://www.semantic-measures-library.org/) | Information content and pairwise/set similarity over the GO DAG |


## Quick run :arrow_forward:
Generate a pairwise functional similarity matrix from protein GO term annotations.

1. Download reference data (`.gaf` and `.obo`)
```bash
bash go_data/download.sh            # DEFAULT (r2025-07-22)
```
2. Install dependencies
```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```
3. Generate functional similarity matrix (Lin's + BMA)
```bash
cd Lins_Similarity
python calc_MF_Lin_set_sim.py {INPUT}.tsv {OUTPUT}.tsv
```

## Technical Details :gear:

### Reference data

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

#### More information
##### Why these specific data files?
- Both files come from the same release directory http://release.geneontology.org/, and that is not incidental: information content is computed from the UniProt database annotations `.gaf` against the GO DAG `go.obo`. This ensures that the annotations and the DAG will include the exact same terms, and thus all terms will have a valid IC value.
- The `noiea` gaf does not include `IEA` (electronically inferred) annotations, as there are not "high quality".

##### What if I want to use another predictor/GO version?
If you use a tool/model to generate the annotations, use the *same* release used by the tool/model, which you generally find in their Github repo or publication. If you use database annotations, try to closely match the date of the database release. Databases like UniProt align their releases with GO releases. After identifying the correct version:

1. Download data from the required release:

```bash
bash go_data/download.sh {YYYY-MM-DD} # Date of other GO release
```

2. Either edit `GO_RELEASE` in `calc_MF_Lin_set_sim.py` OR pass `--go-obo` / `--gaf` explicitly.

---

### IC table
To generate a GO term → information content (IC) table:

```bash
cd Lins_Similarity
python build_ic_table.py {OUTPUT}.tsv                              # all three sub-ontologies, every term
python build_ic_table.py {OUTPUT}.tsv --terms ../terms_deepFRI2    # only terms predictable by deepFRI2
python build_ic_table.py {OUTPUT}.tsv --namespace BP --covered-only
```

| Option | Effect |
|---|---|
| *(none)* | Every non-obsolete term in MF, BP and CC (39,906 terms) |
| `--terms {DIR}` | Restrict to the GO IDs in the specified directory's `*.json` files |
| `--namespace {MF,BP,CC}` | Restrict to one sub-ontology |
| `--covered-only` | Drop terms without annotations |

Columns: `GO_ID`, `name`, `namespace`, `depth`, `n_annotated`, `IC`, `IC_topology`.

**Notes:**
- `IC` = `-log(n_annotated / corpus_size)` where the corpus is the reference data (`.gaf`). The same counts `calc_MF_Lin_set_sim.py` uses.
    > [!NOTE]
    > Each sub-ontology has its own root and corpus size, so IC is computed per sub-ontology.
- `IC = 0` means *maximally general* (the sub-ontology root), but `goatools` also returns `0.0` for terms absent from the corpus. The two are distinguished by `n_annotated`, where unsupported terms get a blank `IC`, not a zero.
- `IC_topology` = `1 - log(n_descendants + 1) / log(sub-ontology size)` ([Seco et al. 2004](https://dl.acm.org/doi/10.5555/3000001.3000272)). Derived from the DAG structure alone — how much of its sub-ontology a term subsumes — so it needs no corpus. The root scores 0, a leaf scores 1.
    > [!NOTE]
    > Because it ignores annotations, it is defined for **every** term, including the 23,070 that have no corpus support and therefore a blank `IC`.
- The two IC columns measure different things and disagree usefully (Pearson *r* ≈ 0.61–0.70, Spearman ρ ≈ 0.48–0.56 — correlated, not redundant):
    - **`IC` low, `IC_topology` high** — a well-studied leaf. `GO:0042803 protein homodimerization activity` has `IC=3.00` (very frequently annotated) but `IC_topology=1.00` (no descendants).
    - **`IC` high, `IC_topology` low** — a broad term nobody has annotated. `GO:0140852 histone ubiquitin ligase activity` hits the `IC` ceiling of 10.14 on a *single* annotated protein, yet `IC_topology=0.71` reveals it subsumes plenty of the DAG. Use `IC_topology` to spot these; raw `IC` makes them look maximally informative.
- Coverage in the `noiea` corpus is thin for rarely studied terms — 57% of MF, 57% of BP and 63% of CC terms have no annotation support.
    > [!NOTE] 
    > Restricting to the DeepFRI2 term space is much better: 851/858 MF, 3,796/3,994 BP, 572/615 CC.
- `ND` ("no biological data") and `NOT`-qualified annotations are excluded; this is also `goatools`' default for `get_id2gos`.

---

### Python usage

Described in [Quick run](#quick-run)
**Note:** `GUIDE_semantic_similarity.ipynb` walks through the method.

#### Input
`.tsv` of `Protein<TAB>GO1<TAB>GO2...` (an optional header row may
start with `Protein`).

#### Output
`.tsv` of pairwise similarity scores between proteins.

---

### Groovy usage

The scripts use `@Grab` to pull `slib-sml` and `gpars` on first run, so Groovy
and a JVM are the only prerequisites:

```bash
groovy groovy_semantics/Sim.groovy
```

- `IC.groovy`, `ICVectorSim.groovy` — information content, IC-vector similarity
- `Sim.groovy`, `Set_Sim.groovy`, `Set_Sim_Lin.groovy` — pairwise and set-level similarity
