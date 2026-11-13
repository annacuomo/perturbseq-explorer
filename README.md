# Perturb-seq explorer

An interactive, browser-based explorer for the genome-scale CRISPRi Perturb-seq
data of [Replogle et al. 2022, *Cell*](https://doi.org/10.1016/j.cell.2022.05.013).

Pick a perturbation → see which genes it changes most, and what those genes have in common.

New to the web side? Read [HOW_IT_WORKS.md](HOW_IT_WORKS.md) first.

## How the repo is organised

```
perturbseq-explorer/
├── pipeline/        Python. Downloads the published data, computes compact
│                    summaries, writes small JSON into app/data/.
├── app/             The website. Plain HTML + JS + Plotly; reads app/data/.
│   └── data/        The "contract" between the two halves (see below).
├── raw/             Downloaded .h5ad files. Git-ignored, never committed.
└── .github/         Auto-deploys app/ to GitHub Pages on every push.
```

The two halves never talk to each other directly. Python writes files; the
browser reads them. If you change the JSON format, change `app/app.js` in the
same commit.

## Data contract (`app/data/`)

`index.json`: one entry per perturbation, used for the search box

```json
{
  "dataset": "K562_essential",
  "source": "Replogle et al. 2022, normalized pseudobulk (…_normalized_bulk_01.h5ad)",
  "units": "z-score relative to non-targeting controls",
  "n_genes_measured": 8563,
  "perturbations": [
    {"id": "GATA1", "file": "perturbations/GATA1.json", "n_cells": 180,
     "knockdown_pct": 76, "n_sig": 412}
  ]
}
```

`perturbations/<id>.json`: one perturbation, with its most affected genes

```json
{
  "id": "GATA1",
  "target": "GATA1",
  "target_id": "ENSG00000102145",
  "tss": "P1P2",
  "n_cells": 108,
  "knockdown_pct": 62,
  "knockdown_z": -1.4,
  "n_sig": 5152,
  "paper_de_genes": 5996,
  "genes": [
    {"gene": "LTB", "z": 40.48, "lfc": 5.76, "abs": 12.23, "before": 0.13, "after": 12.36, "sig": true},
    …
  ],
  "top_lists": {"up": ["LTB", …], "down": […]},
  "enrichment": {
    "up": [{"lib": "Hallmark", "term": "Heme Metabolism", "k": 18, "n_list": 100,
            "set_size": 140, "q": 2.1e-09, "genes": ["ALAS2", …]}],
    "down": […]
  }
}
```

- `id`: the gene symbol, or `GENE_P1` / `GENE_P2` when a gene was targeted at more than one TSS
- `knockdown_pct`: from the file's `fold_expr` (target expression vs controls), as % reduction
- `knockdown_z`: the target gene's own value in the matrix
- `n_sig`: genes significantly changed at `--fdr` (default 5%), see below
- `paper_de_genes`: the paper's own DE count (`anderson_darling_counts`), for comparison
- `genes`: the top `--top` (default 30) **significant** genes up and down, ranked by `lfc` (the chart)
  - `lfc`: approximate log2 fold change: `after = mean + z * std`,
    `lfc = log2((after + pc) / (mean + pc))`, using `var['mean']`, `var['std']` from the file
    and pseudocount `--pseudocount` (default 0.1)
  - `z`: the matrix value, change relative to control cell-to-cell noise (hover text)
  - `abs`: approximate absolute change per cell, `after - mean`; `before` / `after`: approximate expression per cell
- `top_lists`: the top `--enrich-top` (default 100) significant genes each way, by `lfc`:
  the lists tested for enrichment and sent to Enrichr by the buttons
- `enrichment`: only present if gene set files were found. Up to 5 gene sets per direction with
  q < `--fdr` and at least 3 genes in common: library, term, overlap `k` of `n_list`, set size
  (within the measured genes), BH q-value, and up to 12 of the shared genes

### Gene set enrichment

`python pipeline/01_download.py --genesets` saves two Enrichr libraries to `raw/genesets/`
(MSigDB Hallmark 2020, GO Biological Process 2023). `03_summarise.py` then tests each `top_lists`
list against every set with 10 to 500 measured genes: hypergeometric test, background = all
measured genes, Benjamini-Hochberg across sets within each list.

### Significance: control-based null

Each row is a pseudobulk (mean over n cells). For every gene, the non-targeting
control pseudobulks give a single-cell variance estimate
`s2_g = sum_c n_c (x_cg - mu_g)^2 / (C - 1)`. A perturbation with n cells gets
`t = (x - mu_g) / sqrt(s2_g / n)`, normal p-values, then Benjamini-Hochberg across
genes within that perturbation. `03_summarise.py` prints a leave-one-out
calibration on the controls: the number of genes it calls in unperturbed
pseudobulks, which should be close to zero.

## Running it

```bash
# 1. environment (once)
conda env create -f pipeline/environment.yml
conda activate perturbseq-explorer

# 2. download the essential-screen summary file (~300 MB-1 GB) into raw/
python pipeline/01_download.py --dataset K562_essential

# 3. look inside it before trusting anything
python pipeline/02_inspect.py raw/K562_essential_normalized_bulk_01.h5ad

# 4. gene sets for the enrichment panel (once; two small text files into raw/genesets/)
python pipeline/01_download.py --genesets

# 5. build the JSON for the app
python pipeline/03_summarise.py raw/K562_essential_normalized_bulk_01.h5ad \
    --dataset K562_essential --top 30 --fdr 0.05

# the first run caches the matrix next to the .h5ad (*.cache.pkl), so
# re-running with different options takes seconds; --no-cache forces a re-read

# 6. view locally (browsers block fetch() from file://, so use a tiny server)
# run this from the perturbseq-explorer folder; Ctrl+C to stop
python -m http.server 8000 --directory app
# open http://localhost:8000
# "Address already in use" = an old server is still running: use 8001, or
#   lsof -ti :8000 | xargs kill
```

To try the app without downloading anything, `python pipeline/make_demo_data.py`
writes a fake dataset in the same format.

## Deploying

Push to GitHub, then in the repo go to **Settings → Pages → Source: GitHub
Actions**. The workflow in `.github/workflows/pages.yml` publishes `app/` to
`https://annacuomo.github.io/perturbseq-explorer/`.

## Data source and licence

Code: MIT (see `LICENSE`).

Data: the processed datasets are by Replogle & Weissman, Figshare+,
[doi:10.25452/figshare.plus.20029387](https://doi.org/10.25452/figshare.plus.20029387),
licensed [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The JSON in `app/data/` is
derived from them (summarised and reformatted). Please cite
[Replogle et al. 2022, *Cell*](https://doi.org/10.1016/j.cell.2022.05.013).
