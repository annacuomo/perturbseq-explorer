"""Turn the normalized-bulk .h5ad into the small JSON files the app reads.

Structure:
  load_matrix()   .h5ad -> DataFrames. Cached next to the .h5ad after the first
                  run, so re-runs (e.g. trying thresholds) skip the slow read.
  summarise()     DataFrames -> index + one record per perturbation (pure logic)
  write_json()    records -> app/data/

Row names in the Replogle bulk files look like
    10023_ZC3H18_P1P2_ENSG00000158545
    <id>_<gene symbol>_<TSS>_<Ensembl gene id>

Three ways to rank genes, all stored so the app can switch between them:
  z    the matrix value: change relative to the gene's noise in control cells
       ("how surprising"). Favours quiet genes switching on.
  lfc  approximate log2 fold change ("how big"). Uses the per-gene mean and std
       stored in the file's var table to undo the z-scaling:
           after ~= mean + z * std,   lfc = log2((after + pc) / (mean + pc))
       Relative, so near-silent genes switching on still rank high.
  abs  approximate absolute change per cell: after - mean. Favours abundant genes.
  lfc and abs are APPROXIMATE: they assume var['mean'], var['std'] describe
  control cells on one scale.
"""

import argparse
import json
import re
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd

T0 = time.perf_counter()
CONTROL_PATTERN = re.compile(r"non[-_]?targeting|^control$|^NTC", re.IGNORECASE)


def tick(msg):
    print(f"  [{time.perf_counter() - T0:6.1f}s] {msg}")


# ---------------------------------------------------------------- loading
def parse_row_name(name: str):
    """'10023_ZC3H18_P1P2_ENSG00000158545' -> ('ZC3H18', 'P1P2', 'ENSG00000158545').

    Symbol = everything between the first and the last two fields, so symbols
    that contain underscores survive.
    """
    parts = name.split("_")
    if len(parts) >= 4:
        return "_".join(parts[1:-2]), parts[-2], parts[-1]
    return name, "", ""


def load_matrix(path: Path, cells_col: str = "num_cells_filtered", use_cache: bool = True):
    """Return (expr, meta, genes).

    expr:  DataFrame, index = row name, columns = Ensembl gene ids (float32)
    meta:  DataFrame per row: target, tss, target_id, n_cells, fold_expr, paper_de_genes
    genes: DataFrame per gene id: symbol, mean, std
    """
    cache = path.with_suffix(".v3.cache.pkl")  # bump the version when the cached shape changes
    if use_cache and cache.exists() and cache.stat().st_mtime > path.stat().st_mtime:
        expr, meta, genes = pd.read_pickle(cache)
        tick(f"loaded cache {cache.name}")
        return expr, meta, genes

    import anndata as ad  # only needed for the first read
    tick("imported anndata")
    adata = ad.read_h5ad(path)
    tick(f"read {path.name}")

    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    expr = pd.DataFrame(X.astype(np.float32), index=adata.obs_names, columns=adata.var_names)

    var = adata.var
    genes = pd.DataFrame(index=adata.var_names)
    genes["symbol"] = (var["gene_name"] if "gene_name" in var else adata.var_names.to_series()).astype(str)
    genes["mean"] = var["mean"] if "mean" in var else np.nan
    genes["std"] = var["std"] if "std" in var else np.nan

    parsed = [parse_row_name(n) for n in adata.obs_names]
    meta = pd.DataFrame(parsed, index=adata.obs_names, columns=["target", "tss", "target_id"])
    meta["n_cells"] = adata.obs[cells_col] if cells_col in adata.obs else np.nan
    meta["fold_expr"] = adata.obs["fold_expr"] if "fold_expr" in adata.obs else np.nan
    # the paper's own count of DE genes (Anderson-Darling test), for comparison
    meta["paper_de_genes"] = (adata.obs["anderson_darling_counts"]
                              if "anderson_darling_counts" in adata.obs else np.nan)

    pd.to_pickle((expr, meta, genes), cache)
    tick(f"cached to {cache.name}")
    return expr, meta, genes


# ---------------------------------------------------------------- control-based null
#
# Each row is a pseudobulk: the mean over n cells of per-cell expression,
# z-normalised against control cells. For gene g, if a single cell has
# variance s2_g, a pseudobulk of n cells has variance s2_g / n.
#
# The non-targeting rows are pseudobulks of unperturbed cells, so they give an
# empirical estimate of s2_g for every gene:
#     s2_g = sum_c n_c (x_cg - mu_g)^2 / (C - 1)
# A perturbation p with n_p cells then gets
#     t_pg = (x_pg - mu_g) / sqrt(s2_g / n_p)
# two-sided normal p-values, and Benjamini-Hochberg within each perturbation.


def control_null(x_ctrl: np.ndarray, n_ctrl: np.ndarray):
    """Per-gene control mean and single-cell variance from control pseudobulks."""
    ok = np.isfinite(n_ctrl) & (n_ctrl > 0)
    x, n = x_ctrl[ok], n_ctrl[ok][:, None]
    mu = np.nanmean(x, axis=0)
    s2 = np.nansum(n * (x - mu) ** 2, axis=0) / (len(x) - 1)
    floor = np.nanquantile(s2[s2 > 0], 0.01)  # guard against ~zero-variance genes
    return mu, np.maximum(s2, floor)


def bh(P: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values, row by row on a matrix (NaN = untested)."""
    P = np.atleast_2d(P)
    m = np.isfinite(P).sum(axis=1, keepdims=True)           # tests per row
    order = np.argsort(np.where(np.isfinite(P), P, np.inf), axis=1)
    p_sorted = np.take_along_axis(P, order, axis=1)
    rank = np.arange(1, P.shape[1] + 1)
    q_sorted = p_sorted * m / rank
    q_sorted = np.where(np.isfinite(q_sorted), q_sorted, np.inf)
    q_sorted = np.minimum(1, np.minimum.accumulate(q_sorted[:, ::-1], axis=1)[:, ::-1])
    Q = np.empty_like(q_sorted)
    np.put_along_axis(Q, order, q_sorted, axis=1)
    return np.where(np.isfinite(P), Q, np.nan)


def test_matrix(X: np.ndarray, n: np.ndarray, mu: np.ndarray, s2: np.ndarray):
    """t-statistics and BH q-values for every row of X at once."""
    from scipy.special import ndtr  # standard normal CDF; scipy comes with anndata

    T = (X - mu) / np.sqrt(s2 / n[:, None])
    P = 2 * ndtr(-np.abs(T))
    return T, bh(P)


def calibrate(x_ctrl, n_ctrl, fdr):
    """Leave-one-out on the controls: how many genes does the test call
    'changed' in an unperturbed pseudobulk? Should be close to zero."""
    hits = []
    for c in np.flatnonzero(np.isfinite(n_ctrl)):
        keep = np.arange(len(x_ctrl)) != c
        mu, s2 = control_null(x_ctrl[keep], n_ctrl[keep])
        _, q = test_matrix(x_ctrl[c:c + 1], n_ctrl[c:c + 1], mu, s2)
        hits.append(int((q < fdr).sum()))
    return np.array(hits)


# ---------------------------------------------------------------- summarise
def safe_id(label: str) -> str:
    """Make a label safe to use as a filename and URL fragment."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", label)


def r(x, nd):
    return None if x is None or not np.isfinite(x) else round(float(x), nd)


def top_k(M: np.ndarray, k: int, largest: bool) -> np.ndarray:
    """Column indices of the k largest (or smallest) finite values in each row."""
    key = np.where(np.isfinite(M), -M if largest else M, np.inf)
    return np.argsort(key, axis=1)[:, :k]


def summarise(expr, meta, genes, top=30, fdr=0.05, pseudocount=0.1):
    """Per perturbation: knockdown, significance, and a list of the most
    affected genes (top `top` up and down by each of z and lfc, merged).

    Returns (index_entries, records, info).
    """
    # --- null from the controls
    is_control = meta["target"].str.contains(CONTROL_PATTERN).to_numpy()
    x_ctrl = expr.to_numpy(dtype=np.float64)[is_control]
    n_ctrl = meta["n_cells"].to_numpy(dtype=float)[is_control]
    mu, s2 = control_null(x_ctrl, n_ctrl)
    calib = calibrate(x_ctrl, n_ctrl, fdr)
    info = {
        "n_controls": int(is_control.sum()),
        "median_single_cell_sd": float(np.median(np.sqrt(s2))),
        "control_false_hits_median": float(np.median(calib)),
        "control_false_hits_max": int(calib.max()),
    }
    expr, meta = expr.loc[~is_control], meta.loc[~is_control]

    # --- labels: genes targeted at more than one TSS become GENE_P1, GENE_P2
    labels = meta["target"].copy()
    dup = labels.duplicated(keep=False)
    labels[dup] = labels[dup] + "_" + meta.loc[dup, "tss"].where(meta.loc[dup, "tss"] != "", "x")
    labels[labels.duplicated(keep=False)] += "_" + labels.groupby(labels).cumcount().astype(str)

    # --- per-gene info, aligned with the matrix columns
    gene_ids = expr.columns.to_numpy()
    g = genes.reindex(gene_ids)
    symbols = g["symbol"].fillna(pd.Series(gene_ids, index=gene_ids)).to_numpy()
    g_mean = g["mean"].to_numpy(dtype=float)
    g_std = g["std"].to_numpy(dtype=float)
    has_fc = bool(np.isfinite(g_mean).any() and np.isfinite(g_std).any())
    col_of = {gid: j for j, gid in enumerate(gene_ids)}
    col_of_symbol = {s: j for j, s in enumerate(symbols)}

    Z = expr.to_numpy(dtype=np.float64, copy=True)
    rows = meta.to_dict("records")

    # --- the target's own column: record it, then blank it out
    knockdown_z = np.full(len(rows), np.nan)
    for i, m in enumerate(rows):
        j = col_of.get(m["target_id"], col_of_symbol.get(m["target"]))
        if j is not None:
            knockdown_z[i] = Z[i, j]
            Z[i, j] = np.nan

    # --- significance, all perturbations at once
    n = np.array([m["n_cells"] for m in rows], dtype=float)
    _, Q = test_matrix(Z, n, mu, s2)
    SIG = (Q < fdr) & np.isfinite(n)[:, None]
    n_sig_all = SIG.sum(axis=1)

    # --- approximate expression before/after and log2 fold change
    if has_fc:
        AFTER = np.maximum(g_mean + Z * g_std, 0)
        LFC = np.log2((AFTER + pseudocount) / (g_mean + pseudocount))
        ABS = AFTER - g_mean
    else:
        AFTER = LFC = ABS = np.full_like(Z, np.nan)

    # --- candidate genes: top up/down by each ranking, merged
    picks = [top_k(Z, top, False), top_k(Z, top, True)]
    if has_fc:
        picks += [top_k(LFC, top, False), top_k(LFC, top, True),
                  top_k(ABS, top, False), top_k(ABS, top, True)]
    picks = np.concatenate(picks, axis=1)

    records, index_entries = [], []
    for i, (m, label) in enumerate(zip(rows, labels)):
        cols = [k for k in dict.fromkeys(picks[i]) if np.isfinite(Z[i, k])]
        gene_list = [{
            "gene": str(symbols[k]),
            "z": r(Z[i, k], 3),
            "lfc": r(LFC[i, k], 3),
            "abs": r(ABS[i, k], 2),
            "before": r(g_mean[k], 2),
            "after": r(AFTER[i, k], 2),
            "sig": bool(SIG[i, k]),
        } for k in cols]
        gene_list.sort(key=lambda d: d["z"])

        fold = m["fold_expr"]
        knockdown_pct = None if pd.isna(fold) else round(100 * (1 - float(fold)))
        n_cells = None if pd.isna(m["n_cells"]) else int(m["n_cells"])
        n_sig = int(n_sig_all[i]) if n_cells else None
        paper = m.get("paper_de_genes")
        paper = None if paper is None or pd.isna(paper) else int(paper)

        records.append({
            "id": label, "target": m["target"], "target_id": m["target_id"], "tss": m["tss"],
            "n_cells": n_cells, "knockdown_pct": knockdown_pct, "knockdown_z": r(knockdown_z[i], 3),
            "n_sig": n_sig, "paper_de_genes": paper, "genes": gene_list,
        })
        index_entries.append({
            "id": label, "file": f"perturbations/{safe_id(label)}.json",
            "n_cells": n_cells, "knockdown_pct": knockdown_pct, "n_sig": n_sig,
        })

    index_entries.sort(key=lambda e: e["id"])
    info["has_fold_change"] = has_fc
    return index_entries, records, info


# ---------------------------------------------------------------- output
def write_json(index_entries, records, out_dir: Path, header: dict):
    pert_dir = out_dir / "perturbations"
    if pert_dir.exists():
        shutil.rmtree(pert_dir)  # avoid stale files from a previous run
    pert_dir.mkdir(parents=True)

    for rec in records:
        with open(pert_dir / f"{safe_id(rec['id'])}.json", "w") as f:
            json.dump(rec, f, separators=(",", ":"))
    with open(out_dir / "index.json", "w") as f:
        json.dump({**header, "perturbations": index_entries}, f, separators=(",", ":"))

    size_mb = sum(p.stat().st_size for p in out_dir.rglob("*.json")) / 1e6
    tick(f"wrote {len(records)} perturbations to {out_dir} ({size_mb:.1f} MB total)")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("h5ad", type=Path)
    p.add_argument("--dataset", default="K562_essential")
    p.add_argument("--top", type=int, default=30, help="genes kept per direction, per ranking")
    p.add_argument("--fdr", type=float, default=0.05, help="BH FDR for calling a gene changed")
    p.add_argument("--pseudocount", type=float, default=0.1, help="added before taking log2 fold change")
    p.add_argument("--cells-col", default="num_cells_filtered")
    p.add_argument("--no-cache", action="store_true", help="re-read the .h5ad even if a cache exists")
    p.add_argument("--out", type=Path, default=Path(__file__).resolve().parent.parent / "app" / "data")
    args = p.parse_args()

    expr, meta, genes = load_matrix(args.h5ad, args.cells_col, use_cache=not args.no_cache)
    tick(f"{expr.shape[0]} rows x {expr.shape[1]} genes")
    index_entries, records, info = summarise(expr, meta, genes, top=args.top,
                                             fdr=args.fdr, pseudocount=args.pseudocount)
    tick("summarised")

    n_sig = np.array([rec["n_sig"] for rec in records if rec["n_sig"] is not None])
    print(f"\n  control-based null from {info['n_controls']} non-targeting pseudobulks")
    print(f"    median single-cell SD across genes: {info['median_single_cell_sd']:.3f}")
    print(f"    calibration (leave-one-out on controls): genes called at FDR {args.fdr}:"
          f" median {info['control_false_hits_median']:.0f}, max {info['control_false_hits_max']}")
    print(f"    perturbations: genes called, median {np.median(n_sig):.0f},"
          f" 90th pct {np.quantile(n_sig, 0.9):.0f}, max {n_sig.max()}")
    pairs = np.array([(rec["n_sig"], rec["paper_de_genes"]) for rec in records
                      if rec["n_sig"] is not None and rec["paper_de_genes"] is not None], dtype=float)
    if len(pairs):
        rho = pd.Series(pairs[:, 0]).corr(pd.Series(pairs[:, 1]), method="spearman")
        print(f"    paper's Anderson-Darling DE counts: median {np.median(pairs[:, 1]):.0f},"
              f" 90th pct {np.quantile(pairs[:, 1], 0.9):.0f}, max {pairs[:, 1].max():.0f}")
        print(f"    Spearman correlation, ours vs paper's: {rho:.2f}")
    print(f"    fold change available: {info['has_fold_change']}\n")

    write_json(index_entries, records, args.out, {
        "dataset": args.dataset,
        "source": f"Replogle et al. 2022, normalized pseudobulk ({args.h5ad.name})",
        "units": "z-normalised expression relative to non-targeting controls",
        "n_genes_measured": int(expr.shape[1]),
        "fdr": args.fdr,
        "pseudocount": args.pseudocount,
        "has_fold_change": info["has_fold_change"],
        "null": info,
    })
