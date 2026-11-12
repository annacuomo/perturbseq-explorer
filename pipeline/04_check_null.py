"""Check the control-based null used in 03_summarise.py.

Question: is the test too liberal for pseudobulks made of few cells?

Evidence so far (from the published JSON): our DE-gene counts agree with the
paper's overall (Spearman 0.96), but perturbations with < 50 cells get far more
calls from us (median 55) than from the paper (median 10.5). One non-targeting
control also gets 41 false calls in the leave-one-out calibration.

This script, run on the same .h5ad (it reuses the cache from 03_summarise.py):
  1. Leave-one-out on every control: cell count, number of false calls, and the
     called genes with their control mean expression. If false calls pile up in
     controls with few cells, the 1/n variance scaling is the problem.
  2. For the worst control, the distribution of its t-statistics against N(0, 1).
  3. Among perturbations, calls vs cell count, split by gene expression level:
     sparse genes have skewed per-cell values, so a mean over few cells is far
     from normal and the normal tail is too thin.

Usage:
  python pipeline/04_check_null.py raw/K562_essential_normalized_bulk_01.h5ad
"""

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

# reuse the exact functions the app's numbers came from
_spec = importlib.util.spec_from_file_location(
    "summarise", Path(__file__).resolve().parent / "03_summarise.py")
S = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S)


def loo_controls(x_ctrl, n_ctrl, fdr):
    """Leave-one-out over controls: false calls and t-statistics per control."""
    out = []
    for c in np.flatnonzero(np.isfinite(n_ctrl)):
        keep = np.arange(len(x_ctrl)) != c
        mu, s2 = S.control_null(x_ctrl[keep], n_ctrl[keep])
        T, Q = S.test_matrix(x_ctrl[c:c + 1], n_ctrl[c:c + 1], mu, s2)
        out.append((c, int((Q[0] < fdr).sum()), T[0], Q[0]))
    return out


def main(path, fdr, cells_col):
    expr, meta, genes = S.load_matrix(path, cells_col)
    is_ctrl = meta["target"].str.contains(S.CONTROL_PATTERN).to_numpy()
    X = expr.to_numpy(dtype=np.float64)
    x_ctrl, n_ctrl = X[is_ctrl], meta["n_cells"].to_numpy(dtype=float)[is_ctrl]
    ctrl_names = expr.index[is_ctrl]
    g_mean = genes.reindex(expr.columns)["mean"].to_numpy(dtype=float)
    symbols = genes.reindex(expr.columns)["symbol"].to_numpy()

    # 1. leave-one-out over controls
    res = loo_controls(x_ctrl, n_ctrl, fdr)
    tab = pd.DataFrame({
        "control": [ctrl_names[c] for c, *_ in res],
        "n_cells": [n_ctrl[c] for c, *_ in res],
        "false_calls": [k for _, k, *_ in res],
    }).sort_values("false_calls", ascending=False)
    print("\n1. Leave-one-out on controls (worst 10)")
    print(tab.head(10).to_string(index=False))
    rho = tab["n_cells"].corr(tab["false_calls"], method="spearman")
    print(f"   controls: n_cells median {tab.n_cells.median():.0f}, "
          f"range {tab.n_cells.min():.0f}-{tab.n_cells.max():.0f}")
    print(f"   Spearman(n_cells, false calls) = {rho:.2f}")
    for lo, hi in [(0, 50), (50, 100), (100, 200), (200, 1e9)]:
        m = (tab.n_cells >= lo) & (tab.n_cells < hi)
        if m.any():
            print(f"   n_cells {lo:>3}-{hi if hi < 1e9 else 'max'}: {m.sum():>3} controls, "
                  f"mean false calls {tab.false_calls[m].mean():.1f}")

    # 2. the worst control in detail
    c, k, T, Q = max(res, key=lambda r: r[1])
    print(f"\n2. Worst control: {ctrl_names[c]} ({n_ctrl[c]:.0f} cells, {k} false calls)")
    called = np.flatnonzero(Q < fdr)
    detail = pd.DataFrame({"gene": symbols[called], "t": T[called].round(2),
                           "control_mean_expr": g_mean[called].round(3)}).sort_values("t")
    print(detail.to_string(index=False))
    t = T[np.isfinite(T)]
    print(f"   t-statistics: sd {t.std():.2f} (null: 1.00), "
          f"|t| > 4: {(np.abs(t) > 4).mean():.5f} (null: {2 * 3.17e-5:.5f})")
    print(f"   median control mean of called genes {np.median(g_mean[called]):.3f} "
          f"vs all genes {np.nanmedian(g_mean):.3f}")

    # 3. perturbations: calls vs cell count, by expression level of the gene
    Xp, n_p = X[~is_ctrl], meta["n_cells"].to_numpy(dtype=float)[~is_ctrl]
    mu, s2 = S.control_null(x_ctrl, n_ctrl)
    _, Qp = S.test_matrix(Xp, n_p, mu, s2)
    sparse = g_mean < np.nanquantile(g_mean, 0.25)
    print("\n3. Perturbations: share of tested genes called, by cell count and gene level")
    print("   (if the null is fine, the sparse/abundant ratio should not jump at low n)")
    for lo, hi in [(0, 50), (50, 100), (100, 200), (200, 1e9)]:
        m = (n_p >= lo) & (n_p < hi)
        if m.any():
            a = (Qp[m][:, sparse] < fdr).mean()
            b = (Qp[m][:, ~sparse] < fdr).mean()
            print(f"   n_cells {lo:>3}-{hi if hi < 1e9 else 'max'}: {m.sum():>4} perts, "
                  f"sparse genes {a:.4f}, other genes {b:.4f}, ratio {a / max(b, 1e-9):.2f}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("h5ad", type=Path)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--cells-col", default="num_cells_filtered")
    a = p.parse_args()
    main(a.h5ad, a.fdr, a.cells_col)
