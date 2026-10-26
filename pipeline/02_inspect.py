"""Print what is inside an .h5ad file before we build anything on it.

Prints how long each step takes, so we can see what is slow.
"""

import sys
import time

t0 = time.perf_counter()


def tick(msg):
    print(f"  [{time.perf_counter() - t0:6.1f}s] {msg}")


import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402

tick("imported anndata")

path = sys.argv[1]
adata = ad.read_h5ad(path)  # full read: these bulk files fit comfortably in RAM
tick("read file")

print(f"\n{path}\n  shape: {adata.n_obs} rows (obs) x {adata.n_vars} columns (var)")
print(f"  matrix type: {type(adata.X).__name__}, dtype {adata.X.dtype}\n")

print("obs (rows) - first 5:")
print(adata.obs.head().to_string(), "\n")
print("obs columns:", list(adata.obs.columns), "\n")
print("var (columns) - first 5:")
print(adata.var.head().to_string(), "\n")

X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
a = np.abs(X[np.isfinite(X)])
print("distribution of |value| over the whole matrix:")
for q in (0.5, 0.9, 0.99, 0.999, 0.9999):
    print(f"  {q:>7.2%} quantile: {np.quantile(a, q):.3f}")
for thr in (0.1, 0.25, 0.5, 1.0):
    print(f"  fraction |x| > {thr}: {(a > thr).mean():.4%}")

# rows that look like controls
names = adata.obs_names.to_series()
ctrl = names[names.str.contains("non-targeting|control", case=False)]
print(f"\nrows that look like controls: {len(ctrl)}", list(ctrl[:3]))
if "core_control" in adata.obs:
    print("core_control == True:", int(adata.obs["core_control"].sum()))
tick("done")
