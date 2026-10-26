"""Write a small FAKE dataset into app/data/ so the app can be tried without
downloading anything. Row names mimic the real Replogle format; the numbers
are random. Re-run 03_summarise.py on the real file to replace it.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

here = Path(__file__).parent
spec = importlib.util.spec_from_file_location("summarise", here / "03_summarise.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

rng = np.random.default_rng(0)
symbols = ("GATA1 TAL1 KLF1 HBG1 HBG2 HBZ HBA1 MYC MAX RPL3 RPL7 RPS6 RPS19 EIF3A EIF4G1 "
           "POLR2A POLR2B MED12 MED14 SF3B1 SRSF1 U2AF1 NDUFA4 COX4I1 ATP5F1A SDHB "
           "TP53 CDKN1A MDM2 BAX GADD45A ATF4 DDIT3 XBP1 HSPA5 CHAC1 TRIB3 ASNS "
           "SLC7A11 SESN2 CCNB1 CDK1 PLK1 AURKB MKI67 TOP2A E2F1 MCM2 PCNA").split()
symbols += [f"GENE{i}" for i in range(400)]
gene_ids = [f"ENSG{i:011d}" for i in range(len(symbols))]
genes = pd.DataFrame({
    "symbol": symbols,
    "mean": np.exp(rng.normal(0, 1.5, len(symbols))),   # fake per-cell expression levels
    "std": 0.0,
}, index=gene_ids)
genes["std"] = np.sqrt(genes["mean"]) * rng.uniform(0.6, 1.5, len(symbols))

rows = [f"{1000 + i}_{s}_P1P2_{gene_ids[i]}" for i, s in enumerate(symbols[:40])]
rows.append("9999_MYC_P2_" + gene_ids[symbols.index("MYC")])  # a second TSS
n_targeting = len(rows)
rows += [f"{20000 + c}_non-targeting_non-targeting_non-targeting" for c in range(60)]

n_cells = rng.integers(40, 400, len(rows))
cell_sd = rng.uniform(0.6, 1.4, len(symbols))  # per-gene single-cell noise
X = rng.normal(0, 1, size=(len(rows), len(symbols))) * cell_sd / np.sqrt(n_cells)[:, None]
for i, row in enumerate(rows[:n_targeting]):
    _, _, gid = m.parse_row_name(row)
    X[i, gene_ids.index(gid)] = rng.uniform(-2, -0.8)          # knockdown of the target
    hits = rng.choice(len(symbols), size=rng.integers(5, 60), replace=False)
    X[i, hits] += rng.normal(0, 0.6, size=hits.size)          # downstream effects

expr = pd.DataFrame(X, index=rows, columns=gene_ids)
meta = pd.DataFrame([m.parse_row_name(r) for r in rows], index=rows, columns=["target", "tss", "target_id"])
meta["n_cells"] = n_cells
meta["fold_expr"] = rng.uniform(0.05, 0.4, len(rows))

index_entries, records, info = m.summarise(expr, meta, genes, top=25)
print("  null:", info)
m.write_json(index_entries, records, here.parent / "app" / "data", {
    "dataset": "DEMO (random numbers)",
    "source": "make_demo_data.py: fake data for testing the app",
    "units": "fake",
    "n_genes_measured": len(symbols),
    "fdr": 0.05,
    "pseudocount": 0.1,
    "has_fold_change": info["has_fold_change"],
    "null": info,
})
