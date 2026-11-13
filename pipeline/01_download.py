"""Download the per-perturbation summary files from Replogle et al. 2022.

The Figshare+ release (doi:10.25452/figshare.plus.20029387) has, for each of the
three screens, four AnnData files. The one we want for the app is the
*normalized bulk* file: one row per perturbation, one column per measured gene,
values are pseudobulk expression z-scored against non-targeting control cells.

Screens:
  K562_essential  ~2,000 essential genes, K562, day 6   <- start here (smallest)
  rpe1            ~2,400 essential genes, RPE1, day 7
  K562_gwps       ~9,900 expressed genes, K562, day 8   (genome-wide)

The file IDs below were collected from public code that uses these files.
If a download fails or the size looks wrong, open the Figshare page and
check the file list there.
"""

import argparse
from pathlib import Path

import requests
from tqdm import tqdm

FIGSHARE_FILES = {
    "K562_essential": ("K562_essential_normalized_bulk_01.h5ad", 35780870),
    "K562_gwps": ("K562_gwps_normalized_bulk_01.h5ad", 35773217),
    "rpe1": ("rpe1_normalized_bulk_01.h5ad", 35775512),
}


# Gene set libraries for enrichment, as plain-text files from Enrichr's library
# download endpoint (one set per line: TERM<tab><tab>GENE<tab>GENE...).
# Hallmark: 50 well-defined processes (MSigDB). GO BP: thousands of finer terms.
GENESET_LIBRARIES = ["MSigDB_Hallmark_2020", "GO_Biological_Process_2023"]
ENRICHR_LIBRARY_URL = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={}"


def download_genesets(out_dir: Path) -> None:
    out_dir = out_dir / "genesets"
    out_dir.mkdir(parents=True, exist_ok=True)
    for name in GENESET_LIBRARIES:
        dest = out_dir / f"{name}.txt"
        if dest.exists():
            print(f"{dest} already exists, skipping")
            continue
        r = requests.get(ENRICHR_LIBRARY_URL.format(name), timeout=120)
        r.raise_for_status()
        dest.write_text(r.text)
        print(f"Saved {name}: {len(r.text.splitlines())} gene sets -> {dest}")


def download(dataset: str, out_dir: Path) -> Path:
    filename, file_id = FIGSHARE_FILES[dataset]
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / filename
    if dest.exists():
        print(f"{dest} already exists, skipping")
        return dest

    url = f"https://ndownloader.figshare.com/files/{file_id}"
    print(f"Downloading {filename} from {url}")
    tmp = dest.with_suffix(".part")
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with open(tmp, "wb") as f, tqdm(total=total, unit="B", unit_scale=True) as bar:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
                bar.update(len(chunk))
    tmp.rename(dest)
    print(f"Saved to {dest}")
    return dest


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", choices=FIGSHARE_FILES, default="K562_essential")
    p.add_argument("--out", type=Path, default=Path("raw"))
    p.add_argument("--genesets", action="store_true",
                   help="download the gene set libraries for enrichment instead of a dataset")
    args = p.parse_args()
    if args.genesets:
        download_genesets(args.out)
    else:
        download(args.dataset, args.out)
