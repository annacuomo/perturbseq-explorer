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
    args = p.parse_args()
    download(args.dataset, args.out)
