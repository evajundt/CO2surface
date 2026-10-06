"""
Cut a large SOCAT download down to the NW Gulf box so it is small enough to
put on GitHub (GitHub rejects files over 100 MB).

    python scripts/subset_socat.py <big_socat_file> data/raw/socat_nwgom.tsv

* Reads the file in chunks, so it works on multi-GB downloads.
* Keeps every column and writes text straight through: do NOT open the
  result in Excel before uploading, or Expocodes like 316420090512 get turned
  into '3.16E+11' and separate cruises merge into one.
* A 0.5-degree buffer is kept around DOMAIN so neighbours just outside the
  edge are still available for the spatial weights later.
"""

import sys

import pandas as pd

from loaders import _find_col, _read_header_index, _SOCAT_CANDIDATES
from nwgom_coverage_check import DOMAIN

BUFFER_DEG = 0.5


def subset(src, dst, chunksize=200_000):
    header_idx, sep = _read_header_index(src)
    reader = pd.read_csv(src, sep=sep, skiprows=header_idx, dtype=str,
                         chunksize=chunksize, low_memory=False)
    n_in = n_out = 0
    first = True
    for chunk in reader:
        chunk = chunk.dropna(how="all")
        lon_c = _find_col(chunk.columns, _SOCAT_CANDIDATES["lon"])
        lat_c = _find_col(chunk.columns, _SOCAT_CANDIDATES["lat"])
        lon = pd.to_numeric(chunk[lon_c], errors="coerce")
        lon = lon.where(lon <= 180, lon - 360)
        lat = pd.to_numeric(chunk[lat_c], errors="coerce")
        keep = (lon.between(DOMAIN["lon_min"] - BUFFER_DEG, DOMAIN["lon_max"] + BUFFER_DEG) &
                lat.between(DOMAIN["lat_min"] - BUFFER_DEG, DOMAIN["lat_max"] + BUFFER_DEG))
        n_in += len(chunk)
        n_out += int(keep.sum())
        chunk[keep].to_csv(dst, sep="\t", index=False, mode="w" if first else "a",
                           header=first)
        first = False
    print(f"kept {n_out:,} of {n_in:,} rows -> {dst}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    subset(sys.argv[1], sys.argv[2])
