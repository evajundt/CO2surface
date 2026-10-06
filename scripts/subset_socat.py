"""
Cut a large SOCAT download down to the NW Gulf box so it is small enough to
put on GitHub. Uploading through the GitHub website is limited to 25 MB per
file (100 MB via git push), so write a gzip-compressed .tsv.gz:

    python scripts/subset_socat.py <big_socat_file> data/raw/socat_nwgom.tsv.gz
    python scripts/subset_socat.py --info <big_socat_file>   # print layout only

* Works on text (.tsv/.csv) and NetCDF (.nc) downloads; output is always
  tab-separated text, gzip-compressed when the output name ends in .gz
  (typically 5-10x smaller; the loaders read .gz directly).
* Text is read in chunks; for NetCDF only latitude/longitude are read in
  full, then just the rows inside the box. Either way multi-GB files work.
* Keeps every column and writes text straight through: do NOT open the
  result in Excel before uploading, or Expocodes like 316420090512 get turned
  into '3.16E+11' and separate cruises merge into one.
* A 0.5-degree buffer is kept around DOMAIN so neighbours just outside the
  edge are still available for the spatial weights later.
"""

import sys

import numpy as np
import pandas as pd

from loaders import (_SOCAT_CANDIDATES, _find_col, _is_netcdf, _netcdf_layout,
                     _netcdf_to_frame, _read_header_index)
from nwgom_coverage_check import DOMAIN

BUFFER_DEG = 0.5


def _in_box(lon, lat):
    lon = lon.where(lon <= 180, lon - 360)
    return (lon.between(DOMAIN["lon_min"] - BUFFER_DEG, DOMAIN["lon_max"] + BUFFER_DEG) &
            lat.between(DOMAIN["lat_min"] - BUFFER_DEG, DOMAIN["lat_max"] + BUFFER_DEG))


def subset_netcdf(src, dst):
    import xarray as xr
    with xr.open_dataset(src) as ds:
        obs_dim, _, _ = _netcdf_layout(ds)
        lon_v = _find_col(list(ds.variables), _SOCAT_CANDIDATES["lon"])
        lat_v = _find_col(list(ds.variables), _SOCAT_CANDIDATES["lat"])
        lon = pd.Series(ds[lon_v].values.astype(float))
        lat = pd.Series(ds[lat_v].values.astype(float))
        idx = np.flatnonzero(_in_box(lon, lat).to_numpy())
        print(f"kept {len(idx):,} of {ds.sizes[obs_dim]:,} measurements")
        _netcdf_to_frame(ds, obs_index=idx).to_csv(dst, sep="\t", index=False)
    print(f"-> {dst}")
    _report_size(dst)


def subset(src, dst, chunksize=200_000):
    if _is_netcdf(src):
        return subset_netcdf(src, dst)
    header_idx, sep = _read_header_index(src)
    reader = pd.read_csv(src, sep=sep, skiprows=header_idx, dtype=str,
                         chunksize=chunksize, low_memory=False)
    n_in = n_out = 0
    first = True
    for chunk in reader:
        chunk = chunk.dropna(how="all")
        lon = pd.to_numeric(chunk[_find_col(chunk.columns, _SOCAT_CANDIDATES["lon"])],
                            errors="coerce")
        lat = pd.to_numeric(chunk[_find_col(chunk.columns, _SOCAT_CANDIDATES["lat"])],
                            errors="coerce")
        keep = _in_box(lon, lat)
        n_in += len(chunk)
        n_out += int(keep.sum())
        chunk[keep].to_csv(dst, sep="\t", index=False, mode="w" if first else "a",
                           header=first)
        first = False
    print(f"kept {n_out:,} of {n_in:,} rows -> {dst}")
    _report_size(dst)


def _report_size(dst):
    import os
    mb = os.path.getsize(dst) / 1e6
    note = ("OK for website upload" if mb < 25 else
            "too big for website upload (25 MB); use a .gz name or git push" if mb < 100 else
            "too big for GitHub (100 MB); tell Claude")
    print(f"file size: {mb:.1f} MB - {note}")


def info(src):
    """Print the file's variables and dimensions (small enough to paste into chat)."""
    if _is_netcdf(src):
        import xarray as xr
        with xr.open_dataset(src) as ds:
            print(ds)
            print("\nlayout (obs_dim, trajectory_dim, rowSize):", _netcdf_layout(ds))
    else:
        header_idx, sep = _read_header_index(src)
        print(pd.read_csv(src, sep=sep, skiprows=header_idx, nrows=5).T)


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--info":
        info(sys.argv[2])
    elif len(sys.argv) == 3:
        subset(sys.argv[1], sys.argv[2])
    else:
        sys.exit(__doc__)
