"""
Loaders that turn each raw source into one common surface-point table:

    Date | Latitude | Longitude | Temp | Sal | pCO2 | Source | Station

All sources end up with the same column names, so the coverage check and
the later Moran's I / LISA steps can be run on any one of them or on all
of them stacked together.
"""

import gzip
import re

import numpy as np
import pandas as pd

COMMON_COLS = ["Date", "Latitude", "Longitude", "Temp", "Sal", "pCO2", "Source", "Station"]
FILL_VALUES = [-999, -999.0, 999, 999.0, -9999, -99]


# ---------------------------------------------------------------------------
# Ship discrete samples
# ---------------------------------------------------------------------------
def load_ship(path, depth="S"):
    """
    Shipboard bottle data (shipto2023.xlsx), surface samples only by default.

    Fixes applied:
      * the 2022-12-03 WFG rows have Longitude = +93.85 (sign error), so any
        positive longitude is flipped to west
      * -999 fill values become NaN
    """
    raw = pd.read_excel(path)
    raw = raw.replace(FILL_VALUES, np.nan)
    if depth is not None:
        raw = raw[raw["Depth"].astype(str).str.strip() == depth]

    n_flipped = int((raw["Longitude"] > 0).sum())
    if n_flipped:
        print(f"[ship] flipped sign on {n_flipped} positive longitude(s)")
    lon = -raw["Longitude"].abs()

    out = pd.DataFrame({
        "Date": pd.to_datetime(raw["Date"]),
        "Latitude": raw["Latitude"],
        "Longitude": lon,
        "Temp": raw["CTD temp"],
        "Sal": raw["CTD Sal"],
        "pCO2": raw["pCO2"],
        "Source": "Ship",
        "Station": raw["Station"].astype(str).str.strip(),
    })
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# HOBO / SBE reef-cap temperature loggers
# ---------------------------------------------------------------------------
# Approximate East/West FGB reef-cap positions (same as ship EFG / WFG stations).
LOGGER_SITES = {
    "East": {"Latitude": 27.91, "Longitude": -93.60, "date_col": "Date",
             "temp_col": "Etemp_reef cap SBE"},
    "West": {"Latitude": 27.87, "Longitude": -93.82, "date_col": "Date.1",
             "temp_col": "Wtemp_reef cap SBE"},
}


def load_reef_loggers(path):
    """
    Daily-average reef-cap temperature (Temp_E_W_1989-2024_1sheet.xlsx).

    The East and West columns each have their OWN date column, and the two do
    not line up row-by-row (e.g. row 2922 is 2002-01-01 for East but
    2003-01-01 for West). Each pair is read separately so they are never
    matched up by row.

    Note: these are reef-cap loggers at roughly 20 m depth, not surface
    sensors. Treat them as a separate layer, not as surface temperature.
    """
    raw = pd.read_excel(path)
    parts = []
    for site, meta in LOGGER_SITES.items():
        d = raw[[meta["date_col"], meta["temp_col"]]].dropna()
        # Fix: West has 2023-01-01..10 entered twice (identical values);
        # keep one value per day.
        d = d.groupby(meta["date_col"], as_index=False)[meta["temp_col"]].mean()
        parts.append(pd.DataFrame({
            "Date": pd.to_datetime(d[meta["date_col"]]),
            "Latitude": meta["Latitude"],
            "Longitude": meta["Longitude"],
            "Temp": d[meta["temp_col"]].astype(float),
            "Sal": np.nan,
            "pCO2": np.nan,
            "Source": "ReefLogger",
            "Station": site,
        }))
    return pd.concat(parts, ignore_index=True)


# ---------------------------------------------------------------------------
# SOCAT
# ---------------------------------------------------------------------------
# Column names differ between the SOCAT synthesis files, the Data Set Viewer
# export and the ERDDAP export. The first match in each list is used, and
# matching ignores upper/lower case.
_SOCAT_CANDIDATES = {
    "lon":  [r"^longitude", r"^lon"],
    "lat":  [r"^latitude", r"^lat"],
    "sst":  [r"^SST", r"^temp$", r"^temp\b", r"^sea_surface_temp"],
    "sal":  [r"^sal", r"^SSS"],
    # (?![_a-z]) stops fCO2rec from matching fCO2rec_flag / fCO2rec_src
    "fco2": [r"^fCO2rec(?![_a-z])", r"^fco2_recommended", r"^fCO2_rec(?![_a-z])"],
    "flag": [r"^fCO2rec_flag", r"^WOCE_CO2_water", r"^fCO2_flag"],
    "qc":   [r"^QC_flag"],
    "expo": [r"^expocode"],
    "time": [r"^time$", r"^datetime$", r"^date$"],
}

# SOCAT writes missing values as -1E+34 (NetCDF fill). Anything below this
# threshold is treated as missing.
SOCAT_FILL_THRESHOLD = -1e30

# SOCAT's own recommendation: cruise QC flags A-D are fit for most uses
# (accuracy better than 5 uatm). E (accuracy 5-10 uatm) is left out by default.
SOCAT_GOOD_QC = ("A", "B", "C", "D")


def _find_col(columns, patterns):
    for pat in patterns:
        for c in columns:
            if re.search(pat, str(c), flags=re.IGNORECASE):
                return c
    return None


def _read_header_index(path, max_lines=10_000):
    """
    Find the header line (SOCAT files can have a long metadata block above
    it) and the delimiter, streaming line by line so large files are never
    loaded whole.
    """
    # .gz files (written by subset_socat.py to stay under GitHub's upload
    # limit) are read transparently.
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        for i, line in enumerate(fh):
            low = line.lower()
            if low.startswith("expocode") or ("latitude" in low and "longitude" in low):
                sep = "\t" if line.count("\t") >= line.count(",") else ","
                return i, sep
            if i >= max_lines:
                break
    raise ValueError(f"No header line with latitude/longitude found in {path}")


NETCDF_SUFFIXES = (".nc", ".nc4", ".netcdf", ".cdf")


def _is_netcdf(path):
    return str(path).lower().endswith(NETCDF_SUFFIXES)


def _decode_text(values):
    """NetCDF char/bytes variables (Expocode, QC flag) -> plain strings."""
    arr = np.asarray(values)
    if arr.dtype.kind == "S":
        return np.char.decode(arr, "utf-8", errors="replace").astype(str)
    if arr.dtype.kind == "O":
        return np.array([v.decode("utf-8", "replace") if isinstance(v, bytes) else str(v)
                         for v in arr])
    return arr


def _netcdf_layout(ds):
    """
    Work out how a SOCAT NetCDF file is laid out.

    Returns (obs_dim, traj_dim, row_size_var). SOCAT's full-resolution NetCDF
    from ERDDAP is a CF 'contiguous ragged array': cruise-level variables
    (Expocode, QC flag) sit on a 'trajectory' dimension and each cruise's
    measurements are stored back to back on an 'obs' dimension, with
    rowSize saying how many belong to each cruise. A flat table (one
    dimension) is also handled.
    """
    lat_name = _find_col(list(ds.variables), _SOCAT_CANDIDATES["lat"])
    lon_name = _find_col(list(ds.variables), _SOCAT_CANDIDATES["lon"])
    if lat_name is None or lon_name is None:
        raise ValueError(f"No latitude/longitude variables in NetCDF: {list(ds.variables)}")
    lat_dims = ds[lat_name].dims
    # Gridded product: latitude and longitude are separate axes of a grid.
    if len(lat_dims) == 1 and lat_dims[0] == lat_name and ds[lon_name].dims == (lon_name,):
        raise ValueError(
            "This looks like the GRIDDED SOCAT product (1-degree monthly means). "
            "It is too coarse for shelf-scale clustering; download the "
            "full-resolution (per-measurement) data instead.")
    if len(lat_dims) != 1:
        raise ValueError(f"Unexpected latitude dimensions {lat_dims}")
    obs_dim = lat_dims[0]
    row_size = next((v for v in ds.variables
                     if v.lower() in ("rowsize", "row_size")), None)
    traj_dim = ds[row_size].dims[0] if row_size is not None else None
    return obs_dim, traj_dim, row_size


def _netcdf_to_frame(ds, obs_index=None):
    """
    Flatten a SOCAT NetCDF dataset to one row per measurement.

    obs_index: optional integer array of measurements to keep. Only those are
    read from disk, so a subset of a multi-GB file stays small in memory.
    """
    obs_dim, traj_dim, row_size = _netcdf_layout(ds)
    n_obs = ds.sizes[obs_dim]
    idx = np.arange(n_obs) if obs_index is None else np.asarray(obs_index)

    traj_of_obs = None
    if row_size is not None:
        # cruise number for each measurement, from rowSize
        traj_of_obs = np.repeat(np.arange(ds.sizes[traj_dim]),
                                ds[row_size].values.astype(int))[idx]

    cols = {}
    for name, var in ds.variables.items():
        if var.ndim == 1 and var.dims[0] == obs_dim and name != obs_dim:
            cols[name] = _decode_text(var.isel({obs_dim: idx}).values)
        elif (traj_of_obs is not None and var.ndim == 1 and var.dims[0] == traj_dim
              and name != row_size):
            cols[name] = _decode_text(var.values)[traj_of_obs]
        elif var.ndim == 2 and var.dims[0] in (obs_dim, traj_dim) and var.dtype.kind == "S":
            # fixed-width char arrays, e.g. expocode(trajectory, string_length)
            joined = np.array([b"".join(r).strip(b"\x00 ") for r in var.values])
            vals = _decode_text(joined)
            cols[name] = vals[idx] if var.dims[0] == obs_dim else vals[traj_of_obs]
    return pd.DataFrame(cols)


def _read_socat_netcdf(path):
    import xarray as xr  # only needed for NetCDF input
    with xr.open_dataset(path, mask_and_scale=True) as ds:
        return _netcdf_to_frame(ds)


def _read_socat_table(path):
    """Read a SOCAT text file, skipping the metadata block above the header."""
    if _is_netcdf(path):
        return _read_socat_netcdf(path)
    header_idx, sep = _read_header_index(path)
    # Read everything as text first: Expocodes and QC flags must stay strings.
    df = pd.read_csv(path, sep=sep, skiprows=header_idx, dtype=str, low_memory=False)
    # Fix: files saved from Excel carry many fully blank trailing rows
    # (",,,,,"); the 'socat head.csv' sample had ~142k of them.
    df = df.dropna(how="all")
    # ERDDAP CSVs carry a units row directly under the header
    if len(df) and pd.to_numeric(df.iloc[0], errors="coerce").isna().all():
        df = df.iloc[1:]
    return df.reset_index(drop=True)


def _socat_dates(df):
    """
    Build timestamps from the separate year/month/day/hour/minute/second
    columns when present (case-insensitive), because the combined DATETIME
    column is often reformatted by Excel (e.g. '5/12/2009 18:17').
    """
    lower = {str(c).lower(): c for c in df.columns}
    parts = {"year": ["year", "yr"], "month": ["month", "mon"], "day": ["day"],
             "hour": ["hour", "hh"], "minute": ["minute", "mm"], "second": ["second", "ss"]}
    found = {}
    for key, names in parts.items():
        for n in names:
            if n in lower:
                found[key] = pd.to_numeric(df[lower[n]], errors="coerce")
                break
    if {"year", "month", "day"} <= found.keys():
        return pd.to_datetime(pd.DataFrame(found), errors="coerce")
    tcol = _find_col(df.columns, _SOCAT_CANDIDATES["time"])
    if tcol is None:
        raise ValueError("No date information found in SOCAT file")
    t = pd.to_datetime(df[tcol], errors="coerce", utc=True)
    return t.dt.tz_localize(None)


def fco2_to_pco2(fco2, sst_c, p_atm=1.0):
    """
    Fugacity -> partial pressure (Weiss 1974; SOCAT convention).
    pCO2 is about 0.3% higher than fCO2 at surface temperatures.
    """
    T = np.asarray(sst_c, dtype=float) + 273.15
    B = -1636.75 + 12.0408 * T - 3.27957e-2 * T**2 + 3.16528e-5 * T**3
    delta = 57.7 - 0.118 * T
    R = 82.0578  # cm3 atm / (mol K)
    return np.asarray(fco2, dtype=float) / np.exp(p_atm * (B + 2 * delta) / (R * T))


def load_socat(path, good_flags=(2,), good_qc=SOCAT_GOOD_QC):
    """
    Load a SOCAT extract: Data Set Viewer .tsv, synthesis file, ERDDAP .csv,
    the upper-case CSV export used for 'socat head.csv', or NetCDF (.nc).

      * -1E+34 fill values become NaN
      * SOCAT stores longitude as 0-360 E in some exports; converted to -180..180
      * keeps only WOCE flag 2 (good) fCO2 and cruise QC flags A-D
      * drops rows with no fCO2 or fCO2 <= 0
      * converts fCO2 to pCO2 with fco2_to_pco2()
    """
    df = _read_socat_table(path)
    cols = {k: _find_col(df.columns, v) for k, v in _SOCAT_CANDIDATES.items()}
    missing = [k for k in ("lon", "lat", "sst", "fco2") if cols[k] is None]
    if missing:
        raise ValueError(f"Could not find SOCAT columns {missing}. "
                         f"Columns present: {list(df.columns)}")

    for k in ("lon", "lat", "sst", "sal", "fco2", "flag"):
        if cols[k] is not None:
            x = pd.to_numeric(df[cols[k]], errors="coerce")
            df[cols[k]] = x.where(x > SOCAT_FILL_THRESHOLD)

    n0 = len(df)
    log = [f"{n0} rows read"]
    if cols["flag"] is not None and good_flags:
        df = df[df[cols["flag"]].isin(good_flags)]
        log.append(f"{len(df)} with WOCE flag in {good_flags}")
    if cols["qc"] is not None and good_qc:
        df = df[df[cols["qc"]].astype(str).str.strip().str.upper().isin(good_qc)]
        log.append(f"{len(df)} with QC flag in {''.join(good_qc)}")
    df = df[df[cols["fco2"]] > 0]
    log.append(f"{len(df)} with fCO2 > 0")
    print("[socat] " + " -> ".join(log))

    if cols["expo"] is not None:
        station = df[cols["expo"]].astype(str).str.strip()
        # Fix: Excel turns numeric-looking Expocodes (e.g. 316420090512) into
        # scientific notation like '3.16E+11', which merges separate cruises
        # under one ID. They can't be recovered from the file, so flag them.
        mangled = station.str.contains(r"E\+\d+$", regex=True)
        if mangled.any():
            print(f"[socat] WARNING: {int(mangled.sum())} rows have Expocodes in "
                  f"scientific notation ({sorted(station[mangled].unique())}); "
                  f"Excel mangled them. Re-download without opening in Excel "
                  f"to keep cruise IDs.")
    else:
        station = "SOCAT"

    lon = df[cols["lon"]].where(df[cols["lon"]] <= 180, df[cols["lon"]] - 360)

    out = pd.DataFrame({
        "Date": _socat_dates(df),
        "Latitude": df[cols["lat"]],
        "Longitude": lon,
        "Temp": df[cols["sst"]],
        "Sal": df[cols["sal"]] if cols["sal"] is not None else np.nan,
        "pCO2": fco2_to_pco2(df[cols["fco2"]], df[cols["sst"]]),
        "Source": "SOCAT",
        "Station": station,
    })
    return out.dropna(subset=["Latitude", "Longitude"]).reset_index(drop=True)
