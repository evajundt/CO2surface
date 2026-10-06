"""
Loaders that turn each raw source into one common surface-point table:

    Date | Latitude | Longitude | Temp | pCO2 | Source | Station

All sources end up with the same column names, so the coverage check and
the later Moran's I / LISA steps can be run on any one of them or on all
of them stacked together.
"""

import io
import re

import numpy as np
import pandas as pd

COMMON_COLS = ["Date", "Latitude", "Longitude", "Temp", "pCO2", "Source", "Station"]
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
        parts.append(pd.DataFrame({
            "Date": pd.to_datetime(d[meta["date_col"]]),
            "Latitude": meta["Latitude"],
            "Longitude": meta["Longitude"],
            "Temp": d[meta["temp_col"]].astype(float),
            "pCO2": np.nan,
            "Source": "ReefLogger",
            "Station": site,
        }))
    return pd.concat(parts, ignore_index=True)


# ---------------------------------------------------------------------------
# SOCAT
# ---------------------------------------------------------------------------
# Column names differ between the SOCAT synthesis files, the Data Set Viewer
# export and the ERDDAP export. The first match in each list is used.
_SOCAT_CANDIDATES = {
    "lon":  [r"^longitude", r"^lon"],
    "lat":  [r"^latitude", r"^lat"],
    "sst":  [r"^SST", r"^temp", r"^sea_surface_temp"],
    "sal":  [r"^sal", r"^SSS"],
    "fco2": [r"^fCO2rec", r"^fCO2_rec", r"^fco2_recommended", r"^fCO2"],
    "flag": [r"^fCO2rec_flag", r"^WOCE_CO2_water", r"^fCO2_flag"],
    "time": [r"^time$", r"^date"],
}


def _find_col(columns, patterns):
    for pat in patterns:
        for c in columns:
            if re.search(pat, str(c), flags=re.IGNORECASE):
                return c
    return None


def _read_socat_table(path):
    """Read a SOCAT text file, skipping the metadata block above the header."""
    with open(path, "r", errors="replace") as fh:
        lines = fh.readlines()
    header_idx = 0
    for i, line in enumerate(lines):
        low = line.lower()
        if low.startswith("expocode") or ("latitude" in low and "longitude" in low):
            header_idx = i
            break
    header = lines[header_idx]
    sep = "\t" if header.count("\t") >= header.count(",") else ","
    df = pd.read_csv(io.StringIO("".join(lines[header_idx:])), sep=sep,
                     low_memory=False)
    # ERDDAP CSVs carry a units row directly under the header
    if len(df) and pd.to_numeric(df.iloc[0], errors="coerce").isna().all():
        df = df.iloc[1:].reset_index(drop=True)
    return df


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


def load_socat(path, good_flags=(2,)):
    """
    Load a SOCAT extract (Data Set Viewer .tsv, synthesis file or ERDDAP .csv).

      * SOCAT stores longitude as 0-360 E; converted to -180..180
      * keeps only WOCE flag 2 (good) fCO2 when a flag column is present
      * converts fCO2rec to pCO2 with fco2_to_pco2()
    """
    df = _read_socat_table(path)
    cols = {k: _find_col(df.columns, v) for k, v in _SOCAT_CANDIDATES.items()}
    missing = [k for k in ("lon", "lat", "sst", "fco2") if cols[k] is None]
    if missing:
        raise ValueError(f"Could not find SOCAT columns {missing}. "
                         f"Columns present: {list(df.columns)}")

    for k in ("lon", "lat", "sst", "fco2", "flag"):
        if cols[k] is not None:
            df[cols[k]] = pd.to_numeric(df[cols[k]], errors="coerce")

    if cols["flag"] is not None and good_flags:
        before = len(df)
        df = df[df[cols["flag"]].isin(good_flags)]
        print(f"[socat] kept {len(df)}/{before} rows with flag in {good_flags}")

    if cols["time"] is not None:
        date = pd.to_datetime(df[cols["time"]], errors="coerce", utc=True).dt.tz_localize(None)
    else:
        parts = {k: pd.to_numeric(df[k], errors="coerce")
                 for k in ("yr", "mon", "day") if k in df.columns}
        date = pd.to_datetime(pd.DataFrame({"year": parts["yr"], "month": parts["mon"],
                                            "day": parts["day"]}), errors="coerce")

    lon = df[cols["lon"]].where(df[cols["lon"]] <= 180, df[cols["lon"]] - 360)
    station = df["Expocode"].astype(str) if "Expocode" in df.columns else "SOCAT"

    out = pd.DataFrame({
        "Date": date,
        "Latitude": df[cols["lat"]],
        "Longitude": lon,
        "Temp": df[cols["sst"]],
        "pCO2": fco2_to_pco2(df[cols["fco2"]], df[cols["sst"]]),
        "Source": "SOCAT",
        "Station": station,
    })
    return out.dropna(subset=["Latitude", "Longitude"]).reset_index(drop=True)
