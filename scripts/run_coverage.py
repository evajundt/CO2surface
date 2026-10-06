"""
Run the coverage diagnostic on everything currently available.

    python scripts/run_coverage.py                      # uses data/raw/socat_nwgom.tsv.gz
    python scripts/run_coverage.py --socat data/raw/<full_socat_export>.tsv

Writes a coverage map to outputs/coverage_map.png.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from loaders import load_reef_loggers, load_ship, load_socat
from nwgom_coverage_check import DOMAIN, coverage_report

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "outputs"

# File names as uploaded to the repo (spaces included)
SHIP_FILE = RAW / "shipto2023.xlsx"
LOGGER_FILE = RAW / "Temp E_W_1989-2024_1sheet.xlsx"
SOCAT_DEFAULT = RAW / "socat_nwgom.tsv.gz"   # SOCAT v2026 subset (was "socat head.csv")


def nearest_neighbour_km(df):
    """
    Distance from each unique location to its nearest other location (km).

    Uses a KD-tree on 3-D unit vectors rather than a full n x n distance
    matrix, so it stays fast and small for the full SOCAT file
    (a full matrix for 100k points would need ~80 GB of memory).
    """
    pts = df[["Latitude", "Longitude"]].round(4).drop_duplicates().to_numpy()
    if len(pts) < 2:
        return np.array([])
    lat, lon = np.radians(pts[:, 0]), np.radians(pts[:, 1])
    xyz = np.column_stack([np.cos(lat) * np.cos(lon),
                           np.cos(lat) * np.sin(lon),
                           np.sin(lat)])
    chord, _ = cKDTree(xyz).query(xyz, k=2)
    return 2 * 6371 * np.arcsin(chord[:, 1] / 2)   # chord length -> great-circle km


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--socat", default=str(SOCAT_DEFAULT),
                    help="path to a SOCAT export (.tsv/.csv)")
    args = ap.parse_args()

    datasets = {
        "Ship (surface)": load_ship(SHIP_FILE),
        "Reef loggers (~20 m)": load_reef_loggers(LOGGER_FILE),
    }
    if args.socat and Path(args.socat).exists():
        datasets["SOCAT"] = load_socat(args.socat)
    else:
        print(f"[socat] no file at {args.socat}; skipping SOCAT")

    for label, df in datasets.items():
        d = coverage_report(df, label=label)
        nn = nearest_neighbour_km(d)
        if len(nn):
            print(f"  Nearest-neighbour distance between unique locations (km): "
                  f"median {np.median(nn):.1f}, 90th pct {np.percentile(nn, 90):.1f}, "
                  f"max {nn.max():.1f}\n")
        if "Station" in d.columns and d["Station"].nunique() <= 30:
            # For SOCAT, Station is the cruise Expocode
            print(d.groupby("Station").agg(n=("Date", "size"),
                                            first=("Date", "min"),
                                            last=("Date", "max"),
                                            lon_min=("Longitude", "min"),
                                            lon_max=("Longitude", "max")).to_string(), "\n")
        if d["Sal"].notna().any():
            # Low-salinity points are river-plume / bay water; they can
            # dominate a pCO2 clustering, so it's worth knowing how many there are.
            print(f"  Points with salinity < 25: {int((d['Sal'] < 25).sum())} "
                  f"of {int(d['Sal'].notna().sum())}\n")

    OUT.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 8))
    styles = {"Ship (surface)": dict(marker="o", s=60),
              "Reef loggers (~20 m)": dict(marker="^", s=120),
              "SOCAT": dict(marker=".", s=4, alpha=0.4)}
    for label, df in datasets.items():
        # Fix: legend n counts only points inside DOMAIN (SOCAT files extend past it)
        df = df[df["Longitude"].between(DOMAIN["lon_min"], DOMAIN["lon_max"]) &
                df["Latitude"].between(DOMAIN["lat_min"], DOMAIN["lat_max"])]
        d = df.drop_duplicates(subset=["Latitude", "Longitude"])
        ax.scatter(d["Longitude"], d["Latitude"], label=f"{label} (n={len(df):,})",
                   **styles.get(label, {}))
    ax.plot(-94.80, 29.30, "k*", ms=12)
    ax.annotate("Galveston", (-94.80, 29.30), xytext=(5, 5), textcoords="offset points")
    ax.set_xlim(DOMAIN["lon_min"], DOMAIN["lon_max"])
    ax.set_ylim(DOMAIN["lat_min"], DOMAIN["lat_max"])
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("NW Gulf coverage - unique sampling locations")
    ax.legend(frameon=False)
    plt.tight_layout()
    fig.savefig(OUT / "coverage_map.png", dpi=150)
    print(f"Saved {OUT / 'coverage_map.png'}")


if __name__ == "__main__":
    main()
