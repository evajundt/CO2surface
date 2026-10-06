"""
Run the coverage diagnostic on everything currently available.

    python scripts/run_coverage.py
    python scripts/run_coverage.py --socat data/raw/<your_socat_export>.tsv

Writes a coverage map to outputs/coverage_map.png.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from loaders import load_reef_loggers, load_ship, load_socat
from nwgom_coverage_check import DOMAIN, coverage_report

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "outputs"


def nearest_neighbour_km(df):
    """Distance from each unique location to its nearest other location (km)."""
    pts = df[["Latitude", "Longitude"]].round(4).drop_duplicates().to_numpy()
    if len(pts) < 2:
        return np.array([])
    lat = np.radians(pts[:, 0])[:, None]
    lon = np.radians(pts[:, 1])[:, None]
    a = (np.sin((lat - lat.T) / 2) ** 2 +
         np.cos(lat) * np.cos(lat.T) * np.sin((lon - lon.T) / 2) ** 2)
    dist = 2 * 6371 * np.arcsin(np.sqrt(a))
    np.fill_diagonal(dist, np.inf)
    return dist.min(axis=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--socat", help="path to a SOCAT export (.tsv/.csv)")
    args = ap.parse_args()

    datasets = {
        "Ship (surface)": load_ship(RAW / "shipto2023.xlsx"),
        "Reef loggers (~20 m)": load_reef_loggers(RAW / "Temp_E_W_1989-2024_1sheet.xlsx"),
    }
    if args.socat:
        datasets["SOCAT"] = load_socat(args.socat)

    for label, df in datasets.items():
        d = coverage_report(df, label=label)
        nn = nearest_neighbour_km(d)
        if len(nn):
            print(f"  Nearest-neighbour distance between unique locations (km): "
                  f"median {np.median(nn):.1f}, 90th pct {np.percentile(nn, 90):.1f}, "
                  f"max {nn.max():.1f}\n")
        if "Station" in d.columns and d["Station"].nunique() <= 20:
            print(d.groupby("Station").agg(n=("Date", "size"),
                                            first=("Date", "min"),
                                            last=("Date", "max")).to_string(), "\n")

    OUT.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 8))
    styles = {"Ship (surface)": dict(marker="o", s=60),
              "Reef loggers (~20 m)": dict(marker="^", s=120),
              "SOCAT": dict(marker=".", s=4, alpha=0.4)}
    for label, df in datasets.items():
        d = df.drop_duplicates(subset=["Latitude", "Longitude"])
        ax.scatter(d["Longitude"], d["Latitude"], label=f"{label} (n={len(df)})",
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
