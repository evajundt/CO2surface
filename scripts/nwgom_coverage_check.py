"""
Coverage diagnostic: run this BEFORE committing to a grid resolution for
Moran's I / LISA. Answers: do I actually have enough spatial spread to
support spatial clustering, or am I really just looking at a few repeated
cruise tracks plus some point sensors?
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Revised NW Gulf domain, extended east to include Flower Garden Banks
DOMAIN = {
    "lon_min": -97.3,
    "lon_max": -93.4,
    "lat_min": 25.8,
    "lat_max": 29.5,
}


def coverage_report(df: pd.DataFrame, lon_col="Longitude", lat_col="Latitude",
                     date_col="Date", label="dataset", cell_degs=(0.5, 0.25, 0.1, 0.05)):
    """
    Prints and plots coverage stats for one dataset within DOMAIN.
    cell_degs: resolutions to test for how many cells get populated —
               helps pick a grid size that doesn't leave most cells empty.
    """
    df = df.copy()
    mask = (
        df[lon_col].between(DOMAIN["lon_min"], DOMAIN["lon_max"]) &
        df[lat_col].between(DOMAIN["lat_min"], DOMAIN["lat_max"])
    )
    d = df.loc[mask]

    print("=" * 55)
    print(f"COVERAGE REPORT — {label}")
    print("=" * 55)
    print(f"  Points in domain:        {len(d)}")
    if date_col in d.columns and len(d):
        print(f"  Date range:               {d[date_col].min()} to {d[date_col].max()}")
        print(f"  Distinct years:           {d[date_col].dt.year.nunique()}")
    print(f"  Unique (lon,lat) pairs:   {d[[lon_col, lat_col]].drop_duplicates().shape[0]}")

    for cell in cell_degs:
        lon_bin = (d[lon_col] / cell).round() * cell
        lat_bin = (d[lat_col] / cell).round() * cell
        n_cells_total = int(((DOMAIN["lon_max"] - DOMAIN["lon_min"]) / cell) *
                             ((DOMAIN["lat_max"] - DOMAIN["lat_min"]) / cell))
        n_cells_populated = pd.DataFrame({"lon_bin": lon_bin, "lat_bin": lat_bin}).drop_duplicates().shape[0]
        pct = 100 * n_cells_populated / max(n_cells_total, 1)
        print(f"  Grid {cell:.2f}°: {n_cells_populated} / {n_cells_total} cells populated ({pct:.1f}%)")

    print("=" * 55 + "\n")
    return d


def plot_coverage_map(datasets: dict, lon_col="Longitude", lat_col="Latitude"):
    """
    datasets: {label: dataframe} — overlays all of them so you can SEE
    whether ship stations and SOCAT points actually overlap in space,
    or sit in disjoint parts of the domain.
    """
    fig, ax = plt.subplots(figsize=(9, 8))
    colors = plt.cm.tab10.colors
    for i, (label, df) in enumerate(datasets.items()):
        mask = (
            df[lon_col].between(DOMAIN["lon_min"], DOMAIN["lon_max"]) &
            df[lat_col].between(DOMAIN["lat_min"], DOMAIN["lat_max"])
        )
        d = df.loc[mask]
        ax.scatter(d[lon_col], d[lat_col], s=8, alpha=0.4,
                   color=colors[i % len(colors)], label=f"{label} (n={len(d)})")
    ax.set_xlim(DOMAIN["lon_min"], DOMAIN["lon_max"])
    ax.set_ylim(DOMAIN["lat_min"], DOMAIN["lat_max"])
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("NW Gulf coverage — all sources overlaid")
    ax.legend(frameon=False, markerscale=2)
    plt.tight_layout()
    plt.show()


# Example usage once you've pulled the real SOCAT Data Set Viewer extract:
#
# socat_nwgom = pd.read_csv('socat_nwgom_extract.csv')  # from Data Set Viewer
# socat_nwgom['Date'] = pd.to_datetime(socat_nwgom[['yr','mon','day']])
#
# coverage_report(socat_nwgom, label='SOCAT NW Gulf extract')
# coverage_report(ship,        label='Ship discrete')
#
# plot_coverage_map({
#     'SOCAT': socat_nwgom,
#     'Ship':  ship,
# })
