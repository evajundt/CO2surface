"""
Data-driven surface regions and their long-term rates (SOCAT, NW Gulf).

One variable per analysis (pCO2, then temperature), as in the original plan:

  1. Visits: SOCAT points binned to cells; one value per cruise per cell.
  2. Spatial field: Eva's GAM (seasonal cycle + trend, fitted on visits,
     seasonal.compute_seasonal_trend_gam) removed; each cell's MEAN anomaly is
     its typical offset from the domain-wide seasonal cycle. This is the
     same field the LISA analysis used.
  3. Stage 1 - fresh-water region, fixed in advance by salinity, not by the
     variable being clustered: a cell is "fresh-influenced" (region F) if any
     visit had S < FRESH_SAL (the same criterion as spatial_grouping.py).
     Justification reported: Mann-Whitney test of F vs other cells' offsets.
     Without this step Ward only peels off 1-7-cell plume patches for pCO2
     and never divides the shelf. --single-stage skips it.
     Stage 2 - marine regions: SKATER (spopt) on the remaining cells' means,
     with a neighbour graph so every region is one contiguous patch and a
     minimum size of MIN_REGION_CELLS cells (small extreme patches, e.g. a
     3-cell Galveston Bay outflow patch, are absorbed by a neighbour instead
     of blocking the split). Cross-checked against plain spatially
     constrained Ward clustering; agreement = adjusted Rand index.
  4. Number of regions k, chosen by a rule fixed in advance: the smallest k
     after which one more region adds less than K_GAIN of the explained
     variance (R^2 = between-region / total sum of squares), and every
     region must have >= MIN_REGION_CELLS cells and >= MIN_TREND_YEARS
     usable years for a trend. Override with --k.
  5. Rates per region: seasonal cycle removed with a region-specific cyclic
     day-of-year GAM (trend kept), annual means of visits (years with >=
     MIN_VISITS_PER_YEAR visits), Sen's slope per year and Mann-Kendall p
     (the larger of the plain and Hamed-Rao p, as in ch3_trends.py).
     Sensitivity: the same after removing each cell's mean offset, so a
     year that happened to sample only the high side of a region can't
     look like a trend.

Regions come from each cell's long-term MEAN; trends come from change over
TIME. These are different properties of the data, but the regions were
still drawn from the same observations, which should be stated.

    python scripts/regions.py                 # 0.25-degree cells
    python scripts/regions.py --k 4           # force 4 regions
    python scripts/regions.py --cell 0.1

Outputs: outputs/regions_<cell>deg/
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymannkendall as mk
from libpysal.weights import DistanceBand
from pygam import LinearGAM, s
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score

from loaders import load_socat
from nwgom_coverage_check import DOMAIN
from seasonal import compute_seasonal_trend_gam, decimal_year
from spatial_grouping import SOCAT_FILE, make_visits, project_km

ROOT = Path(__file__).resolve().parents[1]
VARIABLES = {"pCO2": "pCO2 (uatm)", "Temp": "Temperature (C)"}
UNITS = {"pCO2": "uatm", "Temp": "C"}

K_RANGE = range(2, 9)
K_GAIN = 0.05              # stop adding regions when the next adds < 5% of variance
MIN_REGION_CELLS = 8
MIN_TREND_YEARS = 10
MIN_VISITS_PER_YEAR = 3
MIN_CELL_VISITS = 2        # cells with fewer visits are too noisy to shape regions
FRESH_SAL = 25.0           # stage-1 fresh-influenced cells: any visit below this salinity
FRESH_COLOR = "#6b6a66"    # region F drawn in neutral grey: a water-mass class, not a patch

# Validated categorical palette (dataviz reference). Only 3 slots clear the
# all-pairs colour checks, so region numbers are also printed on the map.
REGION_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                 "#e87ba4", "#008300", "#4a3aa7", "#e34948"]


# ---------------------------------------------------------------------------
# Regions
# ---------------------------------------------------------------------------
def cell_field(visits, var):
    """Mean GAM anomaly per cell (the long-term spatial offset)."""
    fitted, _ = compute_seasonal_trend_gam(visits, var)
    v = visits.assign(anom=fitted[f"Anomaly_{var}"].to_numpy())
    return (v.dropna(subset=["anom"])
             .groupby(["cell_x", "cell_y", "lon", "lat"])
             .agg(value=("anom", "mean"), n_visits=("anom", "size"),
                  n_years=("Date", lambda d: d.dt.year.nunique()),
                  min_sal=("Sal", "min"))
             .reset_index())


def neighbour_graph(cells, cell_deg):
    """Queen-style contiguity on the grid: neighbours within ~1.6 cells."""
    band = 1.6 * cell_deg * 111.0
    coords = project_km(cells["lon"], cells["lat"])
    w = DistanceBand(coords, threshold=band, binary=True, silence_warnings=True)
    return w, band


def largest_component(cells, w):
    """Clustering needs one connected graph; keep the largest connected patch."""
    from scipy.sparse.csgraph import connected_components
    n, lab = connected_components(w.sparse, directed=False)
    if n == 1:
        return cells, np.zeros(len(cells), bool)
    main = np.bincount(lab).argmax()
    return cells, lab != main


def ward_regions(values, w, k):
    model = AgglomerativeClustering(n_clusters=k, linkage="ward", connectivity=w.sparse)
    return model.fit_predict(values.reshape(-1, 1))


def skater_regions(cells, w, k, floor=1):
    import geopandas as gpd
    from shapely.geometry import Point
    from spopt.region import Skater
    gdf = gpd.GeoDataFrame(cells[["value"]].copy(),
                           geometry=[Point(x, y) for x, y in zip(cells["lon"], cells["lat"])])
    model = Skater(gdf, w, ["value"], n_clusters=k, floor=floor, trace=False,
                   islands="increase")
    model.solve()
    return np.asarray(model.labels_)


def r2_between(values, labels):
    tss = ((values - values.mean()) ** 2).sum()
    wss = sum(((values[labels == g] - values[labels == g].mean()) ** 2).sum()
              for g in np.unique(labels))
    return 1 - wss / tss


def renumber_by_value(values, labels):
    """Region "1" = lowest mean value, so numbering means something."""
    order = pd.Series(values).groupby(labels).mean().sort_values().index
    remap = {old: str(new + 1) for new, old in enumerate(order)}
    return np.array([remap[l] for l in labels], dtype=object)


def region_color(reg):
    return FRESH_COLOR if reg == "F" else REGION_COLORS[(int(reg) - 1) % len(REGION_COLORS)]


def usable_years(visits_region):
    y = visits_region["Date"].dt.year.value_counts()
    return int((y >= MIN_VISITS_PER_YEAR).sum())


def assign_unclustered(cells_all, cells_core, labels_core):
    """Cells left out (few visits / disconnected) join the nearest region cell."""
    from scipy.spatial import cKDTree
    tree = cKDTree(project_km(cells_core["lon"], cells_core["lat"]))
    _, idx = tree.query(project_km(cells_all["lon"], cells_all["lat"]))
    return labels_core[idx]


# ---------------------------------------------------------------------------
# Trends
# ---------------------------------------------------------------------------
def seasonal_only(df, col):
    """Value minus a cyclic day-of-year GAM (trend left in), level preserved."""
    doy = df["Date"].dt.dayofyear.to_numpy(float)[:, None]
    y = df[col].to_numpy(float)
    gam = LinearGAM(s(0, basis="cp", n_splines=12, edge_knots=[1, 366])).fit(doy, y)
    seas = gam.predict(doy)
    return y - seas + seas.mean()


def annual_trend(years, means):
    """Sen's slope per year (on actual years) + cautious Mann-Kendall p."""
    years = np.asarray(years, float)
    means = np.asarray(means, float)
    i, j = np.triu_indices(len(years), k=1)
    sen = float(np.median((means[j] - means[i]) / (years[j] - years[i])))
    o = mk.original_test(means)
    h = mk.hamed_rao_modification_test(means)
    p = float(np.nanmax([o.p, h.p]))   # Hamed-Rao can return NaN when its variance correction fails
    return {"sen_slope_per_yr": sen, "mk_p": p, "mk_p_original": o.p,
            "mk_p_hamed_rao": h.p, "tau": o.Tau,
            "trend": "no trend" if p >= 0.05 else ("increasing" if o.s > 0 else "decreasing")}


def rates_differ(ann, include_fresh):
    """
    Do the regions' rates differ? ANCOVA on the annual means:
    annual_mean ~ year * region; the year:region interaction F-test asks
    whether the slopes differ, beyond the regions just having different
    levels. Annual means (not visits) are used so samples aren't
    pseudo-replicated.
    """
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm
    a = ann if include_fresh else ann[ann["region"] != "F"]
    if a["region"].nunique() < 2:
        return None
    full = smf.ols("annual_mean ~ year * C(region)", data=a).fit()
    reduced = smf.ols("annual_mean ~ year + C(region)", data=a).fit()
    tab = anova_lm(reduced, full)
    return {"regions": ", ".join(sorted(a["region"].unique())), "n_annual_means": len(a),
            "F": float(tab["F"].iloc[1]), "p_slopes_differ": float(tab["Pr(>F)"].iloc[1]),
            "common_slope_per_yr": float(reduced.params["year"]),
            "common_slope_se": float(reduced.bse["year"]),
            "common_slope_p": float(reduced.pvalues["year"])}


def region_trends(visits, var, region_of_cell):
    rows, ann_rows = [], []
    v = visits.merge(region_of_cell, on=["cell_x", "cell_y"]).dropna(subset=[var])
    for reg, g in v.groupby("region"):
        g = g.sort_values("Date").copy()
        g["des"] = seasonal_only(g, var)
        # Sensitivity: remove each cell's mean offset within the region
        g["des_cell"] = g["des"] - g.groupby(["cell_x", "cell_y"])["des"].transform("mean") \
            + g["des"].mean()
        yr = g["Date"].dt.year
        counts = yr.value_counts()
        ok_years = counts[counts >= MIN_VISITS_PER_YEAR].index
        row = {"variable": var, "region": reg, "cells": g[["cell_x", "cell_y"]].drop_duplicates().shape[0],
               "visits": len(g), "cruises": g["Station"].nunique(),
               "first_year": int(yr.min()), "last_year": int(yr.max()),
               "years_used": len(ok_years),
               "mean_value": g[var].mean(), "median_sal": g["Sal"].median()}
        for key, col in [("", "des"), ("cell_adj_", "des_cell")]:
            ann = g[yr.isin(ok_years)].groupby(yr)[col].mean().sort_index()
            if len(ann) >= 4:
                t = annual_trend(ann.index, ann.values)
                row.update({key + k: val for k, val in t.items()})
            if key == "":
                ann_rows.append(pd.DataFrame({"variable": var, "region": reg,
                                              "year": ann.index, "annual_mean": ann.values,
                                              "visits": counts.reindex(ann.index).values}))
        if row["years_used"] < MIN_TREND_YEARS:
            row["note"] = f"fewer than {MIN_TREND_YEARS} usable years"
        rows.append(row)
    return pd.DataFrame(rows), pd.concat(ann_rows, ignore_index=True)


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
def plot_regions(cells, var, k, trends, cell_deg, path):
    fig, ax = plt.subplots(figsize=(8.5, 7))
    for reg in sorted(cells["region"].unique()):
        c = cells[cells["region"] == reg]
        ax.scatter(c["lon"], c["lat"], marker="s", s=(cell_deg / 0.25) ** 2 * 150,
                   color=region_color(reg), linewidths=0)
        # direct label at the region's centre (colour is not the only cue)
        t = trends[trends["region"] == reg].iloc[0]
        lbl = str(reg)
        cx, cy = c["lon"].median(), c["lat"].median()
        ax.text(cx, cy, lbl, ha="center", va="center", fontsize=12, fontweight="bold",
                color="white", bbox=dict(boxstyle="circle,pad=0.25", fc="#1f1f1d", ec="none"))
    ax.plot(-94.80, 29.30, "k*", ms=10)
    ax.annotate("Galveston", (-94.80, 29.30), xytext=(5, 5), textcoords="offset points",
                fontsize=9)
    ax.set_xlim(DOMAIN["lon_min"], DOMAIN["lon_max"])
    ax.set_ylim(DOMAIN["lat_min"], DOMAIN["lat_max"])
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(f"{VARIABLES[var]}: {k} marine regions (numbered low to high offset)"
                 + (" + F = fresh-influenced" if (cells["region"] == "F").any() else ""),
                 fontsize=11)
    for sname in ("top", "right"):
        ax.spines[sname].set_visible(False)
    # legend-as-table: region, mean anomaly, rate
    lines = []
    for _, t in trends.iterrows():
        rate = (f"{t['sen_slope_per_yr']:+.2f} {UNITS[var]}/yr, p={t['mk_p']:.3f}"
                if pd.notna(t.get("sen_slope_per_yr")) else "too few years")
        lines.append(f"{t['region']}: offset {t['mean_offset']:+.1f}; {rate}")
    ax.text(0.01, 0.01, "\n".join(lines), transform=ax.transAxes, fontsize=8.5,
            va="bottom", ha="left", family="monospace",
            bbox=dict(fc="white", ec="#c9c8c4", alpha=0.9))
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_region_trends(ann, trends, var, path):
    regs = sorted(ann["region"].unique(), key=lambda r: (r != "F", r))
    fig, axes = plt.subplots(1, len(regs), figsize=(3.6 * len(regs), 3.6), sharey=True)
    for ax, reg in zip(np.atleast_1d(axes), regs):
        a = ann[ann["region"] == reg]
        t = trends[trends["region"] == reg].iloc[0]
        col = region_color(reg)
        ax.scatter(a["year"], a["annual_mean"], color=col, s=28, zorder=3)
        if pd.notna(t.get("sen_slope_per_yr")):
            b = np.median(a["annual_mean"] - t["sen_slope_per_yr"] * a["year"])
            xs = np.array([a["year"].min(), a["year"].max()])
            ax.plot(xs, t["sen_slope_per_yr"] * xs + b, color="#1f1f1d", lw=1.5)
            p = "p < 0.001" if t["mk_p"] < 0.001 else f"p = {t['mk_p']:.3f}"
            ax.set_title(f"Region {reg}: {t['sen_slope_per_yr']:+.2f} {UNITS[var]}/yr ({p})",
                         fontsize=9)
        else:
            ax.set_title(f"Region {reg}: too few years", fontsize=9)
        ax.set_xlabel("Year")
        for sname in ("top", "right"):
            ax.spines[sname].set_visible(False)
    np.atleast_1d(axes)[0].set_ylabel(f"{VARIABLES[var]}\n(seasonal cycle removed)")
    fig.suptitle(f"{VARIABLES[var]}: annual means and Sen's slope by region", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--socat", default=str(SOCAT_FILE))
    ap.add_argument("--cell", type=float, default=0.25)
    ap.add_argument("--k", type=int, default=None, help="force the number of marine regions")
    ap.add_argument("--single-stage", action="store_true",
                    help="skip the stage-1 fresh-water region and cluster all cells")
    args = ap.parse_args()
    out = ROOT / "outputs" / f"regions_{args.cell:g}deg"
    out.mkdir(parents=True, exist_ok=True)

    df = load_socat(args.socat)
    df = df[df["Longitude"].between(DOMAIN["lon_min"], DOMAIN["lon_max"]) &
            df["Latitude"].between(DOMAIN["lat_min"], DOMAIN["lat_max"])]
    visits = make_visits(df, args.cell)
    report = [f"# Data-driven regions ({args.cell:g} deg cells)\n",
              f"{len(df):,} SOCAT points -> {len(visits):,} cruise-cell visits.\n"]
    all_trends, all_ann, all_scan, all_diff = [], [], [], []

    for var in VARIABLES:
        cells_all = cell_field(visits, var)
        # Stage 1: fresh-influenced cells (salinity rule fixed in advance)
        fresh = (cells_all["min_sal"] < FRESH_SAL) & (not args.single_stage)
        if fresh.any():
            from scipy.stats import mannwhitneyu
            mw = mannwhitneyu(cells_all.loc[fresh, "value"], cells_all.loc[~fresh, "value"])
            stage1 = (f"Stage 1: {int(fresh.sum())} cells with any visit S < {FRESH_SAL:g} form "
                      f"region F. Their mean offset is {cells_all.loc[fresh, 'value'].mean():+.2f} "
                      f"vs {cells_all.loc[~fresh, 'value'].mean():+.2f} for the other cells "
                      f"(Mann-Whitney p = {mw.pvalue:.2g}), justifying a separate region before "
                      f"looking for marine groupings.\n")
        else:
            stage1 = "Stage 1 skipped (--single-stage).\n"
        marine = cells_all[~fresh]
        core = marine[marine["n_visits"] >= MIN_CELL_VISITS].reset_index(drop=True)
        w, band = neighbour_graph(core, args.cell)
        _, off = largest_component(core, w)
        core = core[~off].reset_index(drop=True)
        w, band = neighbour_graph(core, args.cell)
        vals = core["value"].to_numpy()

        # --- k scan --------------------------------------------------------
        scan = []
        for k in K_RANGE:
            try:
                lab = skater_regions(core, w, k, floor=MIN_REGION_CELLS)
            except Exception:
                break
            lab_all = assign_unclustered(marine, core, lab)
            roc = marine[["cell_x", "cell_y"]].assign(region=lab_all)
            vr = visits.merge(roc, on=["cell_x", "cell_y"])
            scan.append({"variable": var, "k": k, "r2": r2_between(vals, lab),
                         "min_cells": int(np.bincount(lab).min()),
                         "min_usable_years": min(usable_years(g) for _, g in vr.groupby("region"))})
        scan = pd.DataFrame(scan)
        scan["gain_next"] = scan["r2"].shift(-1) - scan["r2"]
        feasible = scan[(scan["min_cells"] >= MIN_REGION_CELLS) &
                        (scan["min_usable_years"] >= MIN_TREND_YEARS)]
        if args.k:
            k = args.k
            rule = "set with --k"
        else:
            elbow = feasible[feasible["gain_next"].fillna(0) < K_GAIN]
            k = int(elbow["k"].min()) if len(elbow) else int(feasible["k"].max())
            rule = (f"smallest feasible k whose next region adds < {K_GAIN:.0%} of variance "
                    f"(regions >= {MIN_REGION_CELLS} cells and >= {MIN_TREND_YEARS} usable years)")
        all_scan.append(scan)

        # --- final regions + cross-check -----------------------------------
        lab = renumber_by_value(vals, skater_regions(core, w, k, floor=MIN_REGION_CELLS))
        # cross-check: plain spatially constrained Ward on the same graph
        ari = adjusted_rand_score(lab, ward_regions(vals, w, k))
        cells_all["region"] = "F"
        cells_all.loc[~fresh, "region"] = assign_unclustered(marine, core, lab)
        offsets = cells_all.groupby("region").apply(
            lambda c: np.average(c["value"], weights=c["n_visits"]), include_groups=False)

        trends, ann = region_trends(visits, var, cells_all[["cell_x", "cell_y", "region"]])
        trends["mean_offset"] = trends["region"].map(offsets)
        diff = [d for d in (rates_differ(ann[ann["region"].isin(trends.loc[trends["years_used"] >= MIN_TREND_YEARS, "region"])], False),
                            rates_differ(ann, True)) if d]
        diff = pd.DataFrame(diff).assign(variable=var)
        all_diff.append(diff)
        all_trends.append(trends)
        all_ann.append(ann)
        cells_all.assign(variable=var).to_csv(out / f"cells_{var}.csv", index=False)
        plot_regions(cells_all, var, k, trends, args.cell, out / f"regions_map_{var}.png")
        plot_region_trends(ann, trends, var, out / f"regions_trends_{var}.png")

        cols = ["region", "cells", "visits", "cruises", "years_used", "mean_offset",
                "median_sal", "sen_slope_per_yr", "mk_p", "trend",
                "cell_adj_sen_slope_per_yr", "cell_adj_mk_p"]
        report += [f"## {var}\n", stage1,
                   f"Stage 2 clustering on {len(core)} marine cells with >= {MIN_CELL_VISITS} visits "
                   f"({len(marine) - len(core)} sparse or disconnected cells assigned to the "
                   f"nearest region afterwards). Neighbours within {band:.0f} km.\n",
                   f"k = {k}: {rule}. Agreement of these SKATER regions with plain Ward "
                   f"clustering: adjusted Rand "
                   f"index {ari:.2f} (1 = identical, 0 = chance).\n",
                   "k scan:\n", scan.to_markdown(index=False, floatfmt=".3f") + "\n",
                   "Rates per region (Sen's slope per year; mk_p = larger of plain and "
                   "Hamed-Rao Mann-Kendall p; cell_adj_ = after removing each cell's mean "
                   "offset):\n",
                   trends[[c for c in cols if c in trends]].to_markdown(index=False, floatfmt=".4g")
                   + "\n",
                   "Do the rates differ between regions? ANCOVA on annual means "
                   "(year x region interaction; first row marine regions with enough years, "
                   "second row including F). common_slope = one shared rate if they don't:\n",
                   diff.to_markdown(index=False, floatfmt=".4g") + "\n"]
        print(f"\n=== {var}: {stage1.strip()}\n    k = {k} marine regions ({rule}); Ward cross-check ARI {ari:.2f}")
        print(scan.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
        print(trends[[c for c in cols if c in trends]].to_string(
            index=False, float_format=lambda x: f"{x:.4g}"))
        print(diff.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    pd.concat(all_trends).to_csv(out / "region_trends.csv", index=False)
    pd.concat(all_ann).to_csv(out / "region_annual_means.csv", index=False)
    pd.concat(all_scan).to_csv(out / "k_scan.csv", index=False)
    pd.concat(all_diff).to_csv(out / "rate_differences.csv", index=False)
    (out / "report.md").write_text("\n".join(report))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
