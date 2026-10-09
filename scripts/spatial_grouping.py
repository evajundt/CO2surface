"""
Spatial grouping of surface pCO2 and temperature in the NW Gulf.

Planned two-stage analysis (run separately for each variable):

  1. Bin SOCAT measurements into grid cells. Each cruise's pass through a cell
     is one "visit" (mean of its points). This removes the along-track
     redundancy: consecutive points ~0.2 km apart are not independent
     samples and would make Moran's I significant almost by default.
  2. Remove the seasonal cycle and long-term trend with Eva's GAM
     (seasonal.compute_seasonal_trend_gam, fitted on visits), so that
     "location" is not just "the month/year a cruise happened to pass
     there". Analysis is on the anomalies. --seasonal harmonic gives the
     earlier simpler fit for comparison.
  3. STAGE 1 - all data: global Moran's I + LISA (local clusters, FDR-corrected).
  4. JUSTIFICATION - do the LISA clusters coincide with low-salinity water?
     (Kruskal-Wallis / Mann-Whitney on salinity by cluster type, chi-square of
     cluster type x fresh/marine.)
  5. STAGE 2a - REMOVAL: drop visits with salinity below a threshold fixed in
     advance (primary 25, sensitivity 30 and 33), rebuild neighbours, rerun.
  5. STAGE 2b - ADJUSTMENT: fit anomaly vs salinity with a LOWESS curve and
     rerun Moran's I / LISA on the residuals (keeps all data, no threshold).
     A straight-line fit is also run for comparison.
     SENSITIVITY - cruise-centred anomalies (each cruise's mean removed), to
     check clusters aren't just one unusually warm/high-CO2 cruise's track.
  6. COMPARE removal vs adjustment: Moran's I side by side, and agreement of
     the LISA cluster labels on cells both keep (% agreement, Cohen's kappa).

    python scripts/spatial_grouping.py                 # 0.1-degree cells
    python scripts/spatial_grouping.py --cell 0.25     # coarser sensitivity run

Outputs go to outputs/spatial_<cell>deg/.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess

import esda
from libpysal.weights import DistanceBand

from loaders import load_socat
from seasonal import compute_seasonal_trend_gam, iqr_outliers
from nwgom_coverage_check import DOMAIN

ROOT = Path(__file__).resolve().parents[1]
SOCAT_FILE = ROOT / "data" / "raw" / "socat_nwgom.tsv.gz"

VARIABLES = {"pCO2": "pCO2 (uatm)", "Temp": "Temperature (deg C)"}
PRIMARY_SAL = 25.0                   # fixed before looking at results
SENSITIVITY_SAL = (30.0, 33.0)
PERMUTATIONS = 999
FDR_ALPHA = 0.05
LOWESS_FRAC = 0.3                    # smoothing span for the salinity adjustment
SEED = 12345                         # LISA permutations are random; fixed for reproducibility

# LISA colours (validated: CVD and normal-vision separation pass; the two
# lighter outlier colours are < 3:1 contrast, so shape + legend + CSV carry it too)
LISA_STYLE = {
    "High-High": dict(color="#b8312f", marker="s"),
    "Low-Low":   dict(color="#1c5aa8", marker="s"),
    "High-Low":  dict(color="#ec7b7a", marker="D"),
    "Low-High":  dict(color="#6aa3e8", marker="D"),
    "ns":        dict(color="#c9c8c4", marker="."),
}


# ---------------------------------------------------------------------------
# 1. Visits and cells
# ---------------------------------------------------------------------------
def make_visits(df, cell):
    """One row per (cell, cruise): mean position, time, Temp, Sal, pCO2."""
    d = df.copy()
    d["cell_x"] = np.floor((d["Longitude"] - DOMAIN["lon_min"]) / cell).astype(int)
    d["cell_y"] = np.floor((d["Latitude"] - DOMAIN["lat_min"]) / cell).astype(int)
    d["t"] = d["Date"].astype("int64") / 1e9          # seconds, for averaging
    v = (d.groupby(["cell_x", "cell_y", "Station"])
           .agg(t=("t", "mean"), Temp=("Temp", "mean"), Sal=("Sal", "mean"),
                pCO2=("pCO2", "mean"), n_points=("pCO2", "size"))
           .reset_index())
    v["Date"] = pd.to_datetime(v["t"], unit="s")
    # cell centre
    v["lon"] = DOMAIN["lon_min"] + (v["cell_x"] + 0.5) * cell
    v["lat"] = DOMAIN["lat_min"] + (v["cell_y"] + 0.5) * cell
    return v


def deseasonalize(v, var):
    """
    Anomaly = value - (intercept + linear trend + 2 annual harmonics).

    Fitted on visits (not points), so a cruise with 10,000 points doesn't
    outweigh one with 200. This is a simple domain-wide climatology; the
    GAM approach in compute_seasonal_trend_gam can replace it later.
    """
    yrs = v["Date"].dt.year + (v["Date"].dt.dayofyear - 1) / 365.25
    w = 2 * np.pi * v["Date"].dt.dayofyear / 365.25
    X = np.column_stack([np.ones(len(v)), yrs - 2015,
                         np.cos(w), np.sin(w), np.cos(2 * w), np.sin(2 * w)])
    y = v[var].to_numpy()
    ok = np.isfinite(y)
    beta, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None)
    fit = X @ beta
    r2 = 1 - np.nansum((y - fit) ** 2) / np.nansum((y - np.nanmean(y)) ** 2)
    info = {"trend_per_yr": beta[1], "seasonal_amplitude": np.hypot(beta[2], beta[3]),
            "r2": r2}
    return y - fit, info


def to_cells(v, value_col):
    """
    Average visits into cells; keep the number of visits per cell and the
    share of visits that were fresh (S < PRIMARY_SAL).

    Fix: a cell's MEAN salinity hides occasional plume water (no 0.1-degree
    cell has mean S < 25), so fresh influence is tracked per visit instead.
    """
    v = v.dropna(subset=[value_col]).assign(fresh=lambda x: x["Sal"] < PRIMARY_SAL)
    return (v.groupby(["cell_x", "cell_y", "lon", "lat"])
             .agg(value=(value_col, "mean"), Sal=("Sal", "mean"),
                  min_Sal=("Sal", "min"), frac_fresh=("fresh", "mean"),
                  n_visits=(value_col, "size"))
             .reset_index())


# ---------------------------------------------------------------------------
# 2. Spatial statistics
# ---------------------------------------------------------------------------
def project_km(lon, lat):
    """Simple equirectangular projection to km around the domain centre."""
    lat0 = (DOMAIN["lat_min"] + DOMAIN["lat_max"]) / 2
    lon0 = (DOMAIN["lon_min"] + DOMAIN["lon_max"]) / 2
    x = (np.asarray(lon) - lon0) * 111.32 * np.cos(np.radians(lat0))
    y = (np.asarray(lat) - lat0) * 110.57
    return np.column_stack([x, y])


def build_weights(cells, band_km):
    """
    Distance-band neighbours (binary, then row-standardised). Cells with no
    neighbour inside the band ("islands") are dropped and reported, because
    LISA can't say anything about a cell with no neighbours.
    """
    coords = project_km(cells["lon"], cells["lat"])
    w = DistanceBand(coords, threshold=band_km, binary=True, silence_warnings=True)
    islands = list(w.islands)
    if islands:
        cells = cells.drop(index=cells.index[islands]).reset_index(drop=True)
        coords = project_km(cells["lon"], cells["lat"])
        w = DistanceBand(coords, threshold=band_km, binary=True, silence_warnings=True)
    w.transform = "r"
    return cells, w, len(islands)


def run_spatial(cells, band_km):
    """Global Moran's I and FDR-corrected LISA on cells['value']."""
    cells, w, n_islands = build_weights(cells, band_km)
    y = cells["value"].to_numpy()
    np.random.seed(SEED)              # esda.Moran has no seed argument
    mi = esda.Moran(y, w, permutations=PERMUTATIONS)
    lisa = esda.Moran_Local(y, w, permutations=PERMUTATIONS, seed=SEED)
    # Benjamini-Hochberg style threshold across all local tests
    p_cut = esda.fdr(lisa.p_sim, FDR_ALPHA)
    labels = np.array(["ns"] * len(y), dtype=object)
    names = {1: "High-High", 2: "Low-High", 3: "Low-Low", 4: "High-Low"}
    sig = lisa.p_sim <= p_cut
    for q, name in names.items():
        labels[sig & (lisa.q == q)] = name
    cells = cells.assign(lisa=labels, lisa_p=lisa.p_sim)
    summary = {"n_cells": len(cells), "islands_dropped": n_islands,
               "morans_I": mi.I, "expected_I": mi.EI, "z": mi.z_sim, "p": mi.p_sim,
               "fdr_p_cut": p_cut}
    for name in names.values():
        summary[f"n_{name}"] = int((labels == name).sum())
    return cells, summary


# ---------------------------------------------------------------------------
# 3. Justification and comparison
# ---------------------------------------------------------------------------
def salinity_vs_clusters(cells, sal_cut):
    """Does salinity differ between LISA cluster types? Are clusters 'fresh'?"""
    out = {}
    groups = {k: g["Sal"].dropna().to_numpy() for k, g in cells.groupby("lisa")}
    present = [g for g in groups.values() if len(g) >= 3]
    if len(present) >= 2:
        out["kruskal_H"], out["kruskal_p"] = stats.kruskal(*present)
    in_cluster = cells.loc[cells["lisa"] != "ns", "Sal"].dropna()
    not_cluster = cells.loc[cells["lisa"] == "ns", "Sal"].dropna()
    if len(in_cluster) >= 3 and len(not_cluster) >= 3:
        out["mannwhitney_U"], out["mannwhitney_p"] = stats.mannwhitneyu(
            in_cluster, not_cluster, alternative="two-sided")
    # cell is "fresh-influenced" if ANY visit had S < sal_cut
    table = pd.crosstab(cells["lisa"], (cells["min_Sal"] < sal_cut).rename(f"any visit S<{sal_cut:g}"))
    if table.shape[1] == 2 and table.shape[0] >= 2:
        out["chi2"], out["chi2_p"], _, _ = stats.chi2_contingency(table)
    out["median_sal_by_cluster"] = {k: float(np.median(g)) for k, g in groups.items() if len(g)}
    rho, p = stats.spearmanr(cells["value"], cells["Sal"], nan_policy="omit")
    out["spearman_value_vs_sal"], out["spearman_p"] = rho, p
    return out, table


def cohens_kappa(a, b):
    cats = sorted(set(a) | set(b))
    m = pd.crosstab(pd.Categorical(a, cats), pd.Categorical(b, cats), dropna=False).to_numpy()
    n = m.sum()
    po = np.trace(m) / n
    pe = (m.sum(0) * m.sum(1)).sum() / n**2
    return po, (po - pe) / (1 - pe) if pe < 1 else np.nan


def compare_labels(c1, c2):
    """Agreement of LISA labels on cells present in both runs."""
    m = c1.merge(c2, on=["cell_x", "cell_y"], suffixes=("_a", "_b"))
    if not len(m):
        return {"shared_cells": 0}
    agree, kappa = cohens_kappa(m["lisa_a"].to_numpy(), m["lisa_b"].to_numpy())
    return {"shared_cells": len(m), "pct_agree": 100 * agree, "kappa": kappa}


# ---------------------------------------------------------------------------
# 4. Plot
# ---------------------------------------------------------------------------
def plot_lisa_panels(results, var, label, path):
    fig, axes = plt.subplots(1, len(results), figsize=(5.2 * len(results), 5.4),
                             sharex=True, sharey=True)
    for ax, (title, (cells, summ)) in zip(np.atleast_1d(axes), results.items()):
        for name in ["ns", "Low-High", "High-Low", "Low-Low", "High-High"]:
            c = cells[cells["lisa"] == name]
            st = LISA_STYLE[name]
            ax.scatter(c["lon"], c["lat"], c=st["color"], marker=st["marker"],
                       s=10 if name == "ns" else 16, linewidths=0)
        ax.set_title(f"{title}\nMoran's I = {summ['morans_I']:.2f} "
                     f"(p = {summ['p']:.3f}), n = {summ['n_cells']}", fontsize=9)
        ax.plot(-94.80, 29.30, "k*", ms=8)
        ax.set_xlim(DOMAIN["lon_min"], DOMAIN["lon_max"])
        ax.set_ylim(DOMAIN["lat_min"], DOMAIN["lat_max"])
        ax.set_xlabel("Longitude")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    np.atleast_1d(axes)[0].set_ylabel("Latitude")
    handles = [Line2D([], [], linestyle="", marker=LISA_STYLE[n]["marker"],
                      color=LISA_STYLE[n]["color"], markersize=7,
                      label=n if n != "ns" else "not significant")
               for n in ["High-High", "Low-Low", "High-Low", "Low-High", "ns"]]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False)
    fig.suptitle(f"LISA clusters of {label} anomaly (FDR {FDR_ALPHA})", fontsize=11)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--socat", default=str(SOCAT_FILE))
    ap.add_argument("--cell", type=float, default=0.1, help="grid cell size (degrees)")
    ap.add_argument("--band-km", type=float, default=None,
                    help="neighbour distance; default 2.5 cells")
    ap.add_argument("--seasonal", choices=["gam", "harmonic"], default="gam",
                    help="gam = Eva's compute_seasonal_trend_gam (default); "
                         "harmonic = the earlier 2-harmonic + linear trend fit")
    ap.add_argument("--iqr-clean", action="store_true",
                    help="drop visits whose anomaly is a 3x IQR outlier (off by default: "
                         "it also removes plume water, and SOCAT is already QC'd)")
    args = ap.parse_args()
    band_km = args.band_km or 2.5 * args.cell * 111.0
    suffix = "" if args.seasonal == "gam" else f"_{args.seasonal}"
    suffix += "_iqrclean" if args.iqr_clean else ""
    out = ROOT / "outputs" / f"spatial_{args.cell:g}deg{suffix}"
    out.mkdir(parents=True, exist_ok=True)

    df = load_socat(args.socat)
    df = df[df["Longitude"].between(DOMAIN["lon_min"], DOMAIN["lon_max"]) &
            df["Latitude"].between(DOMAIN["lat_min"], DOMAIN["lat_max"])]
    visits = make_visits(df, args.cell)
    print(f"{len(df):,} points -> {len(visits):,} visits in "
          f"{visits.groupby(['cell_x', 'cell_y']).ngroups:,} cells "
          f"({args.cell} deg); neighbour band {band_km:.1f} km")

    rows, notes = [], []
    for var, label in VARIABLES.items():
        if args.seasonal == "gam":
            # Eva's GAM (cyclic DOY spline + time spline), fitted on visits
            # so a cruise with 10,000 points doesn't outweigh one with 200
            fitted, fit = compute_seasonal_trend_gam(visits, var)
            anom = fitted[f"Anomaly_{var}"].to_numpy()
            method = "GAM (compute_seasonal_trend_gam, decimal-year time axis)"
        else:
            anom, fit = deseasonalize(visits, var)
            method = "2-harmonic + linear trend"
        n_iqr = 0
        if args.iqr_clean:
            bad = iqr_outliers(anom)
            n_iqr = int(bad.sum())
            anom = np.where(bad, np.nan, anom)
        visits[f"{var}_anom"] = anom
        notes.append(f"## {var}\n\nSeasonal+trend removal: {method}, fitted on visits. "
                     f"Trend {fit['trend_per_yr']:+.2f}/yr, seasonal amplitude "
                     f"{fit['seasonal_amplitude']:.2f}, R^2 {fit['r2']:.2f}"
                     + (f"; {n_iqr} visits dropped as 3x IQR anomaly outliers" if args.iqr_clean else "")
                     + "\n")
        if "monthly_mean_anomaly" in fit:
            notes.append("Residual seasonality check (mean anomaly by month; ~0 means the cycle "
                         "was removed): " + ", ".join(f"{m}: {a:+.1f}" for m, a in
                                                      fit["monthly_mean_anomaly"].items()) + "\n")

        results = {}
        # Stage 1: all data
        cells_all, s_all = run_spatial(to_cells(visits, f"{var}_anom"), band_km)
        results["All data"] = (cells_all, s_all)

        # Sensitivity: cruise-centred anomalies. Each cruise covers a region in
        # a few days, so an unusually warm/high-CO2 cruise shifts every cell it
        # visits and can look like a spatial cluster. Subtracting each cruise's
        # mean anomaly keeps only within-cruise spatial contrasts.
        visits[f"{var}_cc"] = (visits[f"{var}_anom"] -
                               visits.groupby("Station")[f"{var}_anom"].transform("mean"))
        results["Cruise-centred"] = run_spatial(to_cells(visits, f"{var}_cc"), band_km)
        just, table = salinity_vs_clusters(cells_all, PRIMARY_SAL)

        # Stage 2a: removal at primary + sensitivity thresholds
        for thr in (PRIMARY_SAL,) + SENSITIVITY_SAL:
            v = visits[visits["Sal"] >= thr]
            c, s = run_spatial(to_cells(v, f"{var}_anom"), band_km)
            s["visits_removed"] = int((visits["Sal"] < thr).sum())
            results[f"Removed S < {thr:g}"] = (c, s)

        # Stage 2b: adjustment - residuals of anomaly ~ salinity (visit level).
        # Main version is LOWESS (curved). Fix: the pCO2-salinity relation is
        # flat above S~32 and drops steeply below it, so a straight line
        # over-corrects offshore water and flips the map; the linear version
        # is kept only as a comparison.
        ok = visits[[f"{var}_anom", "Sal"]].notna().all(axis=1)
        reg = stats.linregress(visits.loc[ok, "Sal"], visits.loc[ok, f"{var}_anom"])
        visits[f"{var}_resid_lin"] = visits[f"{var}_anom"] - (reg.intercept + reg.slope * visits["Sal"])
        results["Salinity-adjusted (linear)"] = run_spatial(
            to_cells(visits, f"{var}_resid_lin"), band_km)
        smooth = lowess(visits.loc[ok, f"{var}_anom"], visits.loc[ok, "Sal"],
                        frac=LOWESS_FRAC, return_sorted=False)
        visits.loc[ok, f"{var}_resid"] = visits.loc[ok, f"{var}_anom"] - smooth
        c_adj, s_adj = run_spatial(to_cells(visits, f"{var}_resid"), band_km)
        results["Salinity-adjusted"] = (c_adj, s_adj)
        sal_bins = pd.cut(visits.loc[ok, "Sal"], [0, 25, 28, 30, 32, 33, 34, 35, 36, 40])
        shape = visits.loc[ok].groupby(sal_bins, observed=True)[f"{var}_anom"].agg(["mean", "size"])

        # Comparisons
        prim = f"Removed S < {PRIMARY_SAL:g}"
        cmp_rows = {
            f"{prim} vs Salinity-adjusted": compare_labels(results[prim][0], c_adj),
            "All data vs Cruise-centred": compare_labels(cells_all, results["Cruise-centred"][0]),
            "Salinity-adjusted vs Salinity-adjusted (linear)": compare_labels(
                c_adj, results["Salinity-adjusted (linear)"][0]),
            f"All data vs {prim}": compare_labels(cells_all, results[prim][0]),
            "All data vs Salinity-adjusted": compare_labels(cells_all, c_adj),
        }
        for thr in SENSITIVITY_SAL:
            cmp_rows[f"{prim} vs Removed S < {thr:g}"] = compare_labels(
                results[prim][0], results[f"Removed S < {thr:g}"][0])

        for scen, (c, s) in results.items():
            rows.append({"variable": var, "scenario": scen, **s})
            c.to_csv(out / f"lisa_{var}_{scen.replace(' ', '_').replace('<', 'lt')}.csv",
                     index=False)
        plot_lisa_panels({k: results[k] for k in ["All data", prim, "Salinity-adjusted",
                                                  "Cruise-centred"]},
                         var, label, out / f"lisa_map_{var}.png")

        # Report text for this variable
        notes.append(f"Salinity adjustment: LOWESS (frac {LOWESS_FRAC}) of anomaly on salinity, "
                     f"visit level, n {ok.sum()}. For comparison, a straight line gives slope "
                     f"{reg.slope:+.2f} per salinity unit, r^2 {reg.rvalue**2:.2f}.\n")
        notes.append("Mean anomaly by salinity band (why a straight line doesn't fit):\n\n"
                     f"```\n{shape.round(1).T.to_string()}\n```\n")
        notes.append("### Justification: do Stage 1 clusters coincide with fresh water?\n")
        notes.append(f"- Spearman rho(anomaly, salinity) across cells: "
                     f"{just['spearman_value_vs_sal']:+.2f} (p {just['spearman_p']:.2g})")
        notes.append("- Median salinity by LISA class: " + ", ".join(
            f"{k} {v:.1f}" for k, v in just["median_sal_by_cluster"].items()))
        if "kruskal_p" in just:
            notes.append(f"- Kruskal-Wallis salinity across classes: H {just['kruskal_H']:.1f}, "
                         f"p {just['kruskal_p']:.2g}")
        if "mannwhitney_p" in just:
            notes.append(f"- Mann-Whitney salinity, clustered vs not: p {just['mannwhitney_p']:.2g}")
        if "chi2_p" in just:
            notes.append(f"- Chi-square cluster class x (S < {PRIMARY_SAL:g}): "
                         f"chi2 {just['chi2']:.1f}, p {just['chi2_p']:.2g}")
        notes.append("\nCross-tab (columns: cell had at least one visit with salinity < "
                     f"{PRIMARY_SAL:g}):\n\n```\n{table.to_string()}\n```\n")
        notes.append("### Agreement of LISA labels between scenarios\n")
        notes.append("| comparison | shared cells | % agree | Cohen's kappa |\n|---|---|---|---|")
        for k, r in cmp_rows.items():
            if r.get("shared_cells"):
                notes.append(f"| {k} | {r['shared_cells']} | {r['pct_agree']:.1f} | {r['kappa']:.2f} |")
        notes.append("")

    summary = pd.DataFrame(rows)
    summary.to_csv(out / "summary.csv", index=False)
    cols = ["variable", "scenario", "n_cells", "visits_removed", "morans_I", "z", "p",
            "n_High-High", "n_Low-Low", "n_High-Low", "n_Low-High"]
    table_md = summary[cols].to_markdown(index=False, floatfmt=".3f")
    report = (f"# Spatial grouping results ({args.cell:g} deg cells)\n\n"
              f"Input: {len(df):,} SOCAT points in DOMAIN -> {len(visits):,} cruise-cell visits. "
              f"Seasonal removal: {args.seasonal}"
              f"{' + anomaly IQR cleaning' if args.iqr_clean else ''}. "
              f"Neighbours: distance band {band_km:.1f} km, row-standardised. "
              f"{PERMUTATIONS} permutations; LISA FDR alpha {FDR_ALPHA}.\n\n"
              f"## Moran's I by scenario\n\n{table_md}\n\n" + "\n".join(notes))
    (out / "report.md").write_text(report)
    print(summary[cols].to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"\nWrote {out}/report.md, summary.csv, lisa_map_*.png")


if __name__ == "__main__":
    main()
