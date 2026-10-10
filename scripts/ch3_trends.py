"""
Per-year long-term trends for the chapter 3 datasets in the repo
(ship surface bottles, HOBO/SBE reef-cap loggers East and West).

Replaces the trend numbers from plot_anomaly_with_trend in
reference/fianlch3analysis.py, which were (a) per SAMPLE not per year, and
(b) ~0 by construction when run on GAM anomalies. See seasonal.py.

Steps for each dataset, variable and period:
  1. Remove the seasonal cycle only (seasonal.long_term_trend: cyclic
     day-of-year GAM, refitted within each period), keeping the trend.
  2. Average to ANNUAL MEANS. Consecutive daily logger values are strongly
     autocorrelated, so treating ~9,000 days as independent made p-values
     absurdly small (~1e-250). Only reasonably complete years are used
     (HOBO: data in >= 9 of 12 months; ship: >= 3 samples).
  3. Mann-Kendall test on the annual means, with Sen's slope as the trend
     estimate. Both the plain test and the Hamed-Rao modified test (which
     corrects for autocorrelation between years) are run; the reported
     annual_mk_p is the LARGER of the two, because with ~10 years Hamed-Rao
     can estimate negative autocorrelation and make p smaller than it should.
  The ordinary regression on all samples is also reported, for reference
  only; its p-value is NOT autocorrelation-corrected.

Periods: full record, 2007-present (~last 15+ years), and the ship-sampling
window (Nov 2013 - Aug 2023) so the loggers can be compared like-for-like
with the ship bottles.

Ship depth sets (--ship-depths): surface (S), mid (M), bottom (B), and
"all" = one water-column value per cast (mean of S, M and B; complete casts
only). Default runs all four so surface can be compared with the full column.

Cleaning: fill values (-999/999) only. The 3x IQR rule from
clean_oceanography_data is NOT applied: it removes real extremes.

    python scripts/ch3_trends.py

Writes outputs/ch3_trends/trends.csv, annual_means.csv and one figure per
dataset and variable.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymannkendall as mk

from seasonal import decimal_year, long_term_trend, trend_per_year

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "outputs" / "ch3_trends"

FILL_VALUES = [-999, 999]
UNITS = {"temp": "C", "sal": "", "pco2": "uatm", "pH": "", "aragonite": "",
         "total alkalinity": "umol/kg", "Temp_E": "C", "Temp_W": "C"}

SHIP_START, SHIP_END = pd.Timestamp("2013-11-01"), pd.Timestamp("2023-08-31")
PERIODS = {
    "full record": (None, None),
    "2007-present": (pd.Timestamp("2007-01-01"), None),
    "ship window 2013-11 to 2023-08": (SHIP_START, SHIP_END),
}
MIN_MONTHS_LOGGER = 9      # a logger year counts if it has data in >= 9 months
MIN_SAMPLES_SHIP = 3       # a ship year counts if it has >= 3 surface samples

# Period line styles for the figures (one hue family, differentiated by dash)
PERIOD_STYLE = {"full record": dict(color="#1c5aa8", ls="-"),
                "2007-present": dict(color="#b8312f", ls="--"),
                "ship window 2013-11 to 2023-08": dict(color="#3d3c39", ls=":")}


SHIP_DEPTH_SETS = {
    "S": "Ship surface",
    "M": "Ship mid-depth",
    "B": "Ship bottom",
    "all": "Ship all depths (cast mean)",
}


def ship_depth_set(ship, cols, key):
    """
    One ship dataset for a depth choice.

    S / M / B: bottles from that depth only (S ~1.6 dbar, M ~10, B ~19 median).
    all: one value per cast = mean of its S, M and B bottles. Only casts with
    all three depths are used (84 of 92), so a cast with only a bottom bottle
    can't pull a year's mean towards bottom-water values, and one cast counts
    once rather than as three independent samples.
    """
    if key in ("S", "M", "B"):
        return ship[ship["Depth"] == key]
    depths = ship.groupby(["Station", "Date"])["Depth"].transform(lambda d: "".join(sorted(set(d))))
    complete = ship[depths == "BMS"]
    return (complete.groupby(["Station", "Date"], as_index=False)[cols].mean())


def load(ship_depths=("S", "all")):
    """Returns [(dataset name, DataFrame with Date + variables, variables, kind)]."""
    ship = pd.read_excel(RAW / "shipto2023.xlsx").rename(columns={
        "CTD temp": "temp", "CTD Sal": "sal", "pCO2": "pco2",
        "Alkalinity": "total alkalinity"})
    ship["Date"] = pd.to_datetime(ship["Date"], errors="coerce")
    ship_cols = ["temp", "sal", "pco2", "pH", "aragonite", "total alkalinity"]
    for c in ship_cols:
        ship[c] = pd.to_numeric(ship[c], errors="coerce").replace(FILL_VALUES, np.nan)
    ships = [(SHIP_DEPTH_SETS[k], ship_depth_set(ship, ship_cols, k), ship_cols, "ship")
             for k in ship_depths]

    h = pd.read_excel(RAW / "Temp E_W_1989-2024_1sheet.xlsx")
    loggers = []
    for name, dcol, tcol, var in [("HOBO East", "Date", "Etemp_reef cap SBE", "Temp_E"),
                                  ("HOBO West", "Date.1", "Wtemp_reef cap SBE", "Temp_W")]:
        # Each site has its own date column (they don't line up by row).
        d = h[[dcol, tcol]].dropna().rename(columns={dcol: "Date", tcol: var})
        d[var] = d[var].replace(FILL_VALUES, np.nan)
        # Fix: West has 2023-01-01..10 entered twice (identical values);
        # average to one value per day so those days aren't double-weighted.
        d = d.groupby("Date", as_index=False)[var].mean()
        loggers.append((name, d, [var], "logger"))
    return ships + loggers


def in_period(df, start, end):
    m = pd.Series(True, index=df.index)
    if start is not None:
        m &= df["Date"] >= start
    if end is not None:
        m &= df["Date"] <= end
    return df[m]


def annual_means(des, col, kind):
    """Annual means of the deseasoned series, keeping only complete-enough years."""
    g = des.groupby(des["Date"].dt.year)
    ann = pd.DataFrame({"mean": g[col].mean(), "n": g[col].count(),
                        "months": g["Date"].apply(lambda s: s.dt.month.nunique())})
    if kind == "logger":
        keep = ann["months"] >= MIN_MONTHS_LOGGER
    else:
        keep = ann["n"] >= MIN_SAMPLES_SHIP
    return ann, keep


def analyse(name, df, col, kind, period, start, end):
    d = in_period(df.dropna(subset=[col]).sort_values("Date"), start, end)
    if len(d) < 10:
        return None, None
    des, reg = long_term_trend(d, col)             # seasonal cycle refitted within period
    dcol = f"Deseasoned_{col}"
    ann, keep = annual_means(des, dcol, kind)
    used = ann[keep]
    row = {"dataset": name, "variable": col, "period": period,
           "first": d["Date"].min().date(), "last": d["Date"].max().date(),
           "n_samples": reg["n"],
           # reference only - not autocorrelation-corrected
           "sample_regression_slope_per_yr": reg["slope_per_yr"],
           "sample_regression_p_uncorrected": reg["p"],
           "years_used": len(used), "years_dropped_incomplete": int((~keep).sum())}
    if len(used) >= 4:
        r = mk.hamed_rao_modification_test(used["mean"].to_numpy())
        o = mk.original_test(used["mean"].to_numpy())
        # Fix: with ~10 annual means, Hamed-Rao sometimes estimates NEGATIVE
        # year-to-year autocorrelation and SHRINKS the variance, giving a
        # smaller p than the plain test (e.g. ship bottom pCO2: 7e-9 vs 3e-3;
        # bottom temp: 0.023 vs 0.16). The correction is meant to guard
        # against positive autocorrelation, so the reported p is the larger
        # (more cautious) of the two.
        p_final = float(np.nanmax([r.p, o.p]))   # Hamed-Rao can return NaN
        trend_final = ("no trend" if p_final >= 0.05 else
                       "increasing" if o.s > 0 else "decreasing")
        yrs = used.index.to_numpy(float)
        ts = pd.Series(used["mean"].to_numpy(), index=yrs)
        # Sen's slope on (year, value) pairs - correct even if years are missing
        i, j = np.triu_indices(len(yrs), k=1)
        slopes = (ts.values[j] - ts.values[i]) / (yrs[j] - yrs[i])
        sen_yr = float(np.median(slopes))
        ann_ols = trend_per_year(pd.to_datetime(used.index.astype(int).astype(str) + "-07-02"),
                                 used["mean"])
        row.update({"annual_sen_slope_per_yr": sen_yr,
                    "annual_mk_p": p_final, "annual_mk_trend": trend_final,
                    "annual_mk_original_p": o.p, "annual_mk_hamed_rao_p": r.p,
                    "annual_mk_tau": o.Tau,
                    "annual_ols_slope_per_yr": ann_ols["slope_per_yr"],
                    "annual_ols_stderr_per_yr": ann_ols["stderr_per_yr"],
                    "years_list": " ".join(str(y) for y in used.index)})
        if len(used) < 10:
            row["note"] = "fewer than 10 annual means: test has low power"
    else:
        row["note"] = "fewer than 4 complete years: no annual test"
    ann = ann.assign(dataset=name, variable=col, period=period, used=keep)
    return row, (des, ann)


def plot(name, col, runs, path):
    """Deseasoned samples (faint), annual means, and Sen's slope per period."""
    full_des, full_ann = runs["full record"][1]
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.scatter(full_des["Date"], full_des[f"Deseasoned_{col}"], s=4, color="#c9c8c4",
               linewidths=0, label="deseasoned samples", zorder=1)
    used = full_ann[full_ann["used"]]
    ax.scatter(pd.to_datetime(used.index.astype(str) + "-07-02"), used["mean"], s=36,
               color="#1f1f1d", zorder=3, label="annual mean (complete years)")
    for period, (row, (des, ann)) in runs.items():
        if "annual_sen_slope_per_yr" not in row:
            continue
        u = ann[ann["used"]]
        yrs = u.index.to_numpy(float)
        # Sen line through the median point (standard Sen intercept)
        b = np.median(u["mean"].to_numpy() - row["annual_sen_slope_per_yr"] * yrs)
        xs = np.array([yrs.min(), yrs.max()])
        p = row["annual_mk_p"]
        p_str = "p < 0.001" if p < 0.001 else f"p = {p:.3f}"
        ax.plot(pd.to_datetime([f"{int(x)}-07-02" for x in xs]),
                row["annual_sen_slope_per_yr"] * xs + b, lw=2, zorder=4,
                label=f"{period}: {row['annual_sen_slope_per_yr']:+.3f} "
                      f"{UNITS.get(col, '')}/yr ({p_str})", **PERIOD_STYLE[period])
    ax.set_ylabel(f"{col} (seasonal cycle removed)")
    ax.set_title(f"{name}: {col} - Sen's slope on annual means, Mann-Kendall p",
                 fontsize=11)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--ship-depths", nargs="+", default=["S", "M", "B", "all"],
                    choices=list(SHIP_DEPTH_SETS),
                    help="ship depth sets to analyse: S, M, B and/or all "
                         "(all = per-cast mean of S, M and B). Default: all four")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rows, anns = [], []
    for name, df, cols, kind in load(args.ship_depths):
        for col in cols:
            runs = {}
            for period, (start, end) in PERIODS.items():
                if kind == "ship" and period != "full record":
                    # ship data start Nov 2013, so every period equals the full record
                    continue
                row, extra = analyse(name, df, col, kind, period, start, end)
                if row is None:
                    continue
                rows.append(row)
                anns.append(extra[1].reset_index(names="year"))
                runs[period] = (row, extra)
            plot(name, col, runs,
                 OUT / f"trend_{name.replace(' ', '_')}_{col.replace(' ', '_')}.png")
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "trends.csv", index=False)
    pd.concat(anns).to_csv(OUT / "annual_means.csv", index=False)
    show = ["dataset", "variable", "period", "years_used", "annual_sen_slope_per_yr",
            "annual_mk_p", "annual_mk_trend",
            "sample_regression_slope_per_yr", "note"]
    with pd.option_context("display.float_format", "{:.4g}".format, "display.width", 250,
                           "display.max_colwidth", 50):
        print(res[[c for c in show if c in res]].to_string(index=False))
    print(f"\nWrote {OUT}/trends.csv, annual_means.csv and trend_*.png")


if __name__ == "__main__":
    main()
