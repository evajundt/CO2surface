"""
Per-year long-term trends for the chapter 3 datasets in the repo
(ship surface bottles, HOBO/SBE reef loggers).

Replaces the trend numbers from plot_anomaly_with_trend in
reference/fianlch3analysis.py, which were (a) per SAMPLE not per year, and
(b) ~0 by construction when run on GAM anomalies. See seasonal.py for details.

Same cleaning as the original pipeline: -999/999 fill values to NaN, then
3x IQR outliers to NaN (clean_oceanography_data's rule, without the plots).

    python scripts/ch3_trends.py

Writes outputs/ch3_trends/trends.csv and one trend figure per variable.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymannkendall as mk

from seasonal import iqr_outliers, long_term_trend, plot_anomaly_with_trend

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "outputs" / "ch3_trends"

UNITS = {"temp": "C", "sal": "", "pco2": "uatm", "pH": "", "aragonite": "",
         "total alkalinity": "umol/kg", "Temp_E": "C", "Temp_W": "C"}


def clean(df, cols):
    """clean_oceanography_data's fill-value + 3x IQR rule, no plots."""
    df = df.copy()
    n = {}
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").replace([-999, 999], np.nan)
        bad = iqr_outliers(df[c])
        n[c] = int(bad.sum())
        df.loc[bad, c] = np.nan
    return df, n


def load():
    ship = pd.read_excel(RAW / "shipto2023.xlsx").rename(columns={
        "CTD temp": "temp", "CTD Sal": "sal", "pCO2": "pco2",
        "Alkalinity": "total alkalinity"})
    ship["Date"] = pd.to_datetime(ship["Date"], errors="coerce")
    ship_cols = ["temp", "sal", "pco2", "pH", "aragonite", "total alkalinity"]
    ship, n_ship = clean(ship, ship_cols)
    ship = ship[ship["Depth"] == "S"]

    h = pd.read_excel(RAW / "Temp E_W_1989-2024_1sheet.xlsx")
    # East and West each have their own date column (they don't line up by row)
    east = h[["Date", "Etemp_reef cap SBE"]].dropna().rename(
        columns={"Etemp_reef cap SBE": "Temp_E"})
    west = h[["Date.1", "Wtemp_reef cap SBE"]].dropna().rename(
        columns={"Date.1": "Date", "Wtemp_reef cap SBE": "Temp_W"})
    east, n_e = clean(east, ["Temp_E"])
    west, n_w = clean(west, ["Temp_W"])
    return [("Ship surface", ship, ship_cols, n_ship),
            ("HOBO East", east, ["Temp_E"], n_e),
            ("HOBO West", west, ["Temp_W"], n_w)]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, df, cols, n_removed in load():
        for col in cols:
            d = df.dropna(subset=[col]).sort_values("Date")
            des, tr = long_term_trend(d, col)
            row = {"dataset": name, "variable": col, "n": tr["n"],
                   "years": round(tr["years"], 1),
                   "slope_per_yr": tr["slope_per_yr"], "stderr_per_yr": tr["stderr_per_yr"],
                   "p": tr["p"], "r2": tr["r2"], "iqr_removed": n_removed[col]}
            if name.startswith("HOBO"):
                # Daily values are strongly autocorrelated, so the ordinary
                # regression p-value is far too small. Annual means and the
                # Hamed-Rao autocorrelation-corrected Mann-Kendall test are
                # the defensible significance tests here.
                ann = des.groupby(des["Date"].dt.year)[f"Deseasoned_{col}"].agg(["mean", "size"])
                ann = ann[ann["size"] >= 300]["mean"]           # near-complete years only
                r = mk.hamed_rao_modification_test(ann.to_numpy())
                row.update({"annual_years_used": len(ann),
                            "annual_sen_slope_per_yr": r.slope, "annual_mk_hr_p": r.p})
            rows.append(row)
            _, fig = plot_anomaly_with_trend(
                des, f"Deseasoned_{col}", title=f"{name} - {col} (seasonal cycle removed)",
                units=UNITS.get(col, ""))
            fig.savefig(OUT / f"trend_{name.replace(' ', '_')}_{col.replace(' ', '_')}.png", dpi=150)
            plt.close(fig)
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "trends.csv", index=False)
    with pd.option_context("display.float_format", "{:.4g}".format, "display.width", 200):
        print(res.to_string(index=False))
    print(f"\nWrote {OUT}/trends.csv and trend_*.png")


if __name__ == "__main__":
    main()
