"""
Seasonal + trend removal, adapted from Eva's ch3 analysis
(reference/fianlch3analysis.py: compute_seasonal_trend_gam and
clean_oceanography_data).

Same model as the original - pyGAM LinearGAM with a cyclic 12-spline term on
day of year plus a 10-spline term on time - and the same output column names
(SeasonalTrend_<col>, Anomaly_<col>). Changes, each marked "Change:" below:

  1. Time axis is decimal year instead of row number. The original used
     np.arange(len(df)) after sorting, so time was the sample's rank: years
     with many samples got stretched and gaps (no data in 2019 or 2023)
     disappeared. With rank, the trend is in "per sample", not "per year".
  2. The cyclic day-of-year term wraps at 1/366 (Jan 1 / Dec 31). By default
     pyGAM wraps at the first/last day that has data (DOY 16 and 363 here).
  3. Rows are returned in their original order with NaN anomalies where the
     input was missing, instead of a sorted copy with those rows dropped, so
     the result can be assigned straight back to the input table.
  4. The IQR cleaning is applied to ANOMALIES, not raw values, and without
     plots. On raw SOCAT values 3x IQR removed 74 winter temperature visits
     (real Dec-Mar cold water) and 43 low-salinity visits (the plume water
     being tested), because the data are dominated by Aug/Sep.

time_axis="rank" reproduces the original behaviour for comparison.
"""

import numpy as np
import pandas as pd
from pygam import LinearGAM, s


def compute_seasonal_trend_gam(df, col, date_col="Date", time_axis="decimal_year",
                               n_doy_splines=12, n_time_splines=10):
    """
    Adds SeasonalTrend_<col> and Anomaly_<col> to a copy of df.

    Returns (df, info) where info holds the fit statistics and the implied
    long-term trend per year (from the time term of the GAM).
    """
    out = df.copy()
    dates = pd.to_datetime(out[date_col])
    doy = dates.dt.dayofyear.to_numpy(float)

    if time_axis == "decimal_year":
        # Change 1: real time, so the trend is per year and gaps are kept
        t = (dates.dt.year + (dates.dt.dayofyear - 1) / 365.25).to_numpy(float)
    elif time_axis == "rank":
        # Original behaviour: position in the date-sorted table
        t = np.empty(len(out))
        t[np.argsort(dates.to_numpy(), kind="stable")] = np.arange(len(out))
    else:
        raise ValueError("time_axis must be 'decimal_year' or 'rank'")

    y = pd.to_numeric(out[col], errors="coerce").to_numpy(float)
    ok = np.isfinite(y) & np.isfinite(t)
    X = np.column_stack([doy, t])

    gam = LinearGAM(
        # Change 2: wrap the seasonal cycle at the calendar year boundary
        s(0, basis="cp", n_splines=n_doy_splines, edge_knots=[1, 366]) +
        s(1, n_splines=n_time_splines))
    gam.fit(X[ok], y[ok])

    # Change 3: keep original row order; missing inputs stay NaN
    fit = np.full(len(out), np.nan)
    fit[ok] = gam.predict(X[ok])
    out[f"SeasonalTrend_{col}"] = fit
    out[f"Anomaly_{col}"] = y - fit

    info = {"r2": gam.statistics_["pseudo_r2"]["explained_deviance"],
            "n": int(ok.sum()), "time_axis": time_axis}
    if time_axis == "decimal_year":
        # Long-term trend = least-squares slope of the GAM's time term over
        # the observed times. (The difference between the curve's two end
        # points was unstable: spline ends wobble.)
        tt = t[ok]
        pdep = gam.partial_dependence(term=1, X=np.column_stack([np.full(ok.sum(), 182.0), tt]))
        info["trend_per_yr"] = float(np.polyfit(tt, pdep, 1)[0])
    # Seasonal amplitude from the cyclic term alone (peak to trough / 2)
    grid = np.column_stack([np.arange(1, 366), np.full(365, np.nanmedian(t[ok]))])
    seas = gam.partial_dependence(term=0, X=grid)
    info["seasonal_amplitude"] = (seas.max() - seas.min()) / 2
    # Residual seasonality check (Panel 4 of plot_seasonal_diagnostics):
    # mean anomaly by month should be ~0 if the cycle was removed
    by_month = pd.Series(out[f"Anomaly_{col}"].to_numpy()).groupby(dates.dt.month.to_numpy()).mean()
    info["max_abs_monthly_mean_anomaly"] = float(by_month.abs().max())
    info["monthly_mean_anomaly"] = by_month.round(2).to_dict()
    return out, info


def iqr_outliers(values, iqr_multiplier=3.0):
    """
    Boolean mask of IQR outliers - the same rule as clean_oceanography_data
    (Q1 - k*IQR, Q3 + k*IQR, default k = 3), without the plots.
    Change 4: meant to be applied to anomalies, not raw values.
    """
    v = pd.Series(values, dtype=float)
    q1, q3 = v.quantile([0.25, 0.75])
    iqr = q3 - q1
    return ((v < q1 - iqr_multiplier * iqr) | (v > q3 + iqr_multiplier * iqr)).to_numpy()
