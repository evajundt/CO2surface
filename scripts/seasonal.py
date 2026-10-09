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

  5. Trends are in units PER YEAR. The original plot_anomaly_with_trend and
     compute_seasonal_trend regressed on np.arange(n) (row number), giving
     "per sample" slopes that can't be converted by multiplying when
     sampling is irregular (HOBO West: x365.25 gives +0.070 C/yr; the true
     slope is +0.046 C/yr).
  6. Long-term trends come from long_term_trend(), which removes ONLY the
     seasonal cycle and then fits anomaly vs decimal year. Slopes fitted to
     Anomaly_ columns from the full GAM are ~0 by construction (p = 1.00),
     because the GAM's time spline has already absorbed the trend.

time_axis="rank" reproduces the original behaviour for comparison.
"""

import numpy as np
import pandas as pd
from pygam import LinearGAM, s
from scipy import stats


def decimal_year(dates):
    """Dates -> decimal year (e.g. 2015-07-02 -> 2015.50)."""
    d = pd.to_datetime(pd.Series(dates))
    return (d.dt.year + (d.dt.dayofyear - 1) / 365.25).to_numpy(float)


def trend_per_year(dates, values):
    """
    Change 5: least-squares slope of values against decimal year.
    Returns slope (units/yr), intercept, r2, p, stderr, and n.
    """
    t = decimal_year(dates)
    y = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)
    ok = np.isfinite(t) & np.isfinite(y)
    r = stats.linregress(t[ok], y[ok])
    return {"slope_per_yr": r.slope, "intercept": r.intercept, "r2": r.rvalue ** 2,
            "p": r.pvalue, "stderr_per_yr": r.stderr, "n": int(ok.sum()),
            "years": float(t[ok].max() - t[ok].min())}


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


def long_term_trend(df, col, date_col="Date", n_doy_splines=12):
    """
    Change 6: long-term trend with the seasonal cycle removed but the trend
    left in.

    Fits a GAM with only the cyclic day-of-year term (same 12-spline cyclic
    basis as compute_seasonal_trend_gam), takes Deseasoned_<col> = value -
    seasonal cycle, then fits Deseasoned_<col> vs decimal year.

    Returns (df with Deseasoned_<col> added, trend dict from trend_per_year).
    Use this for the per-year trends reported in the chapter.
    """
    out = df.copy()
    dates = pd.to_datetime(out[date_col])
    doy = dates.dt.dayofyear.to_numpy(float)
    y = pd.to_numeric(out[col], errors="coerce").to_numpy(float)
    ok = np.isfinite(y)
    gam = LinearGAM(s(0, basis="cp", n_splines=n_doy_splines, edge_knots=[1, 366]))
    gam.fit(doy[ok, None], y[ok])
    seas = np.full(len(out), np.nan)
    seas[ok] = gam.predict(doy[ok, None])
    out[f"Deseasoned_{col}"] = y - seas + np.nanmean(seas)   # keep original units/level
    return out, trend_per_year(dates, out[f"Deseasoned_{col}"])


def compute_seasonal_trend(df, col, date_col="Date", window=15):
    """
    Eva's day-of-year climatology method (same steps as the original), with
    the linear trend fitted against decimal year (Change 5), so Trend_<col>
    is in units per year. Adds Climatology-based SeasonalTrend_<col>,
    Anomaly_<col> and Trend_<col>; returns (df, trend dict).
    """
    out = df.copy()
    out[date_col] = pd.to_datetime(out[date_col])
    doy = out[date_col].dt.dayofyear
    clim = (out.groupby(doy)[col].mean().reindex(range(1, 367))
               .interpolate(method="linear"))
    pad = window
    padded = pd.concat([clim.iloc[-pad:], clim, clim.iloc[:pad]])
    smooth = padded.rolling(window=window, center=True, min_periods=1).mean().iloc[pad:-pad]
    smooth.index = range(1, 367)
    out[f"SeasonalTrend_{col}"] = doy.map(smooth)
    out[f"Anomaly_{col}"] = out[col] - out[f"SeasonalTrend_{col}"]
    tr = trend_per_year(out[date_col], out[f"Anomaly_{col}"])
    out[f"Trend_{col}"] = tr["slope_per_yr"] * decimal_year(out[date_col]) + tr["intercept"]
    return out, tr


def plot_anomaly_with_trend(df, value_col, date_col="Date", title=None, units=""):
    """
    Eva's plot_anomaly_with_trend with the slope in units per year
    (Change 5). Pass a Deseasoned_<col> column from long_term_trend (not an
    Anomaly_ column from the full GAM, whose trend is ~0 by construction).
    Returns the trend dict.
    """
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    d = df[[date_col, value_col]].dropna().sort_values(date_col)
    tr = trend_per_year(d[date_col], d[value_col])
    t = decimal_year(d[date_col])
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.scatter(d[date_col], d[value_col], color="steelblue", alpha=0.55, s=18,
               linewidths=0, label=value_col, zorder=2)
    ax.plot(d[date_col], tr["slope_per_yr"] * t + tr["intercept"], color="firebrick",
            linewidth=1.5, label="Line of best fit", zorder=3)
    p_str = "p < 0.001" if tr["p"] < 0.001 else f"p = {tr['p']:.3f}"
    ax.annotate(f"slope = {tr['slope_per_yr']:.3f} {units}/yr\n{p_str}\n$R^2$ = {tr['r2']:.3f}",
                xy=(0.98, 0.04), xycoords="axes fraction", ha="right", va="bottom",
                fontsize=9, color="firebrick",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="firebrick", alpha=0.8))
    ax.set_xlabel("Date")
    ax.set_ylabel(value_col)
    ax.set_title(title or f"{value_col} with long-term trend", fontsize=12, fontweight="bold")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=8))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.autofmt_xdate(rotation=30, ha="right")
    fig.tight_layout()
    return tr, fig
