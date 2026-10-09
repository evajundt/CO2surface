# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 13:18:15 2026

@author: ejundt
"""

# -*- coding: utf-8 -*-
"""
Seasonal/Trend/Interpolation/Stats Analysis for Shipboard/Discrete and Surface Data (East, West, Stetson)
@author: you
"""
import os
print(os.getcwd())
os.chdir('C:/Users/ejundt/OneDrive - Texas A&M University-Corpus Christi/Desktop/March 2026/ch3 work/code and data')
print(os.getcwd())
# %% Imports and Functions
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import SplineTransformer
import pymannkendall as mk
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler as stsc
from sklearn.preprocessing import PolynomialFeatures
# CHANGED (Claude): mean_squared_error(..., squared=False) was removed in scikit-learn 1.6;
# root_mean_squared_error (scikit-learn >= 1.4) gives the same RMSE.
from sklearn.metrics import mean_absolute_error, root_mean_squared_error, r2_score
from sklearn.ensemble import RandomForestRegressor
from scipy import stats
import matplotlib.ticker as ticker
from pygam import LinearGAM, s

def compute_seasonal_trend_gam(df, col, date_col='Date', window=15):
    dfog = df.copy()
    dfog[date_col] = pd.to_datetime(dfog[date_col])
    dfog = dfog.sort_values(date_col).reset_index(drop=True)
    dfog['_doy'] = dfog[date_col].dt.dayofyear
    # NOTE (Claude): time is the row number after sorting, not real time, so
    # years with many samples are stretched and gaps (e.g. 2019, 2023) vanish.
    # The time spline also absorbs the long-term trend, so Anomaly_<col> has
    # no trend left in it. Each call also drops rows where <col> is missing,
    # so chained calls keep only rows complete in every variable.
    # Fixed version: scripts/seasonal.py.
    dfog['_time_index'] = np.arange(len(dfog))
    # Remove NaN and Inf before fitting
    valid_mask = dfog[col].notna() & np.isfinite(dfog[col])
    df = dfog[valid_mask].copy()

    X = df[['_doy', '_time_index']].values
    y = df[col].values

    gam = LinearGAM(s(0, basis='cp', n_splines=12) +
                    s(1, n_splines=10))
    gam.fit(X, y)

    df[f'SeasonalTrend_{col}'] = gam.predict(X)
    df[f'Anomaly_{col}'] = y - df[f'SeasonalTrend_{col}']

    df.drop(columns=['_doy', '_time_index'], inplace=True)
    return df  

def compute_seasonal_trend(df, col, date_col='Date', window=15):
    """
    Computes a climatological seasonal cycle and anomaly for oceanographic data.

    Method:
        1. Groups data by day-of-year and smooths with a rolling window
           to build a smooth climatological mean (seasonal cycle).
        2. Anomaly = observed value - climatological mean for that day-of-year.
        3. A long-term linear trend is fit to the anomaly series.

    Parameters:
        df       : DataFrame with a datetime date column
        col      : Name of the variable column to process
        date_col : Name of the datetime column (default 'Date')
        window   : Smoothing window in days applied to the climatology (default 30)

    Returns:
        df with three new columns added:
            Climatology_{col}  : Smoothed seasonal cycle
            Anomaly_{col}      : Deseasonalised anomaly
            Trend_{col}        : Long-term linear trend fitted to the anomaly
    """
    df = df.copy()

    # Ensure datetime
    df[date_col] = pd.to_datetime(df[date_col])
    df['_doy'] = df[date_col].dt.dayofyear  # 1–366

    # ------------------------------------------------------------------ #
    # 1. BUILD CLIMATOLOGICAL SEASONAL CYCLE
    #    Mean value for each day-of-year across all years,
    #    then smooth with a centred rolling window to avoid sharp jumps.
    # ------------------------------------------------------------------ #
    clim = (
        df.groupby('_doy')[col]
        .mean()
        .reindex(range(1, 367))          # ensure all 366 doys present
        .interpolate(method='linear')    # fill any missing doys
    )

    # Wrap-around smoothing: pad edges with the other end of the year
    # so the window doesn't see a hard boundary at doy 1/366
    pad = window
    clim_padded = pd.concat([clim.iloc[-pad:], clim, clim.iloc[:pad]])
    clim_smooth = (
        clim_padded
        .rolling(window=window, center=True, min_periods=1)
        .mean()
        .iloc[pad:-pad]                  # trim padding back off
    )
    clim_smooth.index = range(1, 367)    # restore original doy index

    # ------------------------------------------------------------------ #
    # 2. MAP CLIMATOLOGY BACK TO EACH ROW & COMPUTE ANOMALY
    # ------------------------------------------------------------------ #
    df[f'SeasonalTrend_{col}'] = df['_doy'].map(clim_smooth)
    df[f'Anomaly_{col}']     = df[col] - df[f'SeasonalTrend_{col}']

    # ------------------------------------------------------------------ #
    # 3. FIT LINEAR TREND TO ANOMALY SERIES
    # ------------------------------------------------------------------ #
    valid = df[f'Anomaly_{col}'].notna()
    # NOTE (Claude): x = np.arange(...) is the ROW NUMBER, so this slope is
    # "per sample", not per year. With irregular sampling it cannot be
    # converted by multiplying (HOBO West: x365.25 gives +0.070 C/yr, the
    # true slope is +0.046 C/yr). Per-year versions: scripts/seasonal.py.
    x = np.arange(len(df))
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        x[valid], df.loc[valid, f'Anomaly_{col}']
    )
    df[f'Trend_{col}'] = slope * x + intercept

    # Print trend stats
    p_str = "p < 0.001" if p_value < 0.001 else f"p = {p_value:.3f}"
    print(f"[{col}] Trend fit:")
    print(f"  Slope:   {slope:.6f}  (units/timestep)")
    print(f"  {p_str}")
    print(f"  R²:      {r_value**2:.4f}")
    print(f"  Std err: {std_err:.6f}\n")

    df.drop(columns='_doy', inplace=True)
    return df

def clean_oceanography_data(df, data_cols=None, date_col=None, iqr_multiplier=3.0):
    """
    Oceanography data cleaning pipeline.
    
    Parameters:
        df             : Raw DataFrame (e.g. freshly loaded from spreadsheet)
        data_cols      : List of columns to clean. If None, auto-detects numeric columns.
        date_col       : Optional date column name for x-axis on plots.
        iqr_multiplier : Sensitivity for outlier detection (default 3.0; lower = stricter).
    
    Returns:
        Cleaned DataFrame with row IDs added.
    """

    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 11,
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'figure.dpi': 150,
    })

    # ------------------------------------------------------------------ #
    # 1. ADD ROW IDs
    # ------------------------------------------------------------------ #
    df = df.copy()
    df.insert(0, 'row_id', range(1, len(df) + 1))
    print(f"✔ Added row IDs (1 – {len(df)})\n")

    # Auto-detect numeric columns if not specified
    if data_cols is None:
        data_cols = [c for c in df.select_dtypes(include=[np.number]).columns if c != 'row_id']
        print(f"Auto-detected numeric columns: {data_cols}\n")

    x_axis = df[date_col] if date_col else df['row_id']
    x_label = date_col if date_col else 'Row ID'

    # ------------------------------------------------------------------ #
    # 2. PLOT RAW DATA
    # ------------------------------------------------------------------ #
    n_cols = len(data_cols)
    fig, axes = plt.subplots(n_cols, 1, figsize=(10, 3.5 * n_cols), sharex=True)
    if n_cols == 1:
        axes = [axes]

    for ax, col in zip(axes, data_cols):
        ax.scatter(x_axis, df[col], s=10, color='steelblue', alpha=0.5, linewidths=0)
        ax.axhline(0, color='black', linestyle='--', linewidth=0.6, alpha=0.4)
        ax.set_ylabel(col, labelpad=6)
        ax.set_title(f'{col} — Raw', fontsize=10, fontweight='bold')

    axes[-1].set_xlabel(x_label, labelpad=6)
    axes[-1].xaxis.set_major_locator(ticker.MaxNLocator(nbins=10))
    fig.autofmt_xdate(rotation=30, ha='right')
    fig.suptitle('Raw Data', fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.show()

    # ------------------------------------------------------------------ #
    # 3. REMOVE FILL VALUES (-999 / 999)
    # ------------------------------------------------------------------ #
    fill_values = [-999, 999]
    fill_removed = {}

    for col in data_cols:
        mask = df[col].isin(fill_values)
        count = mask.sum()
        if count > 0:
            fill_removed[col] = df.loc[mask, ['row_id', col]].to_dict('records')
            df.loc[mask, col] = np.nan
        print(f"  [{col}] Fill values removed: {count}")

    print()
    if fill_removed:
        print("Fill value locations:")
        for col, rows in fill_removed.items():
            for r in rows:
                print(f"    Row {r['row_id']:>6}  |  {col}: {r[col]}")
        print()

    # ------------------------------------------------------------------ #
    # 4. DETECT & REMOVE EXTREME OUTLIERS (IQR METHOD)
    # ------------------------------------------------------------------ #
    outlier_log = {}

    for col in data_cols:
        clean = df[col].dropna()
        Q1 = clean.quantile(0.25)
        Q3 = clean.quantile(0.75)
        IQR = Q3 - Q1
        lower = Q1 - iqr_multiplier * IQR
        upper = Q3 + iqr_multiplier * IQR

        mask = (df[col] < lower) | (df[col] > upper)
        outliers = df.loc[mask, ['row_id', col]].dropna()

        if not outliers.empty:
            outlier_log[col] = outliers
            df.loc[mask, col] = np.nan

        print(f"  [{col}]  IQR bounds: [{lower:.4f}, {upper:.4f}]  |  Outliers removed: {len(outliers)}")

    print()
    if outlier_log:
        print("Outlier details:")
        for col, rows in outlier_log.items():
            for _, r in rows.iterrows():
                print(f"    Row {int(r['row_id']):>6}  |  {col}: {r[col]:.4f}")
        print()
    else:
        print("  No extreme outliers detected.\n")

    # ------------------------------------------------------------------ #
    # 5. PLOT CLEANED DATA
    # ------------------------------------------------------------------ #
    fig, axes = plt.subplots(n_cols, 1, figsize=(10, 3.5 * n_cols), sharex=True)
    if n_cols == 1:
        axes = [axes]

    for ax, col in zip(axes, data_cols):
        ax.scatter(x_axis, df[col], s=10, color='seagreen', alpha=0.5, linewidths=0)
        ax.axhline(0, color='black', linestyle='--', linewidth=0.6, alpha=0.4)
        ax.set_ylabel(col, labelpad=6)
        ax.set_title(f'{col} — Cleaned', fontsize=10, fontweight='bold')

        # Annotate total points removed per column
        n_removed = df[col].isna().sum()
        ax.annotate(
            f'{n_removed} points removed',
            xy=(0.99, 0.97), xycoords='axes fraction',
            ha='right', va='top', fontsize=8.5,
            color='dimgrey',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='grey', alpha=0.7)
        )

    axes[-1].set_xlabel(x_label, labelpad=6)
    axes[-1].xaxis.set_major_locator(ticker.MaxNLocator(nbins=10))
    fig.autofmt_xdate(rotation=30, ha='right')
    fig.suptitle('Cleaned Data', fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.show()

    # ------------------------------------------------------------------ #
    # 6. SUMMARY
    # ------------------------------------------------------------------ #
    print("=" * 50)
    print("CLEANING SUMMARY")
    print("=" * 50)
    total_fill = sum(len(v) for v in fill_removed.values())
    total_outliers = sum(len(v) for v in outlier_log.values())
    print(f"  Rows in dataset:       {len(df)}")
    print(f"  Fill values removed:   {total_fill}")
    print(f"  Outliers removed:      {total_outliers}")
    print(f"  Total points cleaned:  {total_fill + total_outliers}")
    print("=" * 50)

    return df
def spline_model_analysis(df, xvars, yvar, datevar, degree=3, n_knots=5, label=None):
    '''
    xvars: list of 1-3 predictor column names, e.g., ['temp'] or ['temp', 'sal']
    '''
    dft = df[xvars + [yvar, datevar]].dropna()
    X = dft[xvars].values
    y = dft[yvar].values
    spline = SplineTransformer(degree=degree, n_knots=n_knots, include_bias=False)
    X_spline = spline.fit_transform(X)
    model = LinearRegression().fit(X_spline, y)
    y_pred = model.predict(X_spline)

    print("-"*30)
    print(f"Spline Regression (deg={degree}, knots={n_knots}) for {yvar} ~ {xvars}")
    print(f"Number of samples: {len(dft)}")
    print(f"R²: {r2_score(y, y_pred):.3f}")
    print(f"MAE: {mean_absolute_error(y, y_pred):.3f}")
    print(f"RMSE: {root_mean_squared_error(y, y_pred):.3f}")  # CHANGED (Claude): was mean_squared_error(..., squared=False)
    print("Coefficients:", model.coef_)
    print(f"Intercept: {model.intercept_:.3f}")
    print("-"*30)
    dft = dft.copy()
    dft['model_pred'] = y_pred

    # X/Y plot (first predictor only)
    plt.figure(figsize=(8,5))
    plt.scatter(dft[xvars[0]], dft[yvar], color='steelblue', alpha=0.6, label='Observed')
    xfit = np.linspace(dft[xvars[0]].min(), dft[xvars[0]].max(), 200)
    # Fill non-primary vars with median if >1 predictors
    X_predict = np.column_stack([xfit if i==0 else np.full_like(xfit, dft[x].median()) for i,x in enumerate(xvars)])
    yfit = model.predict(spline.transform(X_predict))
    plt.plot(xfit, yfit, color='darkorange', lw=2, label='Spline Fit')
    plt.xlabel(xvars[0])
    plt.ylabel(yvar)
    plt.title(f'Spline Regression: {yvar} ~ {xvars}')
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Time series
    plt.figure(figsize=(12,5))
    plt.plot(dft[datevar], dft[yvar], label='Observed', marker='o', linestyle='-', alpha=0.7)
    plt.plot(dft[datevar], dft['model_pred'], label='Predicted', marker='.', linestyle='--', color='red', alpha=0.7)
    plt.xlabel('Date')
    plt.ylabel(yvar)
    plt.title(f'Time Series: Measured vs. Spline Model {yvar}')
    plt.legend()
    plt.tight_layout()
    plt.show()

    return model, spline
def rf_model_analysis(df, xvars, yvar, datevar, n_estimators=100, max_depth=None, random_state=42, label=None):
    '''
    xvars: list of 1-3 predictor column names, e.g., ['temp'] or ['temp', 'sal']
    '''
    dft = df[xvars + [yvar, datevar]].dropna()
    X = dft[xvars].values
    y = dft[yvar].values
    model = RandomForestRegressor(n_estimators=n_estimators, max_depth=max_depth, random_state=random_state)
    model.fit(X, y)
    y_pred = model.predict(X)

    print("-"*30)
    print(f"Random Forest Regression for {yvar} ~ {xvars}")
    print(f"Number of samples: {len(dft)}")
    print(f"R²: {r2_score(y, y_pred):.3f}")
    print(f"MAE: {mean_absolute_error(y, y_pred):.3f}")
    print(f"RMSE: {root_mean_squared_error(y, y_pred):.3f}")  # CHANGED (Claude): was mean_squared_error(..., squared=False)
    print("Feature importances:", model.feature_importances_)
    print("-"*30)
    dft = dft.copy()
    dft['model_pred'] = y_pred

    # X/Y plot (first predictor only)
    plt.figure(figsize=(8,5))
    plt.scatter(dft[xvars[0]], dft[yvar], color='steelblue', alpha=0.6, label='Observed')
    xfit = np.linspace(dft[xvars[0]].min(), dft[xvars[0]].max(), 200)
    # Fill non-primary vars with median if >1 predictors
    X_predict = np.column_stack([xfit if i==0 else np.full_like(xfit, dft[x].median()) for i,x in enumerate(xvars)])
    yfit = model.predict(X_predict)
    plt.plot(xfit, yfit, color='green', lw=2, label='RF Fit')
    plt.xlabel(xvars[0])
    plt.ylabel(yvar)
    plt.title(f'Random Forest Regression: {yvar} ~ {xvars}')
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Time series
    plt.figure(figsize=(12,5))
    plt.plot(dft[datevar], dft[yvar], label='Observed', marker='o', linestyle='-', alpha=0.7)
    plt.plot(dft[datevar], dft['model_pred'], label='Predicted', marker='.', linestyle='--', color='red', alpha=0.7)
    plt.xlabel('Date')
    plt.ylabel(yvar)
    plt.title(f'Time Series: Measured vs. Random Forest Model {yvar}')
    plt.legend()
    plt.tight_layout()
    plt.show()

    return model
def poly_model_analysis(df, xvar, yvar, datevar, degree=2, label=None):
    '''
    Fit polynomial regression, print stats, and plot results.
    df: DataFrame
    xvar, yvar, datevar: column names
    degree: polynomial degree (default 2)
    label: optional string for legend
    '''
    # Drop missing date
    dft = df[[xvar, yvar, datevar]].dropna()
    X = dft[[xvar]].values
    y = dft[yvar].values

    # Fit model
    poly = PolynomialFeatures(degree=degree)
    Xpoly = poly.fit_transform(X)
    model = LinearRegression().fit(Xpoly, y)
    y_pred = model.predict(Xpoly)

    # METRICS
    print("-"*30)
    print(f"Polynomial Regression ({degree=}) for {yvar} ~ {xvar}")
    print(f"Number of samples: {len(dft)}")
    print(f"R²: {r2_score(y, y_pred):.3f}")
    print(f"MAE: {mean_absolute_error(y, y_pred):.3f}")
    print(f"RMSE: {root_mean_squared_error(y, y_pred):.3f}")  # CHANGED (Claude): was mean_squared_error(..., squared=False)
    print("Coefficients:", model.coef_)
    print(f"Intercept: {model.intercept_:.3f}")
    print("-"*30)
    # Add predictions to DataFrame for plotting
    dft = dft.copy()
    dft['model_pred'] = y_pred
    # Visual 1: Scatter with regression fit
    plt.figure(figsize=(8,5))
    plt.scatter(dft[xvar], dft[yvar], color='steelblue', alpha=0.6, label='Observed')
    # Curve for model fit
    xfit = np.linspace(dft[xvar].min(), dft[xvar].max(), 200).reshape(-1,1)
    yfit = model.predict(poly.transform(xfit))
    plt.plot(xfit, yfit, color='darkorange', lw=2, label='Model Fit')
    plt.xlabel(xvar)
    plt.ylabel(yvar)
    plt.title(f'Polynomial Regression: {yvar} ~ {xvar}')
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Visual 2: Time series for measured and predicted
    plt.figure(figsize=(12,5))
    plt.plot(dft[datevar], dft[yvar], label='Observed', marker='o', linestyle='-', alpha=0.7)
    plt.plot(dft[datevar], dft['model_pred'], label='Predicted', marker='.', linestyle='--', color='red', alpha=0.7)
    plt.xlabel('Date')
    plt.ylabel(yvar)
    plt.title(f'Time Series: Measured vs. Modelled {yvar}')
    plt.legend()
    plt.tight_layout()
    plt.show()

    return model, poly
def compute_seasonal_trend_og(df, col, window=60):
    trend_col = f'SeasonalTrend_{col}'
    df[trend_col] = df[col].rolling(window=window, center=True, min_periods=1).mean()
    anomaly_col = f'Anomaly_{col}'
    df[anomaly_col] = df[col] - df[trend_col]
    return df

def retimeyeardaytrend_HOBO(dates, temps, name):
    df = pd.DataFrame({'Date': dates, 'Temp': temps})
    df = df.dropna()
    df = df.set_index('Date').resample('D').mean().reset_index()
    df['Year'] = df['Date'].dt.year
    df['Yearday'] = df['Date'].dt.dayofyear
    trend_col = f'SeasonalTrend{name}'
    df[trend_col] = df['Temp'].rolling(window=60, center=True, min_periods=1).mean()
    anomaly_col = f'Anomaly{name}'
    df[anomaly_col] = df['Temp'] - df[trend_col]
    return df    
def compute_spline_trend(df, col, knots=6, degree=3):
    df = df.dropna(subset=[col])
    X = df['Date'].dt.dayofyear.values.reshape(-1, 1)
    y = df[col].values
    spline = SplineTransformer(n_knots=knots, degree=degree, include_bias=False)
    X_spline = spline.fit_transform(X)
    model = LinearRegression()
    model.fit(X_spline, y)
    trend_col = f'SplineTrend_{col}'
    df[trend_col] = model.predict(X_spline)
    anomaly_col = f'AnomalySpline_{col}'
    df[anomaly_col] = y - df[trend_col]
    return df
def print_mann_kendall(df, col):
    result = mk.original_test(df[col].dropna())
    print(f"{col}: Trend = {result.trend}, p = {result.p}, S = {result.s}")
def standardize_columns(df, rename_map):
    return df.rename(columns=rename_map)
def plot_anomaly_with_trend(df, anomaly_var, date_col='Date', title=None):
    
    
    # Drop missing values
    temp_df = df[[date_col, anomaly_var]].dropna()
    # NOTE (Claude): x = np.arange(...) is the ROW NUMBER, so this slope is
    # "per sample", not per year. With irregular sampling it cannot be
    # converted by multiplying (HOBO West: x365.25 gives +0.070 C/yr, the
    # true slope is +0.046 C/yr). Per-year versions: scripts/seasonal.py.
    # ALSO: on Anomaly_ columns from compute_seasonal_trend_gam the slope is
    # ~0 by construction (p = 1.00), because the GAM's time spline already
    # removed the trend. Use seasonal.long_term_trend for the chapter trends.
    x = np.arange(len(temp_df))
    y = temp_df[anomaly_var].values

    # --- Stats ---
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    y_trend = slope * x + intercept

    # Print stats
    print(f"Slope:   {slope:.4f}")
    print(f"P-value: {p_value:.4e}")
    print(f"R²:      {r_value**2:.4f}")

    # --- Publication-quality styling ---
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 11,
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'figure.dpi': 150,
    })

    fig, ax = plt.subplots(figsize=(8, 4.5))

    # Scatter
    ax.scatter(
        temp_df[date_col], y,
        color='steelblue', alpha=0.55, s=18,
        linewidths=0, label='Anomaly', zorder=2
    )

    # Line of best fit
    line, = ax.plot(
        temp_df[date_col], y_trend,
        color='firebrick', linewidth=1.5,
        label='Line of Best Fit', zorder=3
    )

    # --- Annotate line with slope & p-value ---
    # Place label near the right end of the trend line
    label_x = temp_df[date_col].iloc[-1]
    label_y = y_trend[-1]

    # Format p-value: show exact value or <0.001
    p_str = "p < 0.001" if p_value < 0.001 else f"p = {p_value:.3f}"
    annotation = f"slope = {slope:.3f}\n{p_str}\n$R^2$ = {r_value**2:.3f}"

    ax.annotate(
        annotation,
        xy=(label_x, label_y),
        xytext=(-10, 12),                      # offset in points
        textcoords='offset points',
        ha='right', va='bottom',
        fontsize=9,
        color='firebrick',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='firebrick', alpha=0.8),
    )

    # Zero line
    ax.axhline(0, color='black', linestyle='--', linewidth=0.7, alpha=0.5, zorder=1)

    # Labels & title
    ax.set_xlabel('Date', labelpad=6)
    ax.set_ylabel(anomaly_var, labelpad=6)
    ax.set_title(
        title if title else f'{anomaly_var} Anomaly with Line of Best Fit',
        fontsize=12, fontweight='bold', pad=10
    )

    # Legend
    ax.legend(frameon=False, fontsize=9, loc='upper left')

    # Tidy x-axis
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=8))
    fig.autofmt_xdate(rotation=30, ha='right')

    plt.tight_layout()
    plt.show()
# %% Load Data and Apply Column Standardization for shipdata
rename_map = {
    'sst': 'temp',
    'Etemp_reef cap SBE': 'temp',
    'Wtemp_reef cap SBE': 'temp',
    'CTD temp': 'temp',
    'sss': 'sal',
    'CTD Sal': 'sal',
    'fCO2rec': 'pco2',
    'pCO2 (µatm)': 'pco2',
    'pCO2': 'pco2','latitude': 'Latitude','longitude': 'Longitude','Alkalinity' : 'total alkalinity'
}

ship = pd.read_excel('shipto2023.xlsx')
ship['Date'] = pd.to_datetime(ship['Date'], errors='coerce')
ship = standardize_columns(ship, rename_map)
ship = clean_oceanography_data(ship, data_cols=['temp', 'sal', 'pH', 'DIC', 'pco2'], date_col='Date')
#ship = ship.dropna(subset=['pco2'])
# %% Load E W S
east = pd.read_csv('EFGB_pixel_0.10deg.csv')
east['Date'] = pd.to_datetime(east[['yr','mon','day','hh','mm']].rename(
    columns={'yr':'year','mon':'month','day':'day','hh':'hour','mm':'minute'}), errors='coerce')
east['Station'] = 'East'
east = standardize_columns(east, rename_map)
east = clean_oceanography_data(east, data_cols=['temp', 'sal', 'pco2'], date_col='Date')
west = pd.read_csv('WFGB_pixel_0.10deg.csv')
west['Date'] = pd.to_datetime(west[['yr','mon','day','hh','mm']].rename(
    columns={'yr':'year','mon':'month','day':'day','hh':'hour','mm':'minute'}), errors='coerce')
west['Station'] = 'West'
west = standardize_columns(west, rename_map)
west = clean_oceanography_data(west, data_cols=['temp', 'sal', 'pco2'], date_col='Date')
stetson = pd.read_csv('Stetson_pixel_0.10deg.csv')
stetson['Date'] = pd.to_datetime(stetson[['yr','mon','day','hh','mm']].rename(
    columns={'yr':'year','mon':'month','day':'day','hh':'hour','mm':'minute'}), errors='coerce')
stetson['Station'] = 'Stetson'
stetson = standardize_columns(stetson, rename_map)
stetson = clean_oceanography_data(stetson, data_cols=['temp', 'sal', 'pco2'], date_col='Date')
# %% Load HOBO
data = pd.read_excel('Temp E_W_1989-2024_1sheet.xlsx')
data['Date'] = pd.to_datetime(data['Date'], errors='coerce')
data['Date.1'] = pd.to_datetime(data['Date.1'], errors='coerce')
hobo_east = data[['Date', 'Etemp_reef cap SBE']].dropna()
hobo_west = data[['Date.1', 'Wtemp_reef cap SBE']].dropna()
hobo_east = hobo_east.rename(columns={'Date': 'Datetime', 'Etemp_reef cap SBE': 'Temp_E'})
hobo_west = hobo_west.rename(columns={'Date.1': 'Datetime', 'Wtemp_reef cap SBE': 'Temp_W'})
hobo_west['Date']=hobo_west['Datetime']
hobo_east['Date']=hobo_east['Datetime']
hobo_east = clean_oceanography_data(hobo_east, data_cols=['Temp_E'], date_col='Date')
hobo_west = clean_oceanography_data(hobo_west, data_cols=['Temp_W'], date_col='Date')
# %% Load Data and Apply Column Standardization for shipdata
# Ensure 'Date' is datetime type
for col in [ship, east, west, stetson]:
    print("Oldest datapoint:")
    print(col.loc[col['Date'].idxmin()])
    print("\nNewest datapoint:")
    print(col.loc[col['Date'].idxmax()])
    print("\nColumn names:", list(col.columns))
    print("\nNumber of datapoints:", len(col))
    print("\nLatitude: min =", col['Latitude'].min(), "max =", col['Latitude'].max())
    print("Longitude: min =", col['Longitude'].min(), "max =", col['Longitude'].max())
for col in [hobo_east, hobo_west]:
    print("Oldest datapoint:")
    print(col.loc[col['Datetime'].idxmin()])
    print("\nNewest datapoint:")
    print(col.loc[col['Datetime'].idxmax()])
    print("\nColumn names:", list(col.columns))
    print("\nNumber of datapoints:", len(col))
# %%  Apply Trend Calculations to Each Dataset
# UNCOMMENTED (Claude): these create the SeasonalTrend_*/Anomaly_* columns
# that the Mann-Kendall loop and the East/West SOCAT plots below need;
# without them the script stops with a KeyError.
# NOTE: as written, East temp uses the GAM and everything else here uses the
# day-of-year climatology (compute_seasonal_trend), while ship_surface uses
# the GAM for all variables. Worth making these consistent for the chapter.
east = compute_seasonal_trend_gam(east, 'temp')
east = compute_seasonal_trend(east, 'sal')
east = compute_seasonal_trend(east, 'pco2')

west = compute_seasonal_trend(west, 'temp')
west = compute_seasonal_trend(west, 'sal')
west = compute_seasonal_trend(west, 'pco2')

stetson = compute_seasonal_trend(stetson, 'temp')
stetson = compute_seasonal_trend(stetson, 'sal')
stetson = compute_seasonal_trend(stetson, 'pco2')

# Shipboard surface
ship_surface = ship[ship['Depth'] == 'S'].copy()
ship_surface = compute_seasonal_trend_gam(ship_surface, 'temp')
ship_surface = compute_seasonal_trend_gam(ship_surface, 'sal')
ship_surface = compute_seasonal_trend_gam(ship_surface, 'pco2')
ship_surface = compute_seasonal_trend_gam(ship_surface, 'pH')
ship_surface = compute_seasonal_trend_gam(ship_surface, 'aragonite')
ship_surface = compute_seasonal_trend_gam(ship_surface, 'total alkalinity')

# # HOBO logger trends (example for East)
# UNCOMMENTED (Claude): df_east is used by the Mann-Kendall loop and the
# "East Bank HOBO" plot below; without this line the script stops with a
# NameError. NOTE: East uses the climatology method, West uses the GAM.
df_east = compute_seasonal_trend(hobo_east, 'Temp_E')
df_west = compute_seasonal_trend_gam(hobo_west, 'Temp_W')

# %%  Apply Trend Calculations to Each Dataset
for col in ['SeasonalTrend_temp', 'Anomaly_temp', 'SeasonalTrend_sal', 'Anomaly_sal', 
            'SeasonalTrend_pco2', 'Anomaly_pco2']:
    print('East Bank SOCAT')
    print_mann_kendall(east, col)
    print('West Bank SOCAT')
    print_mann_kendall(west, col)
    print('Stetson Bank SOCAT')
    print_mann_kendall(stetson, col)
    print('Ship Surface')
    print_mann_kendall(ship_surface, col)
for col in ['Temp_W', 'SeasonalTrend_Temp_W', 'Anomaly_Temp_W']:    
    print('West Bank HOBO')
    print_mann_kendall(df_west, col)  
for col in ['Temp_E', 'SeasonalTrend_Temp_E', 'Anomaly_Temp_E']:
    print('East Bank HOBO')
    print_mann_kendall(df_east, col)    
    
# %%  Train polynomial model: pCO2 ~ 
# model, poly = poly_model_analysis(ship_surface, 'temp', 'pco2', datevar='Date')
# model, poly = poly_model_analysis(ship_surface, 'temp', 'pH', datevar='Date')
# model, poly = poly_model_analysis(ship_surface, 'temp', 'aragonite', datevar='Date')
# model, poly = poly_model_analysis(ship_surface, 'temp', 'sal', datevar='Date')
# model, poly = poly_model_analysis(ship_surface, 'temp', 'DIC', datevar='Date')
# model, poly = poly_model_analysis(ship_surface, 'SeasonalTrend_temp', 'pco2', datevar='Date')

# model, poly = poly_model_analysis(east, 'temp', 'pco2', datevar='Date')
# model, poly = poly_model_analysis(east, 'SeasonalTrend_temp', 'pco2', datevar='Date')
# model, poly = poly_model_analysis(east, 'sal', 'pco2', datevar='Date')
# model, poly = poly_model_analysis(west, 'temp', 'pco2', datevar='Date')
# model, poly = poly_model_analysis(west, 'SeasonalTrend_temp', 'pco2', datevar='Date')
# model, poly = poly_model_analysis(west, 'sal', 'pco2', datevar='Date')

# %%  Plot the longterm trends


# plot_anomaly_with_trend(ship_surface, 'Anomaly_temp', date_col='Date',title='Trends in Surface Ship Measurements')
# plot_anomaly_with_trend(ship_surface, 'Anomaly_sal', date_col='Date',title='Trends in Surface Ship Measurement')
# plot_anomaly_with_trend(ship_surface, 'Anomaly_pH', date_col='Date',title='Trends in Surface Ship Measurement')
# plot_anomaly_with_trend(ship_surface, 'Anomaly_aragonite', date_col='Date',title='Trends in Surface Ship Measurement')
# plot_anomaly_with_trend(ship_surface, 'Anomaly_pco2', date_col='Date',title='Trends in Surface Ship Measurement')

# plot_anomaly_with_trend(ship_surface, 'Anomaly_temp', date_col='Date',title='East Bank HOBO')
# plot_anomaly_with_trend(ship_surface, 'Anomaly_sal', date_col='Date',title='West Bank HOBO')
# plot_anomaly_with_trend(ship_surface, 'Anomaly_pco2', date_col='Date',title='East Bank SOCAT')
# plot_anomaly_with_trend(ship_surface, 'Anomaly_pco2', date_col='Date',title='West Bank SOCAT')

# --- Ship Surface: all five variables ---
plot_anomaly_with_trend(ship_surface, 'Anomaly_temp',      date_col='Date', title='Trends in Surface Ship Measurements — Temperature')
plot_anomaly_with_trend(ship_surface, 'Anomaly_sal',       date_col='Date', title='Trends in Surface Ship Measurements — Salinity')
plot_anomaly_with_trend(ship_surface, 'Anomaly_pH',        date_col='Date', title='Trends in Surface Ship Measurements — pH')
plot_anomaly_with_trend(ship_surface, 'Anomaly_aragonite', date_col='Date', title='Trends in Surface Ship Measurements — Aragonite')
plot_anomaly_with_trend(ship_surface, 'Anomaly_pco2',      date_col='Date', title='Trends in Surface Ship Measurements — pCO₂')

# --- East Bank HOBO ---
plot_anomaly_with_trend(df_east, 'Anomaly_Temp_E', date_col='Date', title='East Bank HOBO — Temperature')

# --- West Bank HOBO ---
plot_anomaly_with_trend(df_west, 'Anomaly_Temp_W', date_col='Date', title='West Bank HOBO — Temperature')

# --- East Bank SOCAT ---
plot_anomaly_with_trend(east, 'Anomaly_pco2', date_col='Date', title='East Bank SOCAT — pCO₂')
plot_anomaly_with_trend(east, 'Anomaly_temp', date_col='Date', title='East Bank SOCAT — Temperature')

# --- West Bank SOCAT ---
plot_anomaly_with_trend(west, 'Anomaly_pco2', date_col='Date', title='West Bank SOCAT — pCO₂')
plot_anomaly_with_trend(west, 'Anomaly_temp', date_col='Date', title='West Bank SOCAT — Temperature')

# %% ugh
def plot_seasonal_diagnostics(df, col, date_col='Date'):
    """
    Diagnostic plots to evaluate the performance of compute_seasonal_trend.

    Panels:
        1. Raw data vs fitted climatology overlay
        2. Climatological mean by day-of-year (the seasonal cycle shape)
        3. Anomaly series (should have no seasonal structure remaining)
        4. Residual seasonality check — mean anomaly by day-of-year
           (should be flat around zero if removal was successful)
    """

    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 11,
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'figure.dpi': 150,
    })

    clim_col    = f'SeasonalTrend_{col}'
    anomaly_col = f'Anomaly_{col}'
    trend_col   = f'Trend_{col}'

    # Verify required columns exist
    for c in [clim_col, anomaly_col]:
        if c not in df.columns:
            raise ValueError(f"Column '{c}' not found. Run compute_seasonal_trend() first.")

    df = df.copy()
    df[date_col] = pd.to_datetime(df[date_col])
    df = df.sort_values(date_col).reset_index(drop=True)  # <-- add this line
    df['_doy'] = df[date_col].dt.dayofyear
    df['_year'] = df[date_col].dt.year

    fig, axes = plt.subplots(4, 1, figsize=(11, 16))

    # ------------------------------------------------------------------ #
    # PANEL 1: Raw data vs climatology overlay
    # ------------------------------------------------------------------ #
    ax = axes[0]
    ax.scatter(df[date_col], df[col],
               s=8, color='steelblue', alpha=0.4, linewidths=0, label='Observed', zorder=2)
    ax.plot(df[date_col], df[clim_col],
            color='firebrick', linewidth=1.5, label='Climatology', zorder=3)
    if trend_col in df.columns:
        ax.plot(df[date_col], df[clim_col] + df[trend_col],
                color='darkorange', linewidth=1.2, linestyle='--',
                label='Climatology + Trend', zorder=4)
    ax.set_ylabel(col, labelpad=6)
    ax.set_title('Panel 1 — Observed vs Climatology', fontweight='bold', fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=10))
    fig.autofmt_xdate(rotation=30, ha='right')

    # ------------------------------------------------------------------ #
    # PANEL 2: Climatological seasonal cycle (day-of-year)
    # ------------------------------------------------------------------ #
    ax = axes[1]
    doy_clim = df.groupby('_doy')[col].mean()

    # Plot each individual year as a faint line for context
    for year, grp in df.groupby('_year'):
        ax.plot(grp['_doy'], grp[col],
                color='steelblue', alpha=0.15, linewidth=0.8, zorder=1)

    ax.plot(doy_clim.index, doy_clim.values,
            color='firebrick', linewidth=2, label='Climatological mean', zorder=3)
    ax.plot(df.groupby('_doy')[clim_col].mean(),
            color='darkorange', linewidth=1.5, linestyle='--',
            label='Smoothed climatology', zorder=4)
    ax.set_xlabel('Day of Year', labelpad=6)
    ax.set_ylabel(col, labelpad=6)
    ax.set_title('Panel 2 — Seasonal Cycle Shape (each faint line = one year)', fontweight='bold', fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax.set_xlim(1, 366)

    # ------------------------------------------------------------------ #
    # PANEL 3: Anomaly time series + trend
    # ------------------------------------------------------------------ #
    ax = axes[2]
    ax.scatter(df[date_col], df[anomaly_col],
               s=8, color='steelblue', alpha=0.4, linewidths=0, label='Anomaly', zorder=2)
    ax.axhline(0, color='black', linestyle='--', linewidth=0.7, alpha=0.5)
    if trend_col in df.columns:
        valid = df[trend_col].notna()
        slope, _, r_value, p_value, _ = stats.linregress(
            np.arange(valid.sum()), df.loc[valid, anomaly_col]
        )
        p_str = "p < 0.001" if p_value < 0.001 else f"p = {p_value:.3f}"
        ax.plot(df.loc[valid, date_col], df.loc[valid, trend_col],
                color='firebrick', linewidth=1.5,
                label=f'Trend  |  slope = {slope:.4f},  {p_str},  R² = {r_value**2:.3f}',
                zorder=3)
    ax.set_ylabel(f'Anomaly ({col})', labelpad=6)
    ax.set_title('Panel 3 — Anomaly Series & Linear Trend', fontweight='bold', fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=10))

    # ------------------------------------------------------------------ #
    # PANEL 4: Residual seasonality check
    #          Mean anomaly by DOY — should be ~zero if removal worked
    # ------------------------------------------------------------------ #
    ax = axes[3]
    residual_doy = df.groupby('_doy')[anomaly_col].mean()
    residual_std = df.groupby('_doy')[anomaly_col].std()

    ax.fill_between(residual_doy.index,
                    residual_doy - residual_std,
                    residual_doy + residual_std,
                    color='steelblue', alpha=0.2, label='±1 SD')
    ax.plot(residual_doy.index, residual_doy.values,
            color='steelblue', linewidth=1.5, label='Mean anomaly by DOY')
    ax.axhline(0, color='firebrick', linestyle='--', linewidth=1.0,
               label='Zero (perfect removal)')

    # Flag the worst DOY
    worst_doy = residual_doy.abs().idxmax()
    worst_val = residual_doy[worst_doy]
    ax.annotate(f'Largest residual\nDOY {worst_doy}: {worst_val:.3f}',
                xy=(worst_doy, worst_val),
                xytext=(15, 15), textcoords='offset points',
                fontsize=8.5, color='firebrick',
                arrowprops=dict(arrowstyle='->', color='firebrick', lw=0.8),
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor='firebrick', alpha=0.8))

    ax.set_xlabel('Day of Year', labelpad=6)
    ax.set_ylabel(f'Mean Anomaly ({col})', labelpad=6)
    ax.set_title('Panel 4 — Residual Seasonality (flat ≈ zero means good removal)',
                 fontweight='bold', fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax.set_xlim(1, 366)

    fig.suptitle(f'Seasonal Decomposition Diagnostics — {col}',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.show()

    # ------------------------------------------------------------------ #
    # PRINTED DIAGNOSTICS
    # ------------------------------------------------------------------ #
    print("=" * 50)
    print(f"SEASONAL DIAGNOSTICS — {col}")
    print("=" * 50)
    print(f"  Years in record:       {df['_year'].nunique()}  ({df['_year'].min()}–{df['_year'].max()})")
    print(f"  DOYs with data:        {df['_doy'].nunique()} / 366")
    print(f"  Mean anomaly:          {df[anomaly_col].mean():.4f}  (should be ≈ 0)")
    print(f"  Std of anomaly:        {df[anomaly_col].std():.4f}")
    print(f"  Max residual DOY mean: {worst_val:.4f}  at DOY {worst_doy}")
    print(f"  Skew of anomaly:       {df[anomaly_col].skew():.4f}  (should be ≈ 0)")
    print("=" * 50)
# df_east['SeasonalTrend_temp']=df_east['SeasonalTrend_Temp_E'] 
# df_west['SeasonalTrend_temp']=df_west['SeasonalTrend_Temp_W']     
# plot_seasonal_diagnostics(df_east, col='Temp_E', date_col='Date')
# plot_seasonal_diagnostics(df_west, col='Temp_W', date_col='Date')

plot_seasonal_diagnostics(ship_surface, col='temp', date_col='Date')
plot_seasonal_diagnostics(df_west, col='Temp_W', date_col='Date')

# %%PCA  
# #variables = ['temp', 'sal', 'pco2', 'pH', 'DIC', 'aragonite']
# variables = ['temp', 'pco2', 'sal', 'pH']
# ship_surface1 = ship_surface[variables].dropna()
# scaler = StandardScaler()
# X = scaler.fit_transform(ship_surface1[variables])
# pca = PCA()
# pcs = pca.fit_transform(X)
# print('Explained variance:', pca.explained_variance_ratio_)
# print(pd.DataFrame(pca.components_, columns=variables))

# # Scree plot (explained variance)
# plt.figure(figsize=(7,4))
# plt.bar(range(1, len(pca.explained_variance_ratio_)+1), pca.explained_variance_ratio_*100)
# plt.plot(range(1, len(pca.explained_variance_ratio_)+1), np.cumsum(pca.explained_variance_ratio_)*100, 
#          marker='o', color='orange', label='Cumulative')
# plt.xlabel('Principal Component')
# plt.ylabel('Explained Variance (%)')
# plt.title('PCA Scree Plot')
# plt.legend()
# plt.tight_layout()
# plt.show()
# # Biplot (PC1 vs PC2 with variable loadings)
# plt.figure(figsize=(7,6))
# xs = pcs[:,0]
# ys = pcs[:,1]
# plt.scatter(xs, ys, alpha=0.7)
# for i, var in enumerate(variables):
#     plt.arrow(0, 0, 
#               pca.components_[0,i]*max(xs), 
#               pca.components_[1,i]*max(ys),
#               color='red', head_width=0.1)
#     plt.text(pca.components_[0,i]*max(xs)*1.1, 
#              pca.components_[1,i]*max(ys)*1.1, 
#              var, color='red', ha='center', va='center')
# plt.xlabel('PC1 ({:.1f}%)'.format(pca.explained_variance_ratio_[0]*100))
# plt.ylabel('PC2 ({:.1f}%)'.format(pca.explained_variance_ratio_[1]*100))
# plt.title('PCA Biplot (Ship Surface)')
# plt.grid(True)
# plt.tight_layout()
# plt.show()
# # Optional: Barplot of PC1 loadings
# plt.figure(figsize=(7,4))
# plt.bar(variables, pca.components_[0])
# plt.xlabel('Variable')
# plt.ylabel('PC1 Loading')
# plt.title('PC1 Variable Loadings')
# plt.tight_layout()
# plt.show()

# variables = ['temp', 'sal', 'pco2']
# surface_all = pd.concat([
#     east[variables].assign(Station='East'),
#     west[variables].assign(Station='West'),
#     stetson[variables].assign(Station='Stetson')
# ], ignore_index=True)
# surface_all.dropna(inplace=True)
# scaler = StandardScaler()
# X = scaler.fit_transform(surface_all[variables])
# pca = PCA()
# pcs = pca.fit_transform(X)
# print('Explained variance:', pca.explained_variance_ratio_)
# print(pd.DataFrame(pca.components_, columns=variables))
