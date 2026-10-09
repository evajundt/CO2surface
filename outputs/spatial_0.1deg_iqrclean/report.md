# Spatial grouping results (0.1 deg cells)

Input: 75,610 SOCAT points in DOMAIN -> 4,450 cruise-cell visits. Seasonal removal: gam + anomaly IQR cleaning. Neighbours: distance band 27.8 km, row-standardised. 999 permutations; LISA FDR alpha 0.05.

## Moran's I by scenario

| variable   | scenario                   |   n_cells |   visits_removed |   morans_I |      z |     p |   n_High-High |   n_Low-Low |   n_High-Low |   n_Low-High |
|:-----------|:---------------------------|----------:|-----------------:|-----------:|-------:|------:|--------------:|------------:|-------------:|-------------:|
| pCO2       | All data                   |       863 |          nan     |      0.185 | 15.046 | 0.001 |            17 |          33 |            7 |            3 |
| pCO2       | Cruise-centred             |       863 |          nan     |      0.173 | 14.146 | 0.001 |            11 |          13 |            3 |            1 |
| pCO2       | Removed S < 25             |       863 |           50.000 |      0.162 | 13.139 | 0.001 |            15 |          14 |            8 |            3 |
| pCO2       | Removed S < 30             |       857 |          253.000 |      0.193 | 16.494 | 0.001 |            16 |          18 |            3 |            3 |
| pCO2       | Removed S < 33             |       819 |         1021.000 |      0.258 | 21.226 | 0.001 |            27 |          23 |            5 |            3 |
| pCO2       | Salinity-adjusted (linear) |       863 |          nan     |      0.244 | 20.434 | 0.001 |            44 |          40 |            4 |            5 |
| pCO2       | Salinity-adjusted          |       863 |          nan     |      0.252 | 20.993 | 0.001 |            46 |          39 |            4 |            7 |
| Temp       | All data                   |       863 |          nan     |      0.409 | 35.175 | 0.001 |           205 |         181 |           18 |           29 |
| Temp       | Cruise-centred             |       863 |          nan     |      0.452 | 39.403 | 0.001 |           165 |         191 |           17 |           24 |
| Temp       | Removed S < 25             |       863 |           50.000 |      0.404 | 34.859 | 0.001 |           198 |         178 |           18 |           28 |
| Temp       | Removed S < 30             |       857 |          253.000 |      0.401 | 32.962 | 0.001 |           191 |         181 |           24 |           27 |
| Temp       | Removed S < 33             |       818 |         1021.000 |      0.391 | 32.079 | 0.001 |           188 |         159 |           25 |           24 |
| Temp       | Salinity-adjusted (linear) |       863 |          nan     |      0.381 | 32.671 | 0.001 |           194 |         177 |           22 |           33 |
| Temp       | Salinity-adjusted          |       863 |          nan     |      0.384 | 33.111 | 0.001 |           186 |         178 |           22 |           32 |

## pCO2

Seasonal+trend removal: GAM (compute_seasonal_trend_gam, decimal-year time axis), fitted on visits. Trend +2.43/yr, seasonal amplitude 53.88, R^2 0.65; 63 visits dropped as 3x IQR anomaly outliers

Residual seasonality check (mean anomaly by month; ~0 means the cycle was removed): 1: +21.1, 2: -0.3, 3: -1.4, 4: -2.4, 5: -0.1, 6: +9.2, 7: -3.8, 8: +2.1, 9: -1.6, 10: +1.2, 11: +2.5, 12: -2.9

Salinity adjustment: LOWESS (frac 0.3) of anomaly on salinity, visit level, n 4387. For comparison, a straight line gives slope +3.28 per salinity unit, r^2 0.11.

Mean anomaly by salinity band (why a straight line doesn't fit):

```
Sal   (0, 25]  (25, 28]  (28, 30]  (30, 32]  (32, 33]  (33, 34]  (34, 35]  (35, 36]  (36, 40]
mean    -33.1     -31.1     -19.6      -9.4      -0.0      -0.3       1.5       3.0       7.1
size     30.0      66.0     108.0     440.0     314.0     351.0     492.0     899.0    1687.0
```

### Justification: do Stage 1 clusters coincide with fresh water?

- Spearman rho(anomaly, salinity) across cells: +0.17 (p 4e-07)
- Median salinity by LISA class: High-High 35.5, High-Low 32.4, Low-High 35.1, Low-Low 31.6, ns 35.5
- Kruskal-Wallis salinity across classes: H 52.5, p 1.1e-10
- Mann-Whitney salinity, clustered vs not: p 1.4e-09
- Chi-square cluster class x (S < 25): chi2 105.3, p 7.3e-22

Cross-tab (columns: cell had at least one visit with salinity < 25):

```
any visit S<25  False  True 
lisa                        
High-High          17      0
High-Low            5      2
Low-High            3      0
Low-Low            24      9
ns                792     11
```

### Agreement of LISA labels between scenarios

| comparison | shared cells | % agree | Cohen's kappa |
|---|---|---|---|
| Removed S < 25 vs Salinity-adjusted | 863 | 88.9 | 0.26 |
| All data vs Cruise-centred | 863 | 94.9 | 0.48 |
| Salinity-adjusted vs Salinity-adjusted (linear) | 863 | 97.5 | 0.87 |
| All data vs Removed S < 25 | 863 | 97.3 | 0.76 |
| All data vs Salinity-adjusted | 863 | 87.7 | 0.28 |
| Removed S < 25 vs Removed S < 30 | 857 | 96.4 | 0.60 |
| Removed S < 25 vs Removed S < 33 | 819 | 95.0 | 0.55 |

## Temp

Seasonal+trend removal: GAM (compute_seasonal_trend_gam, decimal-year time axis), fitted on visits. Trend +0.03/yr, seasonal amplitude 5.09, R^2 0.91; 118 visits dropped as 3x IQR anomaly outliers

Residual seasonality check (mean anomaly by month; ~0 means the cycle was removed): 1: +0.1, 2: +0.1, 3: -0.2, 4: +0.6, 5: -0.1, 6: +0.4, 7: -0.2, 8: -0.0, 9: -0.0, 10: +0.3, 11: +0.2, 12: -0.0

Salinity adjustment: LOWESS (frac 0.3) of anomaly on salinity, visit level, n 4332. For comparison, a straight line gives slope +0.05 per salinity unit, r^2 0.02.

Mean anomaly by salinity band (why a straight line doesn't fit):

```
Sal   (0, 25]  (25, 28]  (28, 30]  (30, 32]  (32, 33]  (33, 34]  (34, 35]  (35, 36]  (36, 40]
mean      0.5      -0.1      -0.4      -0.2      -0.1      -0.2      -0.1       0.1       0.2
size     41.0      71.0     111.0     436.0     305.0     335.0     486.0     894.0    1653.0
```

### Justification: do Stage 1 clusters coincide with fresh water?

- Spearman rho(anomaly, salinity) across cells: +0.26 (p 5.6e-15)
- Median salinity by LISA class: High-High 35.8, High-Low 35.1, Low-High 36.1, Low-Low 34.0, ns 35.4
- Kruskal-Wallis salinity across classes: H 113.7, p 1.2e-23
- Mann-Whitney salinity, clustered vs not: p 0.37
- Chi-square cluster class x (S < 25): chi2 7.1, p 0.13

Cross-tab (columns: cell had at least one visit with salinity < 25):

```
any visit S<25  False  True 
lisa                        
High-High         200      5
High-Low           18      0
Low-High           28      1
Low-Low           179      2
ns                409     21
```

### Agreement of LISA labels between scenarios

| comparison | shared cells | % agree | Cohen's kappa |
|---|---|---|---|
| Removed S < 25 vs Salinity-adjusted | 863 | 91.0 | 0.86 |
| All data vs Cruise-centred | 863 | 80.8 | 0.70 |
| Salinity-adjusted vs Salinity-adjusted (linear) | 863 | 93.2 | 0.89 |
| All data vs Removed S < 25 | 863 | 98.6 | 0.98 |
| All data vs Salinity-adjusted | 863 | 91.7 | 0.87 |
| Removed S < 25 vs Removed S < 30 | 857 | 95.7 | 0.93 |
| Removed S < 25 vs Removed S < 33 | 818 | 90.0 | 0.84 |
