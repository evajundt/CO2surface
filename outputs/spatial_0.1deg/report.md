# Spatial grouping results (0.1 deg cells)

Input: 75,610 SOCAT points in DOMAIN -> 4,450 cruise-cell visits. Seasonal removal: gam. Neighbours: distance band 27.8 km, row-standardised. 999 permutations; LISA FDR alpha 0.05.

## Moran's I by scenario

| variable   | scenario                   |   n_cells |   visits_removed |   morans_I |      z |     p |   n_High-High |   n_Low-Low |   n_High-Low |   n_Low-High |
|:-----------|:---------------------------|----------:|-----------------:|-----------:|-------:|------:|--------------:|------------:|-------------:|-------------:|
| pCO2       | All data                   |       864 |          nan     |      0.199 | 17.362 | 0.001 |            24 |          45 |           17 |            4 |
| pCO2       | Cruise-centred             |       864 |          nan     |      0.181 | 15.435 | 0.001 |            13 |          33 |           13 |            4 |
| pCO2       | Removed S < 25             |       864 |           50.000 |      0.156 | 13.664 | 0.001 |            18 |          27 |           10 |            1 |
| pCO2       | Removed S < 30             |       858 |          253.000 |      0.198 | 16.844 | 0.001 |            16 |          16 |            2 |            5 |
| pCO2       | Removed S < 33             |       819 |         1021.000 |      0.258 | 21.226 | 0.001 |            27 |          23 |            5 |            3 |
| pCO2       | Salinity-adjusted (linear) |       864 |          nan     |      0.238 | 20.719 | 0.001 |            63 |          53 |           10 |           13 |
| pCO2       | Salinity-adjusted          |       864 |          nan     |      0.183 | 16.125 | 0.001 |            25 |          22 |            7 |            3 |
| Temp       | All data                   |       864 |          nan     |      0.437 | 37.262 | 0.001 |           190 |         186 |           22 |           22 |
| Temp       | Cruise-centred             |       864 |          nan     |      0.448 | 38.197 | 0.001 |           175 |         201 |           25 |           24 |
| Temp       | Removed S < 25             |       864 |           50.000 |      0.434 | 36.889 | 0.001 |           189 |         188 |           23 |           21 |
| Temp       | Removed S < 30             |       858 |          253.000 |      0.441 | 37.343 | 0.001 |           185 |         184 |           25 |           19 |
| Temp       | Removed S < 33             |       819 |         1021.000 |      0.431 | 36.474 | 0.001 |           187 |         161 |           20 |           20 |
| Temp       | Salinity-adjusted (linear) |       864 |          nan     |      0.365 | 30.905 | 0.001 |           167 |         171 |           28 |           21 |
| Temp       | Salinity-adjusted          |       864 |          nan     |      0.402 | 34.101 | 0.001 |           172 |         177 |           21 |           21 |

## pCO2

Seasonal+trend removal: GAM (compute_seasonal_trend_gam, decimal-year time axis), fitted on visits. Trend +2.43/yr, seasonal amplitude 53.88, R^2 0.65

Residual seasonality check (mean anomaly by month; ~0 means the cycle was removed): 1: +21.1, 2: -0.3, 3: -1.4, 4: -2.4, 5: -0.1, 6: +9.2, 7: -3.8, 8: +2.1, 9: -1.6, 10: +1.2, 11: +2.5, 12: -2.9

Salinity adjustment: LOWESS (frac 0.3) of anomaly on salinity, visit level, n 4450. For comparison, a straight line gives slope +4.82 per salinity unit, r^2 0.20.

Mean anomaly by salinity band (why a straight line doesn't fit):

```
Sal   (0, 25]  (25, 28]  (28, 30]  (30, 32]  (32, 33]  (33, 34]  (34, 35]  (35, 36]  (36, 40]
mean    -82.3     -43.3     -29.0     -10.3       1.3      -0.3       1.5       3.0       7.1
size     50.0      79.0     124.0     450.0     318.0     351.0     492.0     899.0    1687.0
```

### Justification: do Stage 1 clusters coincide with fresh water?

- Spearman rho(anomaly, salinity) across cells: +0.26 (p 2.7e-14)
- Median salinity by LISA class: High-High 35.6, High-Low 31.9, Low-High 35.3, Low-Low 31.0, ns 35.5
- Kruskal-Wallis salinity across classes: H 126.2, p 2.5e-26
- Mann-Whitney salinity, clustered vs not: p 7.2e-20
- Chi-square cluster class x (S < 25): chi2 123.3, p 1.1e-25

Cross-tab (columns: cell had at least one visit with salinity < 25):

```
any visit S<25  False  True 
lisa                        
High-High          24      0
High-Low           14      3
Low-High            4      0
Low-Low            29     16
ns                756     18
```

### Agreement of LISA labels between scenarios

| comparison | shared cells | % agree | Cohen's kappa |
|---|---|---|---|
| Removed S < 25 vs Salinity-adjusted | 864 | 91.3 | 0.31 |
| All data vs Cruise-centred | 864 | 94.0 | 0.64 |
| Salinity-adjusted vs Salinity-adjusted (linear) | 864 | 89.9 | 0.53 |
| All data vs Removed S < 25 | 864 | 94.6 | 0.66 |
| All data vs Salinity-adjusted | 864 | 88.0 | 0.25 |
| Removed S < 25 vs Removed S < 30 | 858 | 91.7 | 0.22 |
| Removed S < 25 vs Removed S < 33 | 819 | 91.3 | 0.32 |

## Temp

Seasonal+trend removal: GAM (compute_seasonal_trend_gam, decimal-year time axis), fitted on visits. Trend +0.03/yr, seasonal amplitude 5.09, R^2 0.91

Residual seasonality check (mean anomaly by month; ~0 means the cycle was removed): 1: +0.1, 2: +0.1, 3: -0.2, 4: +0.6, 5: -0.1, 6: +0.4, 7: -0.2, 8: -0.0, 9: -0.0, 10: +0.3, 11: +0.2, 12: -0.0

Salinity adjustment: LOWESS (frac 0.3) of anomaly on salinity, visit level, n 4450. For comparison, a straight line gives slope +0.11 per salinity unit, r^2 0.06.

Mean anomaly by salinity band (why a straight line doesn't fit):

```
Sal   (0, 25]  (25, 28]  (28, 30]  (30, 32]  (32, 33]  (33, 34]  (34, 35]  (35, 36]  (36, 40]
mean     -0.4      -0.7      -0.8      -0.3      -0.3      -0.3      -0.2       0.1       0.3
size     50.0      79.0     124.0     450.0     318.0     351.0     492.0     899.0    1687.0
```

### Justification: do Stage 1 clusters coincide with fresh water?

- Spearman rho(anomaly, salinity) across cells: +0.34 (p 2.4e-25)
- Median salinity by LISA class: High-High 35.9, High-Low 33.7, Low-High 36.2, Low-Low 33.6, ns 35.5
- Kruskal-Wallis salinity across classes: H 192.5, p 1.5e-40
- Mann-Whitney salinity, clustered vs not: p 0.22
- Chi-square cluster class x (S < 25): chi2 13.3, p 0.0099

Cross-tab (columns: cell had at least one visit with salinity < 25):

```
any visit S<25  False  True 
lisa                        
High-High         190      0
High-Low           20      2
Low-High           22      0
Low-Low           175     11
ns                420     24
```

### Agreement of LISA labels between scenarios

| comparison | shared cells | % agree | Cohen's kappa |
|---|---|---|---|
| Removed S < 25 vs Salinity-adjusted | 864 | 92.2 | 0.88 |
| All data vs Cruise-centred | 864 | 84.0 | 0.75 |
| Salinity-adjusted vs Salinity-adjusted (linear) | 864 | 85.2 | 0.76 |
| All data vs Removed S < 25 | 864 | 98.3 | 0.97 |
| All data vs Salinity-adjusted | 864 | 93.1 | 0.89 |
| Removed S < 25 vs Removed S < 30 | 858 | 95.0 | 0.92 |
| Removed S < 25 vs Removed S < 33 | 819 | 92.1 | 0.88 |
