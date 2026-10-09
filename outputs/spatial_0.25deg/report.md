# Spatial grouping results (0.25 deg cells)

Input: 75,610 SOCAT points in DOMAIN -> 1,669 cruise-cell visits. Seasonal removal: gam. Neighbours: distance band 69.4 km, row-standardised. 999 permutations; LISA FDR alpha 0.05.

## Moran's I by scenario

| variable   | scenario                   |   n_cells |   visits_removed |   morans_I |      z |     p |   n_High-High |   n_Low-Low |   n_High-Low |   n_Low-High |
|:-----------|:---------------------------|----------:|-----------------:|-----------:|-------:|------:|--------------:|------------:|-------------:|-------------:|
| pCO2       | All data                   |       171 |          nan     |      0.098 |  4.103 | 0.001 |             1 |           7 |            6 |            1 |
| pCO2       | Cruise-centred             |       171 |          nan     |      0.114 |  4.811 | 0.001 |             1 |           2 |            4 |            1 |
| pCO2       | Removed S < 25             |       171 |           19.000 |     -0.001 |  0.155 | 0.407 |             0 |           0 |            0 |            0 |
| pCO2       | Removed S < 30             |       171 |           98.000 |      0.104 |  4.356 | 0.001 |             0 |           0 |            0 |            0 |
| pCO2       | Removed S < 33             |       167 |          382.000 |      0.122 |  4.903 | 0.001 |             0 |           0 |            0 |            0 |
| pCO2       | Salinity-adjusted (linear) |       171 |          nan     |      0.114 |  4.704 | 0.001 |             3 |          13 |            4 |            1 |
| pCO2       | Salinity-adjusted          |       171 |          nan     |      0.048 |  2.114 | 0.027 |             0 |           0 |            0 |            0 |
| Temp       | All data                   |       171 |          nan     |      0.502 | 19.476 | 0.001 |            53 |          44 |            7 |            0 |
| Temp       | Cruise-centred             |       171 |          nan     |      0.441 | 18.083 | 0.001 |            50 |          53 |            9 |            3 |
| Temp       | Removed S < 25             |       171 |           19.000 |      0.488 | 19.192 | 0.001 |            53 |          51 |            5 |            0 |
| Temp       | Removed S < 30             |       171 |           98.000 |      0.437 | 17.061 | 0.001 |            49 |          44 |            6 |            2 |
| Temp       | Removed S < 33             |       167 |          382.000 |      0.472 | 18.382 | 0.001 |            46 |          25 |            2 |            1 |
| Temp       | Salinity-adjusted (linear) |       171 |          nan     |      0.327 | 12.755 | 0.001 |            38 |          28 |            7 |            3 |
| Temp       | Salinity-adjusted          |       171 |          nan     |      0.442 | 17.198 | 0.001 |            44 |          38 |            5 |            4 |

## pCO2

Seasonal+trend removal: GAM (compute_seasonal_trend_gam, decimal-year time axis), fitted on visits. Trend +2.30/yr, seasonal amplitude 52.30, R^2 0.66

Residual seasonality check (mean anomaly by month; ~0 means the cycle was removed): 1: +24.4, 2: -0.8, 3: -2.5, 4: -1.5, 5: +1.5, 6: +7.5, 7: -4.5, 8: +1.8, 9: -1.6, 10: +1.9, 11: +3.5, 12: -3.2

Salinity adjustment: LOWESS (frac 0.3) of anomaly on salinity, visit level, n 1669. For comparison, a straight line gives slope +4.91 per salinity unit, r^2 0.22.

Mean anomaly by salinity band (why a straight line doesn't fit):

```
Sal   (0, 25]  (25, 28]  (28, 30]  (30, 32]  (32, 33]  (33, 34]  (34, 35]  (35, 36]  (36, 40]
mean    -76.8     -48.4     -35.1      -7.3       0.5      -3.1       2.7       2.6       7.8
size     19.0      32.0      47.0     174.0     110.0     139.0     182.0     333.0     633.0
```

### Justification: do Stage 1 clusters coincide with fresh water?

- Spearman rho(anomaly, salinity) across cells: +0.35 (p 2.9e-06)
- Median salinity by LISA class: High-High 35.3, High-Low 32.5, Low-High 34.6, Low-Low 30.3, ns 35.4
- Kruskal-Wallis salinity across classes: H 25.2, p 3.5e-06
- Mann-Whitney salinity, clustered vs not: p 1.7e-06
- Chi-square cluster class x (S < 25): chi2 64.8, p 2.8e-13

Cross-tab (columns: cell had at least one visit with salinity < 25):

```
any visit S<25  False  True 
lisa                        
High-High           1      0
High-Low            5      1
Low-High            1      0
Low-Low             1      6
ns                150      6
```

### Agreement of LISA labels between scenarios

| comparison | shared cells | % agree | Cohen's kappa |
|---|---|---|---|
| Removed S < 25 vs Salinity-adjusted | 171 | 100.0 | nan |
| All data vs Cruise-centred | 171 | 95.9 | 0.68 |
| Salinity-adjusted vs Salinity-adjusted (linear) | 171 | 87.7 | 0.00 |
| All data vs Removed S < 25 | 171 | 91.2 | 0.00 |
| All data vs Salinity-adjusted | 171 | 91.2 | 0.00 |
| Removed S < 25 vs Removed S < 30 | 171 | 100.0 | nan |
| Removed S < 25 vs Removed S < 33 | 167 | 100.0 | nan |

## Temp

Seasonal+trend removal: GAM (compute_seasonal_trend_gam, decimal-year time axis), fitted on visits. Trend +0.03/yr, seasonal amplitude 4.98, R^2 0.90

Residual seasonality check (mean anomaly by month; ~0 means the cycle was removed): 1: +0.0, 2: +0.1, 3: -0.1, 4: +0.4, 5: -0.1, 6: +0.5, 7: -0.2, 8: -0.0, 9: -0.0, 10: +0.3, 11: +0.2, 12: -0.0

Salinity adjustment: LOWESS (frac 0.3) of anomaly on salinity, visit level, n 1669. For comparison, a straight line gives slope +0.12 per salinity unit, r^2 0.07.

Mean anomaly by salinity band (why a straight line doesn't fit):

```
Sal   (0, 25]  (25, 28]  (28, 30]  (30, 32]  (32, 33]  (33, 34]  (34, 35]  (35, 36]  (36, 40]
mean     -0.6      -0.7      -0.9      -0.3      -0.4      -0.4      -0.2       0.1       0.4
size     19.0      32.0      47.0     174.0     110.0     139.0     182.0     333.0     633.0
```

### Justification: do Stage 1 clusters coincide with fresh water?

- Spearman rho(anomaly, salinity) across cells: +0.52 (p 5.5e-13)
- Median salinity by LISA class: High-High 35.8, High-Low 31.9, Low-Low 33.2, ns 35.4
- Kruskal-Wallis salinity across classes: H 68.8, p 7.5e-15
- Mann-Whitney salinity, clustered vs not: p 0.18
- Chi-square cluster class x (S < 25): chi2 26.6, p 7e-06

Cross-tab (columns: cell had at least one visit with salinity < 25):

```
any visit S<25  False  True 
lisa                        
High-High          53      0
High-Low            5      2
Low-Low            34     10
ns                 66      1
```

### Agreement of LISA labels between scenarios

| comparison | shared cells | % agree | Cohen's kappa |
|---|---|---|---|
| Removed S < 25 vs Salinity-adjusted | 171 | 86.0 | 0.79 |
| All data vs Cruise-centred | 171 | 84.8 | 0.78 |
| Salinity-adjusted vs Salinity-adjusted (linear) | 171 | 82.5 | 0.73 |
| All data vs Removed S < 25 | 171 | 95.9 | 0.94 |
| All data vs Salinity-adjusted | 171 | 89.5 | 0.84 |
| Removed S < 25 vs Removed S < 30 | 171 | 92.4 | 0.89 |
| Removed S < 25 vs Removed S < 33 | 167 | 79.6 | 0.69 |
