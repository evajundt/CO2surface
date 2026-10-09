# Spatial grouping results (0.25 deg cells)

Input: 75,610 SOCAT points in DOMAIN -> 1,669 cruise-cell visits. Seasonal removal: harmonic. Neighbours: distance band 69.4 km, row-standardised. 999 permutations; LISA FDR alpha 0.05.

## Moran's I by scenario

| variable   | scenario                   |   n_cells |   visits_removed |   morans_I |      z |     p |   n_High-High |   n_Low-Low |   n_High-Low |   n_Low-High |
|:-----------|:---------------------------|----------:|-----------------:|-----------:|-------:|------:|--------------:|------------:|-------------:|-------------:|
| pCO2       | All data                   |       171 |          nan     |      0.108 |  4.496 | 0.001 |             0 |           4 |            6 |            1 |
| pCO2       | Cruise-centred             |       171 |          nan     |      0.121 |  5.077 | 0.001 |             1 |           4 |            5 |            1 |
| pCO2       | Removed S < 25             |       171 |           19.000 |      0.010 |  0.586 | 0.264 |             0 |           0 |            0 |            0 |
| pCO2       | Removed S < 30             |       171 |           98.000 |      0.127 |  5.315 | 0.001 |             0 |           0 |            0 |            0 |
| pCO2       | Removed S < 33             |       167 |          382.000 |      0.144 |  5.739 | 0.001 |             0 |           0 |            0 |            0 |
| pCO2       | Salinity-adjusted (linear) |       171 |          nan     |      0.121 |  4.884 | 0.001 |             4 |          18 |            6 |            4 |
| pCO2       | Salinity-adjusted          |       171 |          nan     |      0.071 |  2.985 | 0.006 |             0 |           0 |            0 |            0 |
| Temp       | All data                   |       171 |          nan     |      0.506 | 19.816 | 0.001 |            50 |          39 |            6 |            3 |
| Temp       | Cruise-centred             |       171 |          nan     |      0.438 | 17.946 | 0.001 |            49 |          52 |            9 |            5 |
| Temp       | Removed S < 25             |       171 |           19.000 |      0.492 | 19.360 | 0.001 |            50 |          44 |            5 |            3 |
| Temp       | Removed S < 30             |       171 |           98.000 |      0.451 | 18.009 | 0.001 |            49 |          40 |            8 |            3 |
| Temp       | Removed S < 33             |       167 |          382.000 |      0.475 | 18.310 | 0.001 |            44 |          26 |            2 |            3 |
| Temp       | Salinity-adjusted (linear) |       171 |          nan     |      0.325 | 12.905 | 0.001 |            36 |          22 |            5 |            2 |
| Temp       | Salinity-adjusted          |       171 |          nan     |      0.442 | 17.425 | 0.001 |            43 |          29 |            4 |            5 |

## pCO2

Seasonal+trend removal: 2-harmonic + linear trend, fitted on visits. Trend +2.19/yr, seasonal amplitude 49.01, R^2 0.60

Salinity adjustment: LOWESS (frac 0.3) of anomaly on salinity, visit level, n 1669. For comparison, a straight line gives slope +5.37 per salinity unit, r^2 0.22.

Mean anomaly by salinity band (why a straight line doesn't fit):

```
Sal   (0, 25]  (25, 28]  (28, 30]  (30, 32]  (32, 33]  (33, 34]  (34, 35]  (35, 36]  (36, 40]
mean    -84.7     -52.7     -39.3      -7.9      -1.0      -2.6       3.4       4.1       7.9
size     19.0      32.0      47.0     174.0     110.0     139.0     182.0     333.0     633.0
```

### Justification: do Stage 1 clusters coincide with fresh water?

- Spearman rho(anomaly, salinity) across cells: +0.35 (p 2.6e-06)
- Median salinity by LISA class: High-Low 32.5, Low-High 34.6, Low-Low 29.6, ns 35.4
- Kruskal-Wallis salinity across classes: H 17.5, p 0.00016
- Mann-Whitney salinity, clustered vs not: p 3.8e-05
- Chi-square cluster class x (S < 25): chi2 27.5, p 4.5e-06

Cross-tab (columns: cell had at least one visit with salinity < 25):

```
any visit S<25  False  True 
lisa                        
High-Low            5      1
Low-High            1      0
Low-Low             1      3
ns                151      9
```

### Agreement of LISA labels between scenarios

| comparison | shared cells | % agree | Cohen's kappa |
|---|---|---|---|
| Removed S < 25 vs Salinity-adjusted | 171 | 100.0 | nan |
| All data vs Cruise-centred | 171 | 98.2 | 0.86 |
| Salinity-adjusted vs Salinity-adjusted (linear) | 171 | 81.3 | 0.00 |
| All data vs Removed S < 25 | 171 | 93.6 | 0.00 |
| All data vs Salinity-adjusted | 171 | 93.6 | 0.00 |
| Removed S < 25 vs Removed S < 30 | 171 | 100.0 | nan |
| Removed S < 25 vs Removed S < 33 | 167 | 100.0 | nan |

## Temp

Seasonal+trend removal: 2-harmonic + linear trend, fitted on visits. Trend +0.02/yr, seasonal amplitude 5.14, R^2 0.88

Salinity adjustment: LOWESS (frac 0.3) of anomaly on salinity, visit level, n 1669. For comparison, a straight line gives slope +0.12 per salinity unit, r^2 0.06.

Mean anomaly by salinity band (why a straight line doesn't fit):

```
Sal   (0, 25]  (25, 28]  (28, 30]  (30, 32]  (32, 33]  (33, 34]  (34, 35]  (35, 36]  (36, 40]
mean     -0.4      -0.8      -0.9      -0.3      -0.5      -0.4      -0.2       0.1       0.4
size     19.0      32.0      47.0     174.0     110.0     139.0     182.0     333.0     633.0
```

### Justification: do Stage 1 clusters coincide with fresh water?

- Spearman rho(anomaly, salinity) across cells: +0.47 (p 1e-10)
- Median salinity by LISA class: High-High 35.6, High-Low 32.4, Low-High 36.2, Low-Low 32.8, ns 35.5
- Kruskal-Wallis salinity across classes: H 67.3, p 8.4e-14
- Mann-Whitney salinity, clustered vs not: p 0.029
- Chi-square cluster class x (S < 25): chi2 40.9, p 2.9e-08

Cross-tab (columns: cell had at least one visit with salinity < 25):

```
any visit S<25  False  True 
lisa                        
High-High          50      0
High-Low            5      1
Low-High            3      0
Low-Low            27     12
ns                 73      0
```

### Agreement of LISA labels between scenarios

| comparison | shared cells | % agree | Cohen's kappa |
|---|---|---|---|
| Removed S < 25 vs Salinity-adjusted | 171 | 86.0 | 0.79 |
| All data vs Cruise-centred | 171 | 79.5 | 0.71 |
| Salinity-adjusted vs Salinity-adjusted (linear) | 171 | 85.4 | 0.76 |
| All data vs Removed S < 25 | 171 | 97.1 | 0.96 |
| All data vs Salinity-adjusted | 171 | 88.3 | 0.82 |
| Removed S < 25 vs Removed S < 30 | 171 | 95.3 | 0.93 |
| Removed S < 25 vs Removed S < 33 | 167 | 85.0 | 0.77 |
