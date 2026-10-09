# Spatial grouping results (0.1 deg cells)

Input: 75,610 SOCAT points in DOMAIN -> 4,450 cruise-cell visits. Seasonal removal: harmonic. Neighbours: distance band 27.8 km, row-standardised. 999 permutations; LISA FDR alpha 0.05.

## Moran's I by scenario

| variable   | scenario                   |   n_cells |   visits_removed |   morans_I |      z |     p |   n_High-High |   n_Low-Low |   n_High-Low |   n_Low-High |
|:-----------|:---------------------------|----------:|-----------------:|-----------:|-------:|------:|--------------:|------------:|-------------:|-------------:|
| pCO2       | All data                   |       864 |          nan     |      0.200 | 17.325 | 0.001 |            35 |          47 |           18 |            4 |
| pCO2       | Cruise-centred             |       864 |          nan     |      0.186 | 15.835 | 0.001 |            17 |          34 |           15 |            4 |
| pCO2       | Removed S < 25             |       864 |           50.000 |      0.151 | 13.082 | 0.001 |            16 |          17 |            7 |            1 |
| pCO2       | Removed S < 30             |       858 |          253.000 |      0.210 | 17.559 | 0.001 |            23 |           7 |            6 |            1 |
| pCO2       | Removed S < 33             |       819 |         1021.000 |      0.265 | 21.678 | 0.001 |            50 |          85 |            6 |            7 |
| pCO2       | Salinity-adjusted (linear) |       864 |          nan     |      0.253 | 21.747 | 0.001 |            65 |          95 |           10 |           12 |
| pCO2       | Salinity-adjusted          |       864 |          nan     |      0.201 | 17.355 | 0.001 |            42 |          34 |            9 |            6 |
| Temp       | All data                   |       864 |          nan     |      0.436 | 38.247 | 0.001 |           197 |         191 |           21 |           33 |
| Temp       | Cruise-centred             |       864 |          nan     |      0.448 | 38.211 | 0.001 |           180 |         202 |           23 |           24 |
| Temp       | Removed S < 25             |       864 |           50.000 |      0.431 | 37.370 | 0.001 |           199 |         197 |           22 |           32 |
| Temp       | Removed S < 30             |       858 |          253.000 |      0.423 | 35.247 | 0.001 |           173 |         195 |           24 |           29 |
| Temp       | Removed S < 33             |       819 |         1021.000 |      0.422 | 35.916 | 0.001 |           162 |         154 |           25 |           28 |
| Temp       | Salinity-adjusted (linear) |       864 |          nan     |      0.354 | 30.938 | 0.001 |           132 |         163 |           23 |           24 |
| Temp       | Salinity-adjusted          |       864 |          nan     |      0.399 | 34.885 | 0.001 |           176 |         180 |           21 |           30 |

## pCO2

Seasonal+trend removal: 2-harmonic + linear trend, fitted on visits. Trend +2.27/yr, seasonal amplitude 48.60, R^2 0.58

Salinity adjustment: LOWESS (frac 0.3) of anomaly on salinity, visit level, n 4450. For comparison, a straight line gives slope +5.33 per salinity unit, r^2 0.21.

Mean anomaly by salinity band (why a straight line doesn't fit):

```
Sal   (0, 25]  (25, 28]  (28, 30]  (30, 32]  (32, 33]  (33, 34]  (34, 35]  (35, 36]  (36, 40]
mean    -91.8     -49.0     -33.2     -10.8       0.1      -0.9       3.2       4.6       7.1
size     50.0      79.0     124.0     450.0     318.0     351.0     492.0     899.0    1687.0
```

### Justification: do Stage 1 clusters coincide with fresh water?

- Spearman rho(anomaly, salinity) across cells: +0.23 (p 4.6e-12)
- Median salinity by LISA class: High-High 35.6, High-Low 32.6, Low-High 36.3, Low-Low 31.0, ns 35.5
- Kruskal-Wallis salinity across classes: H 138.5, p 6e-29
- Mann-Whitney salinity, clustered vs not: p 7.4e-17
- Chi-square cluster class x (S < 25): chi2 141.3, p 1.5e-29

Cross-tab (columns: cell had at least one visit with salinity < 25):

```
any visit S<25  False  True 
lisa                        
High-High          35      0
High-Low           14      4
Low-High            4      0
Low-Low            30     17
ns                744     16
```

### Agreement of LISA labels between scenarios

| comparison | shared cells | % agree | Cohen's kappa |
|---|---|---|---|
| Removed S < 25 vs Salinity-adjusted | 864 | 89.4 | 0.27 |
| All data vs Cruise-centred | 864 | 92.9 | 0.63 |
| Salinity-adjusted vs Salinity-adjusted (linear) | 864 | 88.0 | 0.58 |
| All data vs Removed S < 25 | 864 | 92.2 | 0.52 |
| All data vs Salinity-adjusted | 864 | 85.2 | 0.29 |
| Removed S < 25 vs Removed S < 30 | 858 | 92.2 | 0.10 |
| Removed S < 25 vs Removed S < 33 | 819 | 80.6 | 0.10 |

## Temp

Seasonal+trend removal: 2-harmonic + linear trend, fitted on visits. Trend +0.02/yr, seasonal amplitude 5.11, R^2 0.88

Salinity adjustment: LOWESS (frac 0.3) of anomaly on salinity, visit level, n 4450. For comparison, a straight line gives slope +0.12 per salinity unit, r^2 0.05.

Mean anomaly by salinity band (why a straight line doesn't fit):

```
Sal   (0, 25]  (25, 28]  (28, 30]  (30, 32]  (32, 33]  (33, 34]  (34, 35]  (35, 36]  (36, 40]
mean     -0.3      -0.8      -0.9      -0.3      -0.4      -0.4      -0.2       0.1       0.4
size     50.0      79.0     124.0     450.0     318.0     351.0     492.0     899.0    1687.0
```

### Justification: do Stage 1 clusters coincide with fresh water?

- Spearman rho(anomaly, salinity) across cells: +0.29 (p 5.8e-18)
- Median salinity by LISA class: High-High 35.7, High-Low 33.4, Low-High 36.0, Low-Low 33.6, ns 35.4
- Kruskal-Wallis salinity across classes: H 151.9, p 7.9e-32
- Mann-Whitney salinity, clustered vs not: p 0.041
- Chi-square cluster class x (S < 25): chi2 18.3, p 0.0011

Cross-tab (columns: cell had at least one visit with salinity < 25):

```
any visit S<25  False  True 
lisa                        
High-High         197      0
High-Low           20      1
Low-High           33      0
Low-Low           175     16
ns                402     20
```

### Agreement of LISA labels between scenarios

| comparison | shared cells | % agree | Cohen's kappa |
|---|---|---|---|
| Removed S < 25 vs Salinity-adjusted | 864 | 90.7 | 0.86 |
| All data vs Cruise-centred | 864 | 68.9 | 0.52 |
| Salinity-adjusted vs Salinity-adjusted (linear) | 864 | 85.2 | 0.76 |
| All data vs Removed S < 25 | 864 | 98.0 | 0.97 |
| All data vs Salinity-adjusted | 864 | 91.8 | 0.87 |
| Removed S < 25 vs Removed S < 30 | 858 | 91.3 | 0.87 |
| Removed S < 25 vs Removed S < 33 | 819 | 86.4 | 0.79 |
