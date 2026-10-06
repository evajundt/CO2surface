# Data notes

A running log of problems found in the raw files and what the loaders
(`scripts/loaders.py`) do about them. Raw files are never edited, and every
fix happens at load time.

## shipto2023.xlsx
| Issue | Rows | Handling |
|---|---|---|
| Longitude written as **+93.850667** (missing minus sign) | 2022-12-03, WFG, S and B | Any positive longitude flipped to west. **Worth fixing in the master sheet.** |
| -999 fill values | various | Converted to NaN |
| Only 3 stations (EFG, WFG, STET), all near the Flower Garden Banks | all | Can't support spatial clustering alone. Use as an in-situ check on SOCAT at the banks. |

## Temp E_W_1989-2024_1sheet.xlsx
| Issue | Handling |
|---|---|
| East (`Date`) and West (`Date.1`) date columns **don't line up row by row**: from row 2922 West is one year ahead (2002-01-01 next to 2003-01-01) | Each site read with its own date column, never paired by row |
| Reef-cap loggers at roughly 20 m depth, not surface | Kept as a separate `ReefLogger` source, not mixed with surface temperature |
| Logger positions not in file | Approximate EFG/WFG reef-cap positions used (`LOGGER_SITES`); replace if you have exact ones |

## socat head.csv (older SOCAT export, cut short)
| Issue | Handling |
|---|---|
| Missing values written as **-1E+34** | Anything below -1e30 becomes NaN |
| ~142k fully blank trailing rows (`,,,,`) left by Excel | Dropped |
| Upper-case column names (`LATITUDE`, `FCO2_RECOMMENDED`, `WOCE_CO2_WATER`, `YEAR`...) | Column matching ignores case |
| `DATETIME` reformatted by Excel (`5/12/2009 18:17`) | Dates built from `YEAR/MONTH/DAY/HOUR/MINUTE/SECOND` instead |
| **Expocodes mangled to `3.16E+11`** by Excel (separate cruises merged under one ID; these are fixed moorings near 30N 88-90W, all outside our box) | Can't be recovered; loader prints a warning. Re-download and don't open in Excel. |
| 12,310 good-flag rows with fCO2 = 0 or missing | Dropped |
| fCO2 up to ~4,600 uatm at salinity near 0 (river/bay water) | Kept, but run_coverage reports the count of salinity < 25 (1,041 of 11,026 in the box). Decide before clustering whether to exclude plume/bay water. |
| Sample only reaches west to -94.97 | Expected, since the file was cut short. The full file should fill the western shelf. |

Quality filtering applied to SOCAT: WOCE flag 2 and cruise QC flag A-D
(SOCAT's recommendation for accuracy better than 5 uatm). fCO2 is converted
to pCO2 with Weiss (1974), which reproduces the ship file's own conversion
(364.8 -> 366.0 uatm).

## socat_nwgom.tsv.gz (SOCAT v2026, subset of socat_netcdf.nc)
Made with `subset_socat.py` from the SOCAT v2026 DSG NetCDF (340,448
measurements, 153 cruises, 2003-2025), cut to DOMAIN + 0.5 degrees.
Uploaded to the repo root and moved into `data/raw/` by Claude.

| Item | Value / handling |
|---|---|
| Rows in file | 95,894 (all already WOCE flag 2) |
| Dropped by QC flag A-D filter | 3,270 (cruise QC flag E) |
| Points inside DOMAIN | 75,610 from 63 cruises, 19 distinct years |
| Coverage | 71% of 0.25-degree cells, 60% of 0.1-degree cells |
| Southern edge | No data south of ~26.4N, so the download was likely limited there; the domain goes to 25.8N |
| **Seasonal imbalance** | 47% of points are from September and 20% from August; Jan, Apr, Nov have 2 cruises each. Location and season are partly confounded, so season must be removed before testing location. |
| **Long-term trend** | 2003-2025 span; ocean pCO2 rises ~1.5-2 uatm/yr, so ~40 uatm over the record. Must be removed too. |
| Salinity < 25 | 2,476 of 75,596 points |
| Missing Temp / Sal | 6 / 14 points |
