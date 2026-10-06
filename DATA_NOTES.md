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
