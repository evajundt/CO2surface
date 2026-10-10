# CO2surface

Spatial grouping of surface pCO2 and temperature in the NW Gulf (Galveston to the
Mexican border, extended east to include the Flower Garden Banks).

## Layout
- `scripts/loaders.py` – loads ship, reef-logger and SOCAT data into one common table
  (`Date, Latitude, Longitude, Temp, pCO2, Source, Station`)
- `scripts/nwgom_coverage_check.py` – domain definition + `coverage_report()`
- `scripts/run_coverage.py` – runs the coverage check on all sources, saves `outputs/coverage_map.png`
- `scripts/subset_socat.py` – cuts a large SOCAT download (.tsv/.csv or NetCDF .nc) to the NW Gulf box (+0.5°) so it fits on GitHub (website upload limit 25 MB, so write .tsv.gz); `--info` prints the file layout
- `scripts/spatial_grouping.py` – the planned two-stage analysis: Moran's I + LISA on de-seasoned
  anomalies (all data), salinity justification tests, then freshwater removal (S < 25/30/33) and
  LOWESS salinity adjustment, plus a cruise-centred check; compares them (Cohen's kappa).
  Results in `outputs/spatial_<cell>deg/report.md`
- `scripts/seasonal.py` – Eva's `compute_seasonal_trend_gam` (pyGAM, cyclic DOY + time spline) adapted for
  the spatial analysis; changes from the original are listed at the top of the file
- `scripts/ch3_trends.py` – per-year long-term trends for the ch3 datasets (ship surface, HOBO East/West):
  seasonal cycle removed, annual means, Mann-Kendall (more cautious of plain and Hamed-Rao p) + Sen's
  slope, for the full record, 2007-present and the ship window (Nov 2013 - Aug 2023). Ship depth sets
  via `--ship-depths S M B all` (all = per-cast water-column mean). No IQR cleaning. Writes `outputs/ch3_trends/`
- `reference/fianlch3analysis.py` – Eva's ch3 analysis script. Edits are marked `CHANGED (Claude)`,
  `UNCOMMENTED (Claude)` or `NOTE (Claude)`: RMSE call updated for scikit-learn 1.6+, the
  east/west/stetson and df_east trend lines uncommented, and the trend-unit issues flagged
- `DATA_NOTES.md` – log of problems found in each raw file and how they're handled

Raw data lives in `data/raw/`:
- `shipto2023.xlsx`
- `Temp E_W_1989-2024_1sheet.xlsx`
- `socat head.csv` (older, cut-short SOCAT sample; used by default)

## Run
```
pip install -r requirements.txt
python scripts/run_coverage.py                                   # uses socat head.csv
python scripts/subset_socat.py <full_socat_file> data/raw/socat_nwgom.tsv.gz
python scripts/run_coverage.py
python scripts/spatial_grouping.py                 # 0.1-degree cells
python scripts/spatial_grouping.py --cell 0.25     # coarser sensitivity run
python scripts/spatial_grouping.py --seasonal harmonic   # earlier simple seasonal fit, for comparison
python scripts/spatial_grouping.py --iqr-clean           # also drop 3x-IQR anomaly outliers
python scripts/ch3_trends.py                              # per-year ch3 trends (all ship depth sets)
python scripts/ch3_trends.py --ship-depths S all          # surface vs whole water column only
```
