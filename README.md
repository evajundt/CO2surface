# CO2surface

Spatial grouping of surface pCO2 and temperature in the NW Gulf (Galveston to the
Mexican border, extended east to include the Flower Garden Banks).

## Layout
- `scripts/loaders.py` – loads ship, reef-logger and SOCAT data into one common table
  (`Date, Latitude, Longitude, Temp, pCO2, Source, Station`)
- `scripts/nwgom_coverage_check.py` – domain definition + `coverage_report()`
- `scripts/run_coverage.py` – runs the coverage check on all sources, saves `outputs/coverage_map.png`
- `scripts/subset_socat.py` – cuts a large SOCAT download to the NW Gulf box (+0.5°) so it fits on GitHub
- `DATA_NOTES.md` – log of problems found in each raw file and how they're handled

Raw data lives in `data/raw/`:
- `shipto2023.xlsx`
- `Temp E_W_1989-2024_1sheet.xlsx`
- `socat head.csv` (older, cut-short SOCAT sample; used by default)

## Run
```
pip install pandas numpy scipy matplotlib openpyxl
python scripts/run_coverage.py                                   # uses socat head.csv
python scripts/subset_socat.py <full_socat_file> data/raw/socat_nwgom.tsv
python scripts/run_coverage.py --socat data/raw/socat_nwgom.tsv  # once the full file is in
```
