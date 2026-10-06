# CO2surface

Spatial grouping of surface pCO2 and temperature in the NW Gulf (Galveston to the
Mexican border, extended east to include the Flower Garden Banks).

## Layout
- `scripts/loaders.py` – loads ship, reef-logger and SOCAT data into one common table
  (`Date, Latitude, Longitude, Temp, pCO2, Source, Station`)
- `scripts/nwgom_coverage_check.py` – domain definition + `coverage_report()`
- `scripts/run_coverage.py` – runs the coverage check on all sources, saves `outputs/coverage_map.png`

Raw data goes in `data/raw/` (git-ignored):
- `shipto2023.xlsx`
- `Temp_E_W_1989-2024_1sheet.xlsx`
- SOCAT export (any name; pass with `--socat`)

## Run
```
pip install pandas numpy matplotlib openpyxl
python scripts/run_coverage.py --socat data/raw/<socat_export>.tsv
```
