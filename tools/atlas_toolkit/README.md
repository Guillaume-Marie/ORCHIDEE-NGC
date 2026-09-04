# ORCHIDEE Atlas — portable toolkit

Builds a **self-contained HTML file** (`atlas.html`): a clickable map of Europe
from ORCHIDEE outputs (variable selector, year slider), where each grid cell
opens a "pixel dashboard" (time series by age class, species filter, carried /
area-weighted / flux stand ages, disturbances, forest edge). No server, no
account: the file opens by double-click, offline, and can be shared by email or
USB stick (the data travels inside the file, quantized int16/int8 + gzip).

## Requirements
- Python 3 with `numpy` and `netCDF4` (extraction step only).
- A recent browser (Firefox/Chrome/Edge) for viewing — nothing else.

## Usage
1. Edit the **CONFIGURATION** block at the top of `extract_atlas_cube.py`:
   - `SIM_DIR`: directory of the monthly outputs (`.../SBG/Output/MO`)
   - `RUN`: file prefix (`<RUN>_YYYY0101_YYYY1231_1M_stomate_history.nc`),
     without underscores
   - `OUT`: where to write the intermediate blobs (`.` by default)
2. `./build_atlas.sh` (or run the two python scripts by hand)
3. Open `atlas.html`.

## What the simulation must provide
**Required** field: `VEGET_MAX` per PFT (monthly). Everything else is
optional — a missing field (e.g. `AGE_STAND_AREA`, `AC12_F_REAL`,
`EDGE_LENGTH_DYN`, `DA_*` on a standard trunk) simply leaves the corresponding
panels empty, without failing the extraction. Fields used:
`DIA_INV, HEIGHT, AGE_STAND_BM, RDI, IND, TOTAL_M_c, AC12_F_REAL,
AGE_STAND_AREA, FOREST_MANAGED, CLASS0_VEG, FRAC_CLASS0, DA_FIRE, DA_STORM,
DA_PEST, DA_HARVEST, EDGE_LENGTH_DYN, Areas, CONTFRAC`.

## Adapting to another PFT configuration
The PFT → species × age-class mapping is coded in `extract_atlas_cube.py`:
- `TREE_BLOCKS`: list of **0-based** `(start, end)` blocks of contiguous PFTs,
  one block = one species × its age classes in order 1→4;
- `ESSENCES`: the display names, in the same order.
The shipped values correspond to the 51-PFT / 4-age-class configuration
(11 species groups). With a single age class the "class" concept loses its
meaning: the atlas still works but per-class panels show one curve.

## Known limitations
- A regular lat/lon grid is expected (tested at 1°); the map background covers
  Europe (embedded Natural Earth 50m borders, Lambert azimuthal equal-area
  projection centred 52°N/10°E).
- The `DA_*` rates are interpreted as m²-per-timestep accumulators averaged
  monthly (AED_FEEDBACK convention).
- Page size is bounded in practice by the embedded data: ~150 years × 1000
  pixels ≈ 15 MB; beyond that, reduce the window or the domain.

Origin: developed for the euage/euage2 simulations (NextGenCarbon T4.3/T5.2),
Guillaume Marie / Science-Partners, with the assistance of Claude (Anthropic),
2026.
