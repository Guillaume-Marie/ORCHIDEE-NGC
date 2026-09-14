# ORCHIDEE Atlas — portable toolkit

Builds a **self-contained HTML file** (`atlas.html`): a clickable map of ORCHIDEE
outputs (variable selector, year slider), where each grid cell opens a
"pixel plate" (time series by age class, species filter, disturbances…).
No server, no account: the file opens by double-click, offline, and can be
shared by email or USB stick (the data travels inside the file, quantized
int16/int8 + gzip).

## Requirements
- Python 3 with `numpy` and `netCDF4` (extraction); `cartopy` optional (basemap).
- A recent browser (Firefox/Chrome/Edge) for viewing — nothing else.

## Usage
1. Edit **`atlas_config.py`** — the only file to touch.
2. `./build_atlas.sh`
3. Open `atlas.html`.

## The files
| file | role |
|---|---|
| **`atlas_config.py`** | **the whole configuration**: simulation, variables, map layers, panels |
| `extract_atlas_cube.py` | extraction engine, driven by the config — nothing hard-coded |
| `make_borders.py` | clips Natural Earth coastlines and borders to the domain (optional) |
| `assemble_atlas.py` | injects the blobs, the basemap and pako into the template |
| `atlas_template.html` | template; map and plates are built from the metadata |
| `borders_eu.json` | fallback Europe basemap, used when `make_borders.py` has not run |

## Changing the variables of interest
Three lists, in `atlas_config.py`:

- **`PFT_VARS` / `GRID_VARS`** — what is extracted from the outputs. One line
  per field: file family (`file=`, see `FILES`), netCDF name — or a **list of
  names to sum**, never truncated by a missing term —, short key, label, unit,
  scale factor (negative to flip a sign convention), **spatial** aggregation
  mode (`wmean` weighted by the PFT area — i.e. for a flux given per m² of PFT,
  a value per m² of forest; `wsum` for an extensive quantity; `rate` for a
  disturbance accumulator) and **temporal** aggregation mode (`time="mean"` by
  default, `"max"`, or `"annual"` = annual total of a daily rate, with `per="s"`
  when the rate is per second — this is what GPP, NPP, respiration and NBP need).
- **`MAP_LAYERS`** — the **spatial** side: one line = one chip above the map.
  `value` names a series (`"dia"`) or a ratio of series (`"area_c{last}/area"`,
  where `{last}` is the last age class, whatever their number). Choice of colour
  ramp, decimals, anchoring at zero.
- **`PANELS`** — the **temporal** side: one line = one chart of the pixel plate.
  `kind="classes"` (one curve per class), `"stack"` (stacked areas), `"series"`
  (named curves), `"flux_age"` (flux age Σ1/f̄). Only group in one `"series"`
  quantities with the same unit **and the same denominator**.

A layer mapping a **signed** quantity (a budget, a difference: NBP, NEE, an
anomaly…) must carry `diverging=True`: the scale becomes symmetric around zero
and the ramp goes through a central neutral colour, otherwise a sink and a
source look alike. `reverse=True` flips the ramp — useful when two neighbouring
layers follow opposite sign conventions (NBP positive = sink, NEE positive =
source), so that the same colour keeps the same physical meaning.

Adding a variable = one line in `PFT_VARS`, then one line in `MAP_LAYERS`
and/or `PANELS`. Typos are reported **before** the files are read (validation
takes a fraction of a second, not a full extraction).

`FILES` declares the output file families: `main` carries the grid and the PFT
axis (`VEGET_MAX`), the others only bring variables — this is how the CMIP
diagnostics of `stomate_ipcc_history` (`nep_c`, `rhSoil`, `nbp_c`…) are fetched
on top of `stomate_history`. A year missing from a family leaves its variables
at NaN for that year without interrupting the chain.

`YEARS` restricts the time window, `DOMAIN` crops the spatial domain.
Every optional setting can be commented out without breaking the chain.

## What the simulation must provide
**Required** field: `VEGET_MAX` per PFT (monthly), plus `lat`, `lon`, `Areas`,
`CONTFRAC`. Everything else is optional.

A variable absent from the simulation — or present but without a single value
on the domain — is **dropped from the page** (`DROP_EMPTY = True`, the default):
no greyed chip, no empty panel, no wasted room in the embedded data; the series,
the layer, the panel and the slot in the per-species blob disappear together,
and the fields concerned are listed at build time and at the top of the page.
The lines of sections 3 to 5 of `atlas_config.py` can therefore stay in place:
they come back by themselves as soon as the field exists. `DROP_EMPTY = False`
keeps the empty slots — useful to compare two simulations on an identical
layout.

In every case nothing is filled with zeros instead: a weighted mean with no
finite value is NaN, not 0; a flux age with no realised rate would equal the
residence-time ceiling, so it is masked as well.

A year whose file is unreadable (truncated archive) is reported and skipped:
the atlas is built on the remaining years.

## Adapting to another configuration
- **PFTs**: `PFT_MAPPING = "auto"` reads the `veget` axis and groups the PFTs
  named `<Species>_ageNN` into species × age classes — any number of PFTs,
  species and classes. That axis is a fixed-width character array: long names
  are **truncated** in it (`..._age01`…`_age04` become four times `..._age`),
  the classes are then numbered in order of appearance. The resulting mapping
  is printed at build time, with a warning when species have unequal numbers of
  classes. Otherwise: `PFT_MAPPING = dict(essences=[…], blocks=[(1,5), …])`.
- **Resolution and domain**: read from the file. Cell size, centring and
  framing of the map follow; the projection is Lambert azimuthal equal-area
  centred on the domain, and switches to plate carrée beyond 150° of longitude.
- **Basemap**: `make_borders.py` clips the Natural Earth coastlines and borders
  (50 m, or 110 m beyond 60° of extent) to the domain. The first time a
  resolution is used and not yet cached, **cartopy downloads it** from
  naturalearth.s3.amazonaws.com; without cartopy or network, the shipped Europe
  basemap is kept.

## Known limitations
- A regular lat/lon grid is expected (tested at 0.5°, 1° and 2°).
- The `DA_*` rates are interpreted as m²-per-timestep accumulators averaged
  monthly (AED_FEEDBACK convention).
- Extraction memory: only land cells are kept in RAM, but the per-tree-PFT
  state is held there for every year — expect about
  `n_species_variables × n_years × n_tree_PFTs × n_land_cells × 4` bytes.
- Page size: ~4 MB for 3 900 cells × 19 years × 49 series; beyond ~15 MB,
  reduce the window (`YEARS`) or the domain (`DOMAIN`).

Origin: developed for the euage/euage2 simulations (NextGenCarbon T4.3/T5.2),
Guillaume Marie / Science-Partners, with the assistance of Claude (Anthropic),
2026.
