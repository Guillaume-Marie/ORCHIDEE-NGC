#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 ATLAS CONFIGURATION  —  the only file to edit
================================================================================
Everything else adapts: number of PFTs, number of age classes, number of
species, grid resolution and grid domain are all read from the files.

Five sections:
  1. SIMULATION  — where to read, which years, which domain
  2. STRUCTURE   — how the PFTs group into species × age classes
  3. VARIABLES   — what is extracted from the ORCHIDEE outputs
  4. MAP         — the SPATIAL layers (one chip = one map)
  5. PLATE       — the TEMPORAL panels (one entry = one chart)

Adding a variable of interest = 1 line in section 3, then 1 line in
section 4 (to see it on the map) and/or 1 line in section 5 (to see it as a
time series). A variable missing from the simulation is reported at build
time and leaves its panels empty — it breaks nothing.
"""

# ══════════════════════════════════════════════════════════════════════════════
#  1. SIMULATION
# ══════════════════════════════════════════════════════════════════════════════

SIM_DIR = "/home/mguill/Téléchargements/atlas_toolkit/ISIMIP0p5spinup/SBG/Output/YE"
RUN     = "ISIMIP0p5spinup"   # prefix: <RUN>_YYYY0101_YYYY1231_<suffix>.nc

# Output file families. "main" carries the grid and the PFT axis (VEGET_MAX);
# the others only contribute variables, referenced by `file=` in section 3.
FILES = dict(main="1Y_stomate_history",
             ipcc="1Y_stomate_ipcc_history")
OUT     = "."                  # where to write the intermediate blobs

YEARS   = None        # None = all years found; otherwise (1980, 2020)
DOMAIN  = None        # None = the whole domain of the file;
                      # otherwise dict(lat=(35, 72), lon=(-12, 33)) to crop

TITLE   = None        # None = "ORCHIDEE Atlas — <RUN>"; otherwise a string

# A variable missing from the simulation — or present but without a single value
# over the domain — is dropped from the page: no greyed-out chip, no empty panel,
# no space wasted in the embedded data. The lines of sections 3 to 5 can
# therefore stay in place: they come back by themselves as soon as the field
# exists. False to keep the empty placeholders (useful to compare two
# simulations on an identical layout).
DROP_EMPTY = True


# ══════════════════════════════════════════════════════════════════════════════
#  2. PFT STRUCTURE  →  species × age classes
# ══════════════════════════════════════════════════════════════════════════════
# "auto": the names of the `veget` axis are read from the file. Any PFT named
#          "<Species>_ageNN" forms, with its siblings, a species with N classes.
#          Works for 39 PFTs/4 classes, 51/4, 15/1, 79/6… without any change.
# Otherwise, explicit mapping:
#     PFT_MAPPING = dict(essences=["Temp. pine", "Beech"], blocks=[(1, 5), (5, 9)])
#     (0-based [start, end) blocks of contiguous PFTs, classes in order 1→N)

PFT_MAPPING = "auto"

# Display names for ORCHIDEE PFTs (optional: anything missing is shown as is,
# automatically shortened). Key = PFT name without the _ageNN suffix.
ESSENCE_LABELS = {
    "BroadLeavedEvergreenTropical":   "Trop. evergreen",
    "BroadLeavedRaingreenTropical":   "Trop. deciduous",
    "NeedleleafEvergreenTemperate":   "Temp. conifer",
    "BroadLeavedEvergreenTemperate":  "Temp. evergreen broadl.",
    "BroadLeavedSummergreenTemperate":"Temp. broadleaf",
    "NeedleleafEvergreenBoreal":      "Boreal conifer",
    "BroadLeavedSummergreenBoreal":   "Boreal broadleaf",
    "LarixSpBoreal":                  "Larch",
    "NeedleleafEvergreenTemperateSpruce": "Temp. spruce",
    "NeedleleafEvergreenBorealSpruce":    "Boreal spruce",
    "BroadLeavedEvergreenEucalyptus":     "Eucalyptus",
}

# A pixel only displays the species exceeding this share of forest cover,
# and at most this many (size of the per-species blob).
ESSENCE_MIN_SHARE   = 0.03
ESSENCE_MAX_PER_PIX = 4


# ══════════════════════════════════════════════════════════════════════════════
#  3. VARIABLES OF INTEREST
# ══════════════════════════════════════════════════════════════════════════════
# PFT_VARS: fields of dimension (time, veget, lat, lon).
#   file    file family to read from (default "main", see FILES in section 1)
#   nc      name in the file — or a LIST of names to be summed
#           (e.g. ["MAINT_RESP", "GROWTH_RESP"]; the sum is only computed
#           if all its terms exist, never truncated by a missing term)
#   key     short name, used in sections 4 and 5
#   label   / unit  : what is displayed (unit after scale is applied)
#   scale   factor applied for display (m→cm: 100; g→kg: 1e-3)
#   agg     "wmean" mean weighted by the PFT areas      (intensive quantity)
#           "wsum"  sum of the values × PFT areas         (extensive quantity)
#   time    aggregation of the 12 months: "mean" (default), "sum", "max", "first",
#           "annual" = annual total of a daily rate (Σ months × days in month;
#           add per="s" if the rate is per second)
#   classes True (default) = one series per age class + the pixel aggregate
#   essence True = also stored per species (species filter of the plate)
#   clip    physical ceiling before quantization (prevents an empty slot
#           with an aberrant value from crushing the resolution of all others)
#   step    quantization step of the per-species blob (int8, 250 levels)

PFT_VARS = [
  dict(nc="DIA_INV",        key="dia",  label="Quadratic diameter", unit="cm",
       scale=100,  agg="wmean", essence=True, step=0.8/250),
  dict(nc="HEIGHT",         key="hgt",  label="Height",              unit="m",
       agg="wmean", essence=True, step=45/250),
  #dict(nc="AGE_STAND_BM",   key="age",  label="Carried age (AGE_STAND_BM)", unit="yr",
  #     agg="wmean", essence=True, step=300/250),
  #dict(nc="AGE_STAND_AREA", key="aga",  label="Area-weighted age (AGE_STAND_AREA)", unit="yr",
   #    agg="wmean", essence=True, step=300/250),
  dict(nc="RDI",            key="rdi",  label="RDI",                  unit="–",
       agg="wmean", essence=True, clip=3.0, step=2.5/250),
  dict(nc="IND",            key="ind",  label="Stem density",     unit="ind/ha",
       scale=1e4, agg="wmean", clip=1.0),
  #dict(nc="AC12_F_REAL",    key="ac12", label="Realized transfer rate", unit="1/yr",
  #     agg="wmean", essence=True, step=0.5/250),
  dict(nc="TOTAL_M_c",      key="biomass", label="Living biomass (forest)", unit="kgC/m²",
       scale=1e-3, agg="wsum", classes=False),

  # Carbon fluxes. They are given per m² of PFT: agg="wmean" thus turns them into
  # a productivity PER M² OF FOREST — comparable from one grid cell and one age
  # class to the next. (agg="wsum" would give the share of the grid-cell flux
  # attributable to forest, whose map would mostly redraw the forest cover rate.)
  dict(nc="GPP",       key="gpp", label="Gross productivity (GPP)", unit="gC/m² forest/yr",
       agg="wmean", time="annual"),
  dict(nc="NPP",       key="npp", label="Net productivity (NPP)", unit="gC/m² forest/yr",
       agg="wmean", time="annual"),
  dict(nc=["MAINT_RESP", "GROWTH_RESP"],
       key="ra",  label="Autotrophic respiration", unit="gC/m² forest/yr",
       agg="wmean", time="annual"),
  dict(nc="HET_RESP",  key="rh",  label="Heterotrophic respiration", unit="gC/m² forest/yr",
       agg="wmean", time="annual"),
  dict(nc="LAI_MAX",   key="lai", label="Maximum LAI", unit="m²/m²",
       agg="wmean", time="max", essence=True, step=12/250),
]

# GRID_VARS: fields of dimension (time, lat, lon) — no age classes.
#   agg  "mean" annual mean of the months
#        "rate" m²/time-step accumulator → annual area / land area
#        "sum"  sum of the months
GRID_VARS = [
  dict(nc="DA_FIRE",         key="da_fire",  label="Fire",      unit="%/yr", scale=100, agg="rate", clip=0.5),
  dict(nc="DA_STORM",        key="da_storm", label="Storm",  unit="%/yr", scale=100, agg="rate", clip=0.5),
  dict(nc="DA_PEST",         key="da_pest",  label="Bark beetle",  unit="%/yr", scale=100, agg="rate", clip=0.5),
  dict(nc="DA_HARVEST",      key="da_harv",  label="Harvest",  unit="%/yr", scale=100, agg="rate", clip=0.5),
  #dict(nc="EDGE_LENGTH_DYN", key="edge",     label="Edge length", unit="km",
  #     scale=1e-3, agg="mean"),

  # NEE and soil respiration come from the ipcc file (CMIP names), at the level
  # of the WHOLE GRID CELL like NBP.
  #  · nep_c is counted positive towards the land (sink); NEE follows the opposite
  #    convention, hence the negative scale: positive = emission to the atmosphere.
  #  · rhSoil is the SOIL share of heterotrophic respiration (rhSoil + rhLitter
  #    = rh); root respiration is not diagnosed separately, so this is not
  #    the total soil respiration measured in a chamber.
  dict(nc="nep_c",  key="nee",    file="ipcc", label="NEE (whole grid cell)",
       unit="gC/m²/yr", scale=-1000, agg="annual", per="s"),
  dict(nc="rhSoil", key="rhsoil", file="ipcc", label="Soil respiration (heterotrophic)",
       unit="gC/m²/yr", scale=1000, agg="annual", per="s"),

  # NBP: net balance of the WHOLE GRID CELL (forest, crops, grasslands, wood
  # products) — not attributable to forest alone. Given in kgC/s/m², hence per="s"
  # and scale=1000 to output gC/m²/yr. Positive = sink.
  dict(nc="NBP_pool_c", key="nbp", label="NBP (whole grid cell)", unit="gC/m²/yr",
       scale=1000, agg="annual", per="s"),
]

# The area of the tree PFTs (VEGET_MAX) is always extracted, under the key "area":
#   area_c1..area_cN  area of each class,  area  forest area of the pixel.
AREA_LABEL = dict(label="Forest area", unit="% of grid cell", scale=100)

# Flux age Σ1/f̄: age estimator derived from the realized transfer rates,
# insensitive to cohort mixing. None to disable it.
#   rate  key of the rate variable (PFT_VARS)  ·  window  moving average (yr)
#   t_max ceiling of the residence time  ·  last_half half-residence of the last class
#FLUX_AGE = dict(rate="ac12", window=15, t_max=200.0, last_half=12.0)

# Management intensity pie chart (None to disable it).
MANAGEMENT = dict(nc="FOREST_MANAGED",
                  labels=["1 · unmanaged", "2 · extensive", "3 · moderate",
                          "4 · intensive", "5 · very intensive"])


# ══════════════════════════════════════════════════════════════════════════════
#  4. MAP  —  SPATIAL layers
# ══════════════════════════════════════════════════════════════════════════════
# One entry = one chip above the map. The order is the display order,
# the first one is active on opening.
#   value   what is mapped:
#             "dia"                 the pixel aggregate of variable dia
#             "area_c{last}/area"   a ratio between two series
#             "{last}" = last age class, "{first}" = first one
#   group   heading of the block of chips
#   label / unit / scale : inherited from the variable if omitted
#   ramp    green · carbon · gold · fire · storm · pest · harv · edge
#           balance = diverging (to be used with diverging=True)
#   dec     decimals displayed       ·  zero True = scale forced to start at 0

MAP_LAYERS = [
  #dict(value="age",              group="Age structure", ramp="green", dec=0),
  #dict(value="age_flux",         group="Age structure", ramp="green", dec=0,
  #     label="Flux age Σ1/f̄", unit="yr"),
  #dict(value="aga",              group="Age structure", ramp="green", dec=0,
  #     label="Carried age (area)"),
  #dict(value="area_c{last}/area", group="Age structure", ramp="green", dec=0,
  #     label="Mature class share", unit="% forest", scale=100, zero=True),
  dict(value="dia",              group="Age structure", ramp="green", dec=0,
       label="Mean diameter"),
  dict(value="hgt",              group="Age structure", ramp="green", dec=1,
       label="Mean height"),
  dict(value="area",             group="Age structure", ramp="green", dec=0, zero=True),
  dict(value="biomass",          group="Age structure", ramp="gold",  dec=1, zero=True),

  dict(value="da_fire",  group="Disturbance & harvest", ramp="fire",  dec=2, zero=True),
  dict(value="da_storm", group="Disturbance & harvest", ramp="storm", dec=2, zero=True),
  dict(value="da_pest",  group="Disturbance & harvest", ramp="pest",  dec=2, zero=True),
  dict(value="da_harv",  group="Disturbance & harvest", ramp="harv",  dec=2, zero=True),

  #dict(value="edge",     group="Edge", ramp="edge", dec=0, zero=True),

  dict(value="gpp", group="Carbon & productivity", ramp="carbon", dec=0, zero=True),
  dict(value="npp", group="Carbon & productivity", ramp="carbon", dec=0, zero=True),
  dict(value="ra",  group="Carbon & productivity", ramp="carbon", dec=0, zero=True),
  dict(value="rh",  group="Carbon & productivity", ramp="carbon", dec=0, zero=True),
  dict(value="lai", group="Carbon & productivity", ramp="green",  dec=1, zero=True),
  # signed quantity: symmetric scale and diverging ramp with a neutral centre,
  # otherwise a sink and a source look alike.
  dict(value="rhsoil", group="Carbon & productivity", ramp="carbon", dec=0, zero=True),
  dict(value="nbp", group="Carbon & productivity", ramp="balance", dec=0, diverging=True),
  # positive NEE = source, positive NBP = sink: the ramp is reversed so that
  # brown means "the grid cell loses carbon" on both maps.
  dict(value="nee", group="Carbon & productivity", ramp="balance", dec=0,
       diverging=True, reverse=True),
]


# ══════════════════════════════════════════════════════════════════════════════
#  5. PIXEL PLATE  —  TEMPORAL panels
# ══════════════════════════════════════════════════════════════════════════════
# One entry = one chart of the plate opened by clicking on a grid cell.
#   kind "stack"    stacked areas by age class
#        "classes"  one curve per age class (+ pixel curve if pixel=True)
#        "series"   one curve per named series
#        "flux_age" one curve per class, flux age (requires FLUX_AGE)
#   var     variable key (kinds stack/classes)
#   series  list of keys (kind series)
#   dec     decimals · ymin 0 to anchor the axis · title/unit inherited if omitted

PANELS = [
  dict(kind="stack",    var="area", title="Areas by age class", unit="% of grid cell",
       dec=0),          # extra="<key>" would add a series to the stack
  dict(kind="classes",  var="dia",  dec=0),
  dict(kind="classes",  var="hgt",  dec=1),
  #dict(kind="classes",  var="age",  dec=0, pixel=True),
  #dict(kind="classes",  var="aga",  dec=0, pixel=True),
  #dict(kind="flux_age", title="Flux age (Σ1/f̄)", unit="yr", dec=0, pixel=True),
  dict(kind="classes",  var="rdi",  dec=2),
  dict(kind="classes",  var="ind",  dec=0, per_essence=False),  # not stored per species
  dict(kind="series",   title="Disturbance & harvest", unit="% land/yr", dec=2, ymin=0,
       series=["da_fire", "da_storm", "da_pest", "da_harv"]),
  #dict(kind="series",   title="Edge length", unit="km", dec=0, ymin=0,
  #     series=["edge"]),
  dict(kind="series",   title="Living biomass (forest)", unit="kgC/m²", dec=1, ymin=0,
       series=["biomass"]),
  dict(kind="classes",  var="gpp", dec=0, pixel=True),
  dict(kind="classes",  var="npp", dec=0, pixel=True),
  dict(kind="classes",  var="lai", dec=1, pixel=True),
  # same denominator (m² of forest): the fluxes are comparable on a single axis.
  dict(kind="series",   title="Carbon fluxes (forest)", unit="gC/m² forest/yr", dec=0,
       series=["gpp", "npp", "ra", "rh"]),
  # same unit and denominator (the whole grid cell): comparable on a single axis.
  dict(kind="series",   title="Net balance of the grid cell", unit="gC/m²/yr", dec=0,
       series=["nbp", "nee"]),
  dict(kind="series",   title="Soil respiration (grid cell)", unit="gC/m²/yr", dec=0,
       ymin=0, series=["rhsoil"]),
]
