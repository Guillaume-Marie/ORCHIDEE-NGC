#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Annual cube extraction for the atlas — generic engine driven by atlas_config.py.

Nothing is hard-coded here: the number of PFTs, age classes and species,
resolution and domain come from the files; the variables, map layers and
panels come from the configuration.

Outputs (in OUT):
  atlas_meta.json   metadata + per-pixel statics + map/plate catalogue
  atlas_cube.b64    quantised int16, gzip, base64 — dims (nseries, npix, nyear)
  atlas_ess.b64     int8 per species — [pixel][species][var][class][year]
"""
import glob, gzip, base64, json, os, re, sys
import numpy as np
import netCDF4

import atlas_config as C

NUL = 1e-9


def die(msg):
    sys.exit("configuration error: " + msg)


# Optional settings of atlas_config.py may be commented out: they are read
# with a default value rather than crashing on an AttributeError.
DEFAULTS = dict(OUT=".", YEARS=None, DOMAIN=None, TITLE=None, DROP_EMPTY=True,
                FILES=dict(main="1M_stomate_history"),
                PFT_MAPPING="auto", ESSENCE_LABELS={}, ESSENCE_MIN_SHARE=0.03,
                ESSENCE_MAX_PER_PIX=4, FLUX_AGE=None, MANAGEMENT=None,
                PFT_VARS=[], GRID_VARS=[], PANELS=[],
                AREA_LABEL=dict(label="Forest area", unit="% of grid cell", scale=100))


def cfg(name):
    if hasattr(C, name):
        return getattr(C, name)
    if name in DEFAULTS:
        return DEFAULTS[name]
    die(f"{name} is missing from atlas_config.py")


DIR = cfg("SIM_DIR").rstrip("/") + "/"
OUT = cfg("OUT").rstrip("/") + "/"


# ═══════════════════════════════════════════════════════════ files & years
# A simulation writes several file families per year (stomate_history,
# stomate_ipcc_history…). The "main" family carries the PFT structure and the grid;
# the others only contribute variables (`file` key of atlas_config.py).
FAM = dict(cfg("FILES"))
if "main" not in FAM:
    die('FILES must contain a "main" family (the one carrying VEGET_MAX)')


def family_files(suffix):
    return {int(os.path.basename(p).split("_")[1][:4]): p
            for p in sorted(glob.glob(f"{DIR}{cfg('RUN')}_*_{suffix}.nc"))}


FAMFILES = {k: family_files(v) for k, v in FAM.items()}
files = sorted(FAMFILES["main"].values())
if not files:
    die(f"no file {cfg('RUN')}_*_{FAM['main']}.nc in {DIR}")


def year_of(path):
    return int(os.path.basename(path).split("_")[1][:4])


if cfg("YEARS"):
    y0, y1 = cfg("YEARS")
    files = [f for f in files if y0 <= year_of(f) <= y1]
    if not files:
        die(f"no year between {y0} and {y1} in {DIR}")

# A truncated file (interrupted archive) would make the pass fail halfway:
# unreadable years are discarded up front, with an explicit warning.
ok = []
for f in files:
    try:
        with netCDF4.Dataset(f) as ds:
            ds.variables["VEGET_MAX"][0, 0]
        ok.append(f)
    except Exception as e:
        print(f"  WARNING year skipped (unreadable file) : {os.path.basename(f)}"
              f" [{type(e).__name__}]", flush=True)
files, nyr = ok, len(ok)
if not files:
    die(f"no readable {FAM['main']} file")
years = [year_of(f) for f in files]
print(f"{nyr} years: {years[0]}–{years[-1]}", flush=True)


# ═══════════════════════════════════════════════════════ PFT structure
def pft_names(ds):
    v = ds.variables["veget"][:]
    if v.dtype.kind in "SU":                       # character axis (str_len)
        return [b"".join(np.ma.filled(r, b" ")).decode("utf-8", "replace").strip()
                for r in v]
    return [f"PFT{i + 1:02d}" for i in range(len(v))]


TREE_WORDS = ("broadleaved", "needleleaf", "larix", "eucalyptus", "tree", "forest")
NOTREE_WORDS = ("grass", "agriculture", "crop", "bare", "soil", "pasture")
# The veget axis is a fixed-width character array: long names are TRUNCATED
# there (e.g. "..._age01" .. "..._age04" all become "..._age").
# A partial class suffix is therefore accepted and, without a number, classes
# are numbered in order of appearance — the PFTs of a species are contiguous.
AGE_RE = re.compile(r"^(?P<base>.+?)_ag?e?0*(?P<k>\d*)$", re.I)

with netCDF4.Dataset(files[0]) as ds:
    NPFT = ds.variables["VEGET_MAX"].shape[1]
    PFTNAMES = pft_names(ds)
    NCVARS = {"main": set(ds.variables)}
for fam, byyear in FAMFILES.items():
    if fam == "main":
        continue
    path = byyear.get(year_of(files[0])) or (sorted(byyear.values()) or [None])[0]
    try:
        with netCDF4.Dataset(path) as ds2:
            NCVARS[fam] = set(ds2.variables)
    except Exception:
        NCVARS[fam] = set()
        print(f"  WARNING family {fam} ({FAM[fam]}) not found or unreadable : "
              "its variables will be discarded", flush=True)
if len(PFTNAMES) != NPFT:
    PFTNAMES = [f"PFT{i + 1:02d}" for i in range(NPFT)]


def auto_mapping():
    """Reads the veget axis: '<Species>_ageNN' -> species × class. Without an _age
    suffix, a woody PFT counts as a single-class species."""
    order, slots, guessed, dup = [], {}, False, False   # base -> {class: PFT index}
    for i, nm in enumerate(PFTNAMES):
        m = AGE_RE.match(nm)
        base = m.group("base") if m else nm
        low = base.lower()
        if any(w in low for w in NOTREE_WORDS) or not any(w in low for w in TREE_WORDS):
            continue
        if m and m.group("k"):                  # explicit class suffix
            k = int(m.group("k"))
        elif m:                                 # truncated suffix: rank of appearance
            k = len(slots.get(base, {})) + 1
            guessed = guessed or k > 1
        else:
            # NO class suffix: a repeated name is not an age class, but two
            # distinct PFTs sharing the same MTC (NAGEC=1 config, e.g. two
            # management types of the same vegetation type). Each is one entry.
            k = 1
            if base in slots:
                dup = True
                n = 2
                while f"{base} #{n}" in slots:
                    n += 1
                base = f"{base} #{n}"
        if base not in slots:
            slots[base] = {}
            order.append(base)
        if k in slots[base]:
            die(f"ambiguous veget axis: class {k} of {base!r} appears twice "
                "(PFT {} and {}) — set PFT_MAPPING explicitly".format(slots[base][k], i))
        slots[base][k] = i
    if guessed:
        print("  note: PFT names truncated by the file — classes inferred from "
              "the order of the veget axis", flush=True)
    if dup:
        print("  note: some woody PFTs share a name with no class suffix "
              "(_ageNN) — counted as distinct types, suffixed #2, #3…",
              flush=True)
    if not order:
        die("no woody PFT recognised on the veget axis "
            f"({', '.join(PFTNAMES[:4])}…) — set PFT_MAPPING explicitly")
    return order, slots


def explicit_mapping(cfg):
    names, blocks = cfg["essences"], cfg["blocks"]
    if len(names) != len(blocks):
        die("PFT_MAPPING: essences and blocks have different lengths")
    slots = {}
    for nm, (a, b) in zip(names, blocks):
        if b > NPFT:
            die(f"PFT_MAPPING: block {(a, b)} overruns the veget axis ({NPFT} PFT)")
        slots[nm] = {k + 1: p for k, p in enumerate(range(a, b))}
    return list(names), slots


base_order, slots = (auto_mapping() if cfg("PFT_MAPPING") == "auto"
                     else explicit_mapping(cfg("PFT_MAPPING")))
NESS = len(base_order)
NCLS = max(max(d) for d in slots.values())
def esslabel(b):
    """Presentation: "Name #2" reuses the label of "Name", suffix kept."""
    lab = cfg("ESSENCE_LABELS")
    if b in lab:
        return lab[b]
    root, sep, tag = b.partition(" #")
    return f"{lab[root]} #{tag}" if sep and root in lab else b


ESSENCES = [esslabel(b) for b in base_order]
# PFT_OF[g][k] = index of the PFT of species g, class k+1 (-1 if the class is missing)
PFT_OF = [[slots[b].get(k + 1, -1) for k in range(NCLS)] for b in base_order]
TREE_IDX = sorted(p for row in PFT_OF for p in row if p >= 0)
tpos = {p: i for i, p in enumerate(TREE_IDX)}
NTREE = len(TREE_IDX)
# class of each tree PFT, in TREE_IDX order
tcls = np.zeros(NTREE, dtype=int)
for g, row in enumerate(PFT_OF):
    for k, p in enumerate(row):
        if p >= 0:
            tcls[tpos[p]] = k + 1
print(f"{NPFT} PFT → {NESS} species × {NCLS} classes ({NTREE} tree PFTs)", flush=True)
counts = {len([p for p in row if p >= 0]) for row in PFT_OF}
if len(counts) > 1:
    print(f"  WARNING species of unequal size ({sorted(counts)} classes) : "
          "check the mapping below, or fix it via PFT_MAPPING.", flush=True)
print("  " + " · ".join(f"{e}[{sum(1 for p in row if p >= 0)}]"
                        for e, row in zip(ESSENCES, PFT_OF)), flush=True)


# ═══════════════════════════════════════════════════════ series catalogue
def ncl(v):
    """NetCDF names of a variable: `nc` is one name, or a list to sum."""
    n = v["nc"]
    return [n] if isinstance(n, str) else list(n)


def nclab(v):
    return " + ".join(ncl(v))


def norm(v, kind):
    d = dict(v)
    d.setdefault("file", "main")
    if d["file"] not in FAM:
        die(f"{d['key']} : file={d['file']!r} is not a FILES family "
            f"({', '.join(FAM)})")
    d.setdefault("scale", 1.0)
    d.setdefault("label", d["key"])
    d.setdefault("unit", "")
    d.setdefault("clip", None)
    if kind == "pft":
        d.setdefault("agg", "wmean")
        d.setdefault("time", "mean")
        d.setdefault("per", "day")
        if d["time"] not in ("mean", "sum", "max", "first", "annual"):
            die(f"{d['key']} : time must be mean, sum, max, first or annual")
        d.setdefault("classes", True)
        d.setdefault("essence", False)
        d.setdefault("step", None)
        if d["agg"] not in ("wmean", "wsum"):
            die(f"{d['key']} : agg must be wmean or wsum")
    else:
        d.setdefault("agg", "mean")
        d.setdefault("per", "day")
        if d["agg"] not in ("mean", "sum", "rate", "annual"):
            die(f"{d['key']} : agg must be mean, sum, rate or annual")
    return d


PV = [norm(v, "pft") for v in cfg("PFT_VARS")]
GV = [norm(v, "grid") for v in cfg("GRID_VARS")]
AREA = dict(key="area", classes=True, clip=None, **cfg("AREA_LABEL"))

VARS = {v["key"]: v for v in [AREA] + PV + GV}
if len(VARS) != 1 + len(PV) + len(GV):
    die("two variables share the same key")

SERIES, CLIP = [], {}


def add_series(name, clip=None):
    SERIES.append(name)
    if clip is not None:
        CLIP[name] = clip


for v in [AREA] + PV:
    if v.get("classes", True):
        for k in range(1, NCLS + 1):
            add_series(f"{v['key']}_c{k}", v["clip"])
    add_series(v["key"], v["clip"])
for v in GV:
    add_series(v["key"], v["clip"])

FLUX = dict(cfg("FLUX_AGE")) if cfg("FLUX_AGE") else None
if FLUX:
    if FLUX["rate"] not in {v["key"] for v in PV}:
        die(f"FLUX_AGE.rate = {FLUX['rate']!r} is not a PFT_VARS key")
    add_series("age_flux")
SI = {n: i for i, n in enumerate(SERIES)}
NS = len(SERIES)

ESSV = [v for v in PV if v["essence"]]
for v in ESSV:
    if not v["step"]:
        die(f"{v['key']} : essence=True requires a quantisation step 'step'")
# Per-species blob: the area always opens the list (it serves as weight and
# mask for the other variables), hence EB = area + the essence=True variables.
EB = [dict(AREA, step=1 / 250)] + ESSV
ESS_STEPS = [v["step"] for v in EB]
epos = {v["key"]: i for i, v in enumerate(EB)}   # position in the per-species blob
erow = dict(epos)                                # corresponding row of pftann
if FLUX and FLUX["rate"] not in epos:
    print(f"  note: {FLUX['rate']} is not stored per species - the flux age "
          "of the plate will stay at pixel level", flush=True)

OTHER_FAMS = sorted({v["file"] for v in PV + GV} - {"main"})
MGT = cfg("MANAGEMENT")
MISSING = sorted({n for v in PV + GV for n in ncl(v) if n not in NCVARS[v["file"]]}
                 | ({MGT["nc"]} if MGT and MGT["nc"] not in NCVARS["main"] else set()))
if MISSING:
    print("fields absent from the run (empty panels) : " + ", ".join(MISSING), flush=True)


# ═══════════════════════════════════════════ map & plate catalogue
LASTC, FIRSTC = str(NCLS), "1"


def resolve(expr):
    """'dia' or 'area_c{last}/area' -> (numerator, denominator|None)."""
    e = expr.replace("{last}", LASTC).replace("{first}", FIRSTC)
    parts = [p.strip() for p in e.split("/")]
    if len(parts) > 2:
        die(f"expression {expr!r} : at most one division")
    for p in parts:
        if p not in SI:
            die(f"expression {expr!r} : unknown series {p!r}. "
                f"Available series : {', '.join(SERIES)}")
    return parts[0], (parts[1] if len(parts) == 2 else None)


def base_var(name):
    return VARS.get(re.sub(r"_c\d+$", "", name), {})


layers = []
for L in cfg("MAP_LAYERS"):
    num, den = resolve(L["value"])
    bv = base_var(num) if den is None else {}
    layers.append({
        "key": L.get("key", L["value"]),
        "group": L.get("group", ""),
        "label": L.get("label", bv.get("label", num)),
        "unit": L.get("unit", bv.get("unit", "")),
        "scale": float(L.get("scale", bv.get("scale", 1.0))),
        "ramp": L.get("ramp", "green"),
        "dec": int(L.get("dec", 1)),
        "zero": bool(L.get("zero", False)),
        # signed quantity (sink / source): scale symmetric around 0
        "div": bool(L.get("diverging", False)),
        # reverses the ramp: useful when two neighbouring layers follow
        # opposite sign conventions (positive NBP = sink, positive NEE = source)
        "rev": bool(L.get("reverse", False)),
        "num": num, "den": den,
    })
if not layers:
    die("MAP_LAYERS is empty: at least one map layer is required")

panels = []
for P in cfg("PANELS"):
    kind = P.get("kind", "classes")
    if kind not in ("stack", "classes", "series", "flux_age"):
        die(f"panel {P.get('title', P.get('var'))!r} : unknown kind {kind!r}")
    if kind == "flux_age" and not FLUX:
        continue
    v = {}
    if kind in ("stack", "classes"):
        if P.get("var") not in VARS:
            die(f"panel : var {P.get('var')!r} is not a known variable "
                f"({', '.join(VARS)})")
        v = VARS[P["var"]]
        if not v.get("classes", True):
            die(f"panel {P['var']!r} : kind={kind} requires a per-class variable "
                "(classes=True)")
    for s in P.get("series", []):
        if s not in SI:
            die(f"panel {P.get('title')!r} : unknown series {s!r}")
    panels.append({
        "kind": kind,
        "title": P.get("title", v.get("label", "")),
        "unit": P.get("unit", v.get("unit", "")),
        "dec": int(P.get("dec", 1)),
        "ymin": P.get("ymin", None),
        "var": P.get("var"),
        "scale": float(P.get("scale", v.get("scale", 1.0))),
        "series": [{"key": s, "label": VARS[s]["label"] if s in VARS else s,
                    "scale": float(VARS[s]["scale"]) if s in VARS else 1.0}
                   for s in P.get("series", [])],
        "pixel": bool(P.get("pixel", False)),
        "extra": P.get("extra"),
        # the species filter applies if the variable is in the int8 blob;
        # for the flux age, it is the realised rate that must be there.
        "ess": bool(P.get("per_essence", True)) and (
            (FLUX and FLUX["rate"] in epos) if kind == "flux_age" else P.get("var") in epos),
        "essvar": epos.get(P.get("var")),
    })


# The catalogue is resolved BEFORE reading the files: a typo in
# atlas_config.py must cost one second, not a full extraction.
# ═══════════════════════════════════════════════════════ grid & domain
with netCDF4.Dataset(files[0]) as ds:
    lat = np.ma.filled(ds.variables["lat"][:].astype("f8"), np.nan)
    lon = np.ma.filled(ds.variables["lon"][:].astype("f8"), np.nan)
    areas = np.ma.filled(ds.variables["Areas"][:].astype("f8"), np.nan)
    contfrac = np.ma.filled(ds.variables["CONTFRAC"][:].astype("f8"), np.nan)
    if areas.ndim == 3:
        areas = areas[0]
    if contfrac.ndim == 3:
        contfrac = contfrac[0]

jsel = np.arange(len(lat))
isel = np.arange(len(lon))
if cfg("DOMAIN"):
    if "lat" in cfg("DOMAIN"):
        a, b = cfg("DOMAIN")["lat"]; jsel = jsel[(lat >= a) & (lat <= b)]
    if "lon" in cfg("DOMAIN"):
        a, b = cfg("DOMAIN")["lon"]; isel = isel[(lon >= a) & (lon <= b)]
    if not len(jsel) or not len(isel):
        die(f"DOMAIN {cfg("DOMAIN")} does not intersect the grid "
            f"(lat {lat.min():g}..{lat.max():g}, lon {lon.min():g}..{lon.max():g})")
    lat, lon = lat[jsel], lon[isel]
    areas = areas[np.ix_(jsel, isel)]
    contfrac = contfrac[np.ix_(jsel, isel)]
CROP = (slice(jsel[0], jsel[-1] + 1), slice(isel[0], isel[-1] + 1))
nlat, nlon = len(lat), len(lon)
dlat = float(abs(np.diff(lat)).mean()) if nlat > 1 else 1.0
dlon = float(abs(np.diff(lon)).mean()) if nlon > 1 else 1.0
land_area = areas * contfrac                     # m² of land per grid cell
land = contfrac > 0.05
lj, li = np.where(land)
NL = len(lj)
print(f"grid {nlat}x{nlon} ({dlat:g}°×{dlon:g}°), {NL} land cells", flush=True)

# All computation is done on the land cells only, flattened: this is what
# allows going up in resolution without blowing up the memory.
cube = np.full((NS, nyr, NL), np.nan, dtype="f4")
pftann = np.full((len(EB), nyr, NTREE, NL), np.nan, dtype="f4")
ess_acc = np.zeros((NESS, NL))
ess_n = 0
la_land = land_area[lj, li]


def rd(ds, name, default=None):
    """Reads a field, or the sum of several (`name` may be a list).

    A sum is only computed if ALL its terms exist: a total missing one
    term would be wrong, and would go unnoticed."""
    if ds is None:
        return default
    names = [name] if isinstance(name, str) else list(name)
    if any(n not in ds.variables for n in names):
        return default
    arrs = [np.ma.filled(ds.variables[n][:].astype("f4"), np.nan) for n in names]
    if len(arrs) == 1:
        return arrs[0]
    tot = np.zeros_like(arrs[0])
    fin = np.zeros(arrs[0].shape, dtype=bool)
    for a in arrs:
        f = np.isfinite(a)
        tot += np.where(f, a, 0)
        fin |= f
    return np.where(fin, tot, np.nan)


def tagg(x, how, days=None, per="day"):
    """(time, …) -> (…) according to the requested temporal aggregation.

    "annual" : annual total of a RATE given per day (or per second, per="s"),
    i.e. Σ months (mean monthly value × days in the month). This is what is
    needed for GPP, NPP, respiration, NBP — a mean over the months would give a
    daily rate, not an annual cumulative."""
    with np.errstate(invalid="ignore"):
        if how == "annual":
            k = days * (86400.0 if per == "s" else 1.0)
            k = k.reshape((-1,) + (1,) * (x.ndim - 1))
            return np.where(np.isfinite(x).any(axis=0), np.nansum(x * k, axis=0), np.nan)
        if how == "sum":
            return np.where(np.isfinite(x).any(axis=0), np.nansum(x, axis=0), np.nan)
        if how == "max":
            return np.nanmax(x, axis=0)
        if how == "first":
            return x[0]
        return np.nanmean(x, axis=0)


def wmean(x, w, den, denmin):
    """Weighted mean; NaN — not zero — when no value is finite
    (field absent from the simulation, empty slot)."""
    fx = np.isfinite(x)
    num = np.nansum(np.where(fx, x, 0) * w, axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where((den > denmin) & fx.any(axis=0), num / np.maximum(den, NUL), np.nan)


AMIN = 1e-4
fm = None
for iy, f in enumerate(files):
    with netCDF4.Dataset(f) as ds:
        nt = ds.variables["VEGET_MAX"].shape[0]
        DAYS = (np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], "f8")
                if nt == 12 else np.full(nt, 365.0 / nt))
        vm = rd(ds, "VEGET_MAX")[(slice(None), slice(None)) + CROP]
        vma = tagg(vm, "mean")[:, lj, li]                       # (NPFT, NL)
        nanp = np.full((NPFT, NL), np.nan, "f4")
        st = {}
        for v in PV:
            if v["file"] != "main":
                continue
            x = rd(ds, v["nc"])
            st[v["key"]] = (tagg(x[(slice(None), slice(None)) + CROP], v["time"],
                                 DAYS, v["per"])[:, lj, li] if x is not None else nanp)
        d2 = {}
        for v in GV:
            if v["file"] != "main":
                continue
            x = rd(ds, v["nc"])
            d2[v["key"]] = (x[(slice(None),) + CROP][:, lj, li]
                            if x is not None else np.full((nt, NL), np.nan, "f4"))
        if iy == nyr - 1 and MGT:
            x = rd(ds, MGT["nc"])
            fm = (tagg(x[(slice(None), slice(None)) + CROP], "mean")[:, lj, li]
                  if x is not None else nanp)

    # secondary families: same grid and calendar, additional variables.
    # A year missing from a family leaves its variables at NaN for that year.
    for fam in OTHER_FAMS:
        path = FAMFILES[fam].get(years[iy])
        ds2 = None
        try:
            ds2 = netCDF4.Dataset(path) if path else None
        except Exception:
            ds2 = None
        if ds2 is None and iy == 0:
            print(f"  WARNING {FAM[fam]} missing for {years[iy]}", flush=True)
        for v in PV:
            if v["file"] != fam:
                continue
            x = rd(ds2, v["nc"]) if ds2 else None
            st[v["key"]] = (tagg(x[(slice(None), slice(None)) + CROP], v["time"],
                                 DAYS, v["per"])[:, lj, li] if x is not None else nanp)
        for v in GV:
            if v["file"] != fam:
                continue
            x = rd(ds2, v["nc"]) if ds2 else None
            d2[v["key"]] = (x[(slice(None),) + CROP][:, lj, li]
                            if x is not None else np.full((nt, NL), np.nan, "f4"))
        if ds2 is not None:
            ds2.close()

    w = np.where(np.isfinite(vma), vma, 0.0)                    # weight = PFT area
    wt = w[TREE_IDX]
    pftann[epos["area"], iy] = wt
    for v in ESSV:
        pftann[epos[v["key"]], iy] = st[v["key"]][TREE_IDX]

    # --- per age class, then pixel aggregate (all classes)
    fa = wt.sum(axis=0)
    cube[SI["area"], iy] = fa
    for k in range(1, NCLS + 1):
        sel = tcls == k
        wk = wt[sel]
        ak = wk.sum(axis=0)
        cube[SI[f"area_c{k}"], iy] = ak
        for v in PV:
            if not v["classes"]:
                continue
            x = st[v["key"]][TREE_IDX][sel]
            if v["agg"] == "wmean":
                cube[SI[f"{v['key']}_c{k}"], iy] = wmean(x, wk, ak, AMIN)
            else:
                fx = np.isfinite(x)
                cube[SI[f"{v['key']}_c{k}"], iy] = np.where(
                    fx.any(axis=0), np.nansum(np.where(fx, x, 0) * wk, axis=0), np.nan)
    for v in PV:
        x = st[v["key"]][TREE_IDX]
        if v["agg"] == "wmean":
            cube[SI[v["key"]], iy] = wmean(x, wt, fa, AMIN)
        else:
            fx = np.isfinite(x)
            cube[SI[v["key"]], iy] = np.where(
                fx.any(axis=0), np.nansum(np.where(fx, x, 0) * wt, axis=0), np.nan)

    # --- grid fields
    for v in GV:
        x = d2[v["key"]]
        fin = np.isfinite(x).any(axis=0)
        if v["agg"] == "rate":
            # m²/day accumulator averaged over the month -> annual area / land area
            tot = np.nansum(x * DAYS[:, None], axis=0)
            with np.errstate(invalid="ignore", divide="ignore"):
                cube[SI[v["key"]], iy] = np.where(la_land > 0, tot / np.maximum(la_land, NUL),
                                                  np.nan)
            cube[SI[v["key"]], iy] = np.where(fin, cube[SI[v["key"]], iy], np.nan)
        else:
            cube[SI[v["key"]], iy] = tagg(x, v["agg"] if v["agg"] in ("sum", "annual")
                                          else "mean", DAYS, v["per"])

    if iy >= nyr - 5:
        for g, row in enumerate(PFT_OF):
            ess_acc[g] += w[[p for p in row if p >= 0]].sum(axis=0)
        ess_n += 1
    if iy % 10 == 0:
        print(f"  {iy + 1}/{nyr}", flush=True)

# ═══════════════════════════════════════════ flux age Σ1/f̄ (optional)
if FLUX:
    # Realised rates smoothed over `window` years ZEROS INCLUDED; T capped at t_max
    # when there is no transfer; terminal class = entry age + last_half.
    r = pftann[epos[FLUX["rate"]]] if FLUX["rate"] in epos else None
    if r is None:
        cube[SI["age_flux"]] = np.nan
    else:
        W, TMAX, THALF = int(FLUX["window"]), float(FLUX["t_max"]), float(FLUX["last_half"])
        num = np.zeros((nyr, NL)); den = np.zeros((nyr, NL))
        for row in PFT_OF:
            entry = np.zeros((nyr, NL))
            for k, p in enumerate(row):
                if p < 0:
                    continue
                fr = np.nan_to_num(r[:, tpos[p]])
                fbar = np.empty_like(fr)
                for t in range(nyr):
                    fbar[t] = fr[max(0, t - W + 1):t + 1].mean(axis=0)
                T = np.minimum(np.where(fbar > 1.0 / TMAX, 1.0 / np.maximum(fbar, NUL), TMAX),
                               TMAX)
                age = entry + (THALF if k == NCLS - 1 else T / 2.0)
                wp = np.nan_to_num(pftann[epos["area"], :, tpos[p]])
                num += wp * age
                den += wp
                if k < NCLS - 1:
                    entry = entry + T
        # without a realised rate in the simulation, T saturates at the cap: the
        # resulting "age" would only be that cap — NaN is returned to leave the panel empty.
        has = np.isfinite(r).any(axis=(0, 1))
        with np.errstate(invalid="ignore", divide="ignore"):
            cube[SI["age_flux"]] = np.where(has & (den > AMIN), num / np.maximum(den, NUL),
                                            np.nan)

# ═══════════════════════════════════════════ pixel selection
fmax = np.nanmax(cube[SI["area"]], axis=0)
keep = np.isfinite(fmax) & (fmax > 1e-3)
sel = np.where(keep)[0]
npix = len(sel)
py, px = lj[sel], li[sel]
if not npix:
    die("no pixel carries forest: check SIM_DIR, DOMAIN and the PFT mapping")
print(f"{npix} pixels kept (land + forest) out of {NL} land cells", flush=True)

# ═══════════════════════════════════════════ species & management per pixel
ess = ess_acc[:, sel] / max(ess_n, 1)
ess_share = ess / np.maximum(ess.sum(axis=0), NUL)
pix_ess = []
for i in range(npix):
    sh = ess_share[:, i]
    pix_ess.append([int(g) for g in np.argsort(-sh)
                    if sh[g] >= cfg("ESSENCE_MIN_SHARE")][:cfg("ESSENCE_MAX_PER_PIX")])

fm_share, MI_LABELS = [], []
if MGT and fm is not None and np.isfinite(fm).any():
    MI_LABELS = list(MGT["labels"])
    wlast = pftann[epos["area"], -1][:, sel]                       # tree PFT areas, last year
    cls = np.clip(np.rint(np.where(np.isfinite(fm[TREE_IDX][:, sel]), fm[TREE_IDX][:, sel], 0)),
                  1, len(MI_LABELS)).astype(int)
    sh = np.stack([np.where(cls == c, wlast, 0).sum(axis=0)
                   for c in range(1, len(MI_LABELS) + 1)])
    fm_share = np.round(sh / np.maximum(sh.sum(axis=0), NUL), 4).T.tolist()

# ═══════════════════════════════════════════ pruning of empty variables
# A field absent from the simulation — or present but without a single finite
# value over the domain — has nothing to show. Rather than greyed-out chips and
# "no data" panels, it is discarded: series, layers, panels and its place
# in the per-species blob. DROP_EMPTY = False restores the empty slots.
DROPPED = []
if cfg("DROP_EMPTY"):
    live = {n for i, n in enumerate(SERIES) if np.isfinite(cube[i][:, sel]).any()}
    live |= {"area"} | {f"area_c{k}" for k in range(1, NCLS + 1)}   # structure: never pruned
    gone = [n for n in SERIES if n not in live]
    if gone:
        cube = cube[[SI[n] for n in SERIES if n in live]]
        SERIES = [n for n in SERIES if n in live]
        SI = {n: i for i, n in enumerate(SERIES)}
        NS = len(SERIES)

    # variables with no surviving series left
    def alive(key):
        return key in SI or any(f"{key}_c{k}" in SI for k in range(1, NCLS + 1))

    DROPPED = sorted({nclab(v) for v in PV + GV if not alive(v["key"])})

    # per-species blob: only what carries values is kept
    EB = [v for v in EB if v["key"] == "area"
          or (alive(v["key"]) and np.isfinite(pftann[epos[v["key"]]]).any())]
    erow = {v["key"]: epos[v["key"]] for v in EB}      # row of pftann, unchanged
    epos = {v["key"]: i for i, v in enumerate(EB)}     # position in the blob
    ESS_STEPS = [v["step"] for v in EB]

    layers = [L for L in layers if L["num"] in SI and (L["den"] is None or L["den"] in SI)]
    if not layers:
        die("no map layer carries data: check SIM_DIR and MAP_LAYERS")
    kept = []
    for P in panels:
        if P["kind"] == "series":
            P["series"] = [x for x in P["series"] if x["key"] in SI]
            if not P["series"]:
                continue
        elif P["kind"] == "flux_age":
            if "age_flux" not in SI:
                continue
        else:
            if not alive(P["var"]):
                continue
            P["pixel"] = P["pixel"] and P["var"] in SI
            if P["extra"] and P["extra"] not in SI:
                P["extra"] = None
        P["ess"] = P["ess"] and (
            (FLUX and FLUX["rate"] in epos) if P["kind"] == "flux_age" else P["var"] in epos)
        P["essvar"] = epos.get(P["var"])
        kept.append(P)
    panels = kept
    if FLUX:
        FLUX_ESS = epos.get(FLUX["rate"], -1)
    if DROPPED or gone:
        print(f"pruned: {len(gone)} empty series"
              + (f" ({', '.join(DROPPED)})" if DROPPED else "")
              + f" — remaining {NS} series, {len(layers)} layers, {len(panels)} panels",
              flush=True)
else:
    FLUX_ESS = epos.get(FLUX["rate"], -1) if FLUX else -1


# ═══════════════════════════════════════════ int16 quantisation
data = np.transpose(cube[:, :, sel], (0, 2, 1)).copy()          # (NS, npix, nyr)
q = np.zeros(data.shape, dtype="<i2")
meta_s = []
for s, name in enumerate(SERIES):
    a = data[s]
    if name in CLIP and CLIP[name] is not None:
        a = np.minimum(a, CLIP[name])
    fin = np.isfinite(a)
    lo = float(np.nanmin(a)) if fin.any() else 0.0
    hi = float(np.nanmax(a)) if fin.any() else 1.0
    if hi <= lo:
        hi = lo + 1.0
    sc = (hi - lo) / 32000.0
    q[s] = np.where(fin, np.rint((a - lo) / sc), -32768).astype("<i2")
    meta_s.append({"name": name, "off": lo, "scale": sc})
b64 = base64.b64encode(gzip.compress(q.tobytes(), 9)).decode()

# ═══════════════════════════════════════════ per-species blob (int8)
# Layout: [pixel][kept species][variable][class][year] — each pair
# (species, class) is a single PFT, hence raw series, without weighting.
segs = []
for i, p in enumerate(sel):
    for g in pix_ess[i]:
        for v in EB:
            stp = v["step"]
            row = erow[v["key"]]
            for k in range(NCLS):
                pf = PFT_OF[g][k]
                arr = (pftann[row, :, tpos[pf], p] if pf >= 0
                       else np.full(nyr, np.nan, "f4"))
                segs.append(np.where(np.isfinite(arr),
                                     np.clip(np.rint(arr / stp), 0, 250) - 124,
                                     -128).astype("<i1"))
eraw = np.concatenate(segs) if segs else np.zeros(0, dtype="<i1")
eb64 = base64.b64encode(gzip.compress(eraw.tobytes(), 9)).decode()

meta = {
    "run": cfg("RUN"),
    "title": cfg("TITLE") or f"ORCHIDEE Atlas — {cfg("RUN")}",
    "res": (f"{dlat:g}°" if abs(dlat - dlon) < 1e-6 else f"{dlat:g}°×{dlon:g}°"),
    "missing": MISSING, "dropped": DROPPED,
    "years": years, "nyr": nyr, "npix": npix, "nseries": NS,
    "ncls": NCLS, "series": meta_s, "essences": ESSENCES,
    "vars": {k: {"label": v["label"], "unit": v["unit"], "scale": float(v["scale"])}
             for k, v in VARS.items()
             if k in SI or any(f"{k}_c{j}" in SI for j in range(1, NCLS + 1))},
    "layers": layers, "panels": panels,
    "ess_vars": [{"key": v["key"], "label": v["label"]} for v in EB],
    "ess_steps": ESS_STEPS,
    "flux": ({"window": FLUX["window"], "t_max": FLUX["t_max"],
              "last_half": FLUX["last_half"], "ess": FLUX_ESS}
             if FLUX and any(P["kind"] == "flux_age" for P in panels) else None),
    "mi_labels": MI_LABELS,
    "lat": [round(float(lat[j]), 4) for j in py],
    "lon": [round(float(lon[i]), 4) for i in px],
    "iy": py.tolist(), "ix": px.tolist(),
    "grid": {"lat": [round(float(x), 4) for x in lat],
             "lon": [round(float(x), 4) for x in lon],
             "dlat": round(dlat, 6), "dlon": round(dlon, 6),
             "lat0": round(float(lat.mean()), 3), "lon0": round(float(lon.mean()), 3)},
    "ess_share": np.round(ess_share, 4).T.tolist(),
    "fm_share": fm_share,
    "land_area_km2": [round(float(land_area[j, i]) / 1e6, 1) for j, i in zip(py, px)],
    "bare_iy": lj[~keep].tolist(), "bare_ix": li[~keep].tolist(),
    "pix_ess": pix_ess,
}
with open(OUT + "atlas_meta.json", "w") as f:
    json.dump(meta, f)
with open(OUT + "atlas_cube.b64", "w") as f:
    f.write(b64)
with open(OUT + "atlas_ess.b64", "w") as f:
    f.write(eb64)
print(f"OK: {NS} series × {npix} pixels × {nyr} years — cube {len(b64) / 1e6:.2f} MB b64, "
      f"species {len(eb64) / 1e6:.2f} MB b64", flush=True)
