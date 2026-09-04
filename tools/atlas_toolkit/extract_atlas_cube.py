#!/usr/bin/env python3
"""Extraction du cube annuel euage pour l'artifact carte+planches pixel.

Sortie : euage_cube_meta.json (métadonnées + statiques par pixel)
         euage_cube.b64      (int16 quantifié, gzip, base64) dims (nseries, npix, nyear)
"""
import glob, gzip, base64, json, sys
import numpy as np
import netCDF4

# ===================== CONFIGURATION (à adapter) =====================
SIM_DIR = "/chemin/vers/IGCM_OUT/.../SBG/Output/MO"   # dossier des stomate_history mensuels
RUN     = "masimu"    # préfixe des fichiers <RUN>_YYYY0101_YYYY1231_1M_stomate_history.nc (sans "_")
OUT     = "."         # dossier de sortie des blobs (meta.json + .b64)
# Le mapping PFT->essence×classe (TREE_BLOCKS/ESSENCES ci-dessous) correspond à la
# configuration 51 PFT / 4 classes d'âge ; à adapter pour une autre config (README).
# =====================================================================
DIR = SIM_DIR.rstrip("/") + "/"
OUT = OUT.rstrip("/") + "/"

files = sorted(glob.glob(DIR + RUN + "_*_1M_stomate_history.nc"))
print(f"{len(files)} fichiers", flush=True)

# --- Mapping PFT (0-based) -> (groupe essence 0..10, classe d'age 1..4)
# 11 groupes d'essences arbres, 4 classes chacun.
ESSENCES = ["Trop. BE", "Trop. BR", "Pin temp.", "Épicéa temp.", "Feuillu semp.",
            "Feuillu temp.", "Pin boréal", "Épicéa boréal", "Feuillu boréal",
            "Mélèze", "Eucalyptus"]
TREE_BLOCKS = [(1, 5), (5, 9), (9, 13), (13, 17), (17, 21), (21, 25),
               (25, 29), (29, 33), (33, 37), (37, 41), (47, 51)]  # 0-based [a,b)
pft_grp = -np.ones(51, dtype=int)
pft_cls = -np.ones(51, dtype=int)
for g, (a, b) in enumerate(TREE_BLOCKS):
    for k, p in enumerate(range(a, b)):
        pft_grp[p] = g
        pft_cls[p] = k + 1  # classe 1..4
TREE = pft_grp >= 0

DAYS = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype="f8")

VEG_VARS = ["VEGET_MAX", "DIA_INV", "HEIGHT", "AGE_STAND_BM", "RDI", "IND", "TOTAL_M_c",
            "AC12_F_REAL", "AGE_STAND_AREA"]
V2D = ["CLASS0_VEG", "FRAC_CLASS0", "DA_FIRE", "DA_STORM", "DA_PEST", "DA_HARVEST",
       "EDGE_LENGTH_DYN"]

def rd(ds, name):
    v = ds.variables[name][:]
    return np.ma.filled(v.astype("f4"), np.nan)

# --- Statique grille (premier fichier)
with netCDF4.Dataset(files[0]) as ds:
    lat = np.ma.filled(ds.variables["lat"][:].astype("f8"), np.nan)
    lon = np.ma.filled(ds.variables["lon"][:].astype("f8"), np.nan)
    areas = np.ma.filled(ds.variables["Areas"][:].astype("f8"), np.nan)      # m2
    contfrac = np.ma.filled(ds.variables["CONTFRAC"][:].astype("f8"), np.nan)
    if areas.ndim == 3: areas = areas[0]
    if contfrac.ndim == 3: contfrac = contfrac[0]

nlat, nlon = len(lat), len(lon)
nyr = len(files)
land_area = areas * contfrac  # m2 de terre par maille

SERIES = (["area_c%d" % k for k in (1, 2, 3, 4)] +
          ["dia_c%d" % k for k in (1, 2, 3, 4)] +
          ["hgt_c%d" % k for k in (1, 2, 3, 4)] +
          ["age_c%d" % k for k in (1, 2, 3, 4)] +
          ["rdi_c%d" % k for k in (1, 2, 3, 4)] +
          ["aga_c%d" % k for k in (1, 2, 3, 4)] +
          ["ind_c%d" % k for k in (1, 2, 3, 4)] +
          ["class0", "frac0", "da_fire", "da_storm", "da_pest", "da_harv",
           "edge", "forest_area", "age_porte", "biomass", "age_flux", "age_area"])
NS = len(SERIES)
cube = np.full((NS, nyr, nlat, nlon), np.nan, dtype="f4")
# état annuel par PFT arbre (pour les séries PAR ESSENCE de la planche) : ~160 Mo, ok sur p920
tree_idx = np.where(TREE)[0]
tpos = {int(p): i for i, p in enumerate(tree_idx)}
pftann = np.full((7, nyr, len(tree_idx), nlat, nlon), np.nan, dtype="f4")
years = []
# accumulateur essences (5 dernières années) et FM (dernier fichier)
ess_acc = np.zeros((11, nlat, nlon)); ess_n = 0

for iy, f in enumerate(files):
    years.append(int(f.split("/")[-1].split("_")[1][:4]))
    with netCDF4.Dataset(f) as ds:
        vm = rd(ds, "VEGET_MAX")                       # (12,51,lat,lon)
        vma = np.nanmean(vm, axis=0)                   # (51,lat,lon) moyenne annuelle
        nanv = np.full(vm.shape[1:], np.nan, "f4")               # (nvm,lat,lon)
        nan2 = np.full((vm.shape[0],) + vm.shape[2:], np.nan, "f4")  # (12,lat,lon)
        st = {v: (np.nanmean(rd(ds, v), axis=0) if v in ds.variables else nanv)
              for v in VEG_VARS[1:]}
        d2 = {v: (rd(ds, v) if v in ds.variables else nan2) for v in V2D}
        if iy == nyr - 1:
            fm = (np.nanmean(rd(ds, "FOREST_MANAGED"), axis=0)
                  if "FOREST_MANAGED" in ds.variables else nanv.copy())

    pftann[0, iy] = vma[TREE]
    pftann[1, iy] = st["DIA_INV"][TREE]
    pftann[2, iy] = st["HEIGHT"][TREE]
    pftann[3, iy] = st["AGE_STAND_BM"][TREE]
    pftann[4, iy] = st["RDI"][TREE]
    pftann[5, iy] = st["AC12_F_REAL"][TREE]   # janvier seul valide -> nanmean = valeur annuelle
    pftann[6, iy] = st["AGE_STAND_AREA"][TREE]

    w = np.where(np.isfinite(vma), vma, 0.0)
    for k in (1, 2, 3, 4):
        sel = TREE & (pft_cls == k)
        wk = w[sel]                                    # (npft_k,lat,lon)
        ak = wk.sum(axis=0)
        cube[SERIES.index("area_c%d" % k), iy] = ak
        with np.errstate(invalid="ignore", divide="ignore"):
            for pre, var in (("dia_", "DIA_INV"), ("hgt_", "HEIGHT"),
                             ("age_", "AGE_STAND_BM"), ("rdi_", "RDI"),
                             ("ind_", "IND"), ("aga_", "AGE_STAND_AREA")):
                x = st[var][sel]
                num = np.nansum(np.where(np.isfinite(x), x, 0) * wk, axis=0)
                cube[SERIES.index(pre + "c%d" % k), iy] = np.where(ak > 1e-4, num / ak, np.nan)
    wt = w[TREE]
    fa = wt.sum(axis=0)
    cube[SERIES.index("forest_area"), iy] = fa
    with np.errstate(invalid="ignore", divide="ignore"):
        xa = st["AGE_STAND_BM"][TREE]
        num = np.nansum(np.where(np.isfinite(xa), xa, 0) * wt, axis=0)
        cube[SERIES.index("age_porte"), iy] = np.where(fa > 1e-4, num / fa, np.nan)
        xb = st["TOTAL_M_c"][TREE]
        cube[SERIES.index("biomass"), iy] = np.nansum(np.where(np.isfinite(xb), xb, 0) * wt, axis=0)
        xaa = st["AGE_STAND_AREA"][TREE]
        num = np.nansum(np.where(np.isfinite(xaa), xaa, 0) * wt, axis=0)
        cube[SERIES.index("age_area"), iy] = np.where(fa > 1e-4, num / fa, np.nan)

    cube[SERIES.index("class0"), iy] = np.nanmean(d2["CLASS0_VEG"], axis=0)
    cube[SERIES.index("frac0"), iy] = np.nanmean(d2["FRAC_CLASS0"], axis=0)
    cube[SERIES.index("edge"), iy] = np.nanmean(d2["EDGE_LENGTH_DYN"], axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        for nm, v in (("da_fire", "DA_FIRE"), ("da_storm", "DA_STORM"),
                      ("da_pest", "DA_PEST"), ("da_harv", "DA_HARVEST")):
            # accumulateur m2/jour moyenné par mois -> somme annuelle m2 -> fraction de terre
            tot = np.nansum(d2[v] * DAYS[:, None, None], axis=0)
            cube[SERIES.index(nm), iy] = np.where(land_area > 0, tot / land_area, np.nan)

    if iy >= nyr - 5:
        for g in range(11):
            ess_acc[g] += w[TREE & (pft_grp == g)].sum(axis=0)
        ess_n += 1
    if iy % 10 == 0:
        print(f"  {iy+1}/{nyr}", flush=True)

# --- Âge de FLUX pixel (Sigma 1/f) : par essence puis pondération veget_max.
# Conventions : taux annuels AC12_F_REAL lissés en moyenne glissante W=15 ans ZÉROS
# INCLUS (plot_pixel_cards._flux_age) ; T plafonné à T_MAX=200 quand aucun transfert ;
# classe terminale = âge d'ENTRÉE + T4_HALF=12 (map_flux_age, séjour mûr ~24 yr pxneuU).
W_FLUX, T_MAX, T4_HALF = 15, 200.0, 12.0
fnum = np.zeros((nyr, nlat, nlon)); fden = np.zeros((nyr, nlat, nlon))
kern = None
for g, (a, b) in enumerate(TREE_BLOCKS):
    entry = np.zeros((nyr, nlat, nlon))
    for k in range(4):
        pos = tpos[a + k]
        f = np.nan_to_num(pftann[5, :, pos])              # (nyr,lat,lon), zéros inclus
        fbar = np.empty_like(f)
        for t in range(nyr):
            fbar[t] = f[max(0, t - W_FLUX + 1):t + 1].mean(axis=0)
        T = np.where(fbar > 1.0 / T_MAX, 1.0 / np.maximum(fbar, 1e-9), T_MAX)
        T = np.minimum(T, T_MAX)
        age = (entry + T / 2.0) if k < 3 else (entry + T4_HALF)
        w = np.nan_to_num(pftann[0, :, pos])
        fnum += w * age
        fden += w
        if k < 3:
            entry = entry + T
with np.errstate(invalid="ignore", divide="ignore"):
    cube[SERIES.index("age_flux")] = np.where(fden > 1e-4, fnum / fden, np.nan)

# --- Sélection des pixels : terre avec un minimum de forêt à un moment quelconque
fmax = np.nanmax(cube[SERIES.index("forest_area")], axis=0)
landm = (contfrac > 0.05) & np.isfinite(fmax)
pixm = landm & (fmax > 1e-3)
py, px = np.where(pixm)
npix = len(py)
print(f"{npix} pixels retenus (terre+forêt) sur {int(landm.sum())} terrestres", flush=True)

# --- Essences & FM par pixel
ess = ess_acc[:, py, px] / max(ess_n, 1)               # (11,npix) aires moyennes
ess_share = ess / np.maximum(ess.sum(axis=0), 1e-9)
wfm = np.where(np.isfinite(fm), fm, 0)[TREE][:, py, px]  # (npft_tree,npix) valeurs FM
wvm = np.where(np.isfinite(cube[SERIES.index("forest_area"), -1]), 1, 1)
vm_last = None
fm_share = np.zeros((5, npix))
wt_last = np.zeros_like(wfm)
# poids = aire des PFT arbres du dernier fichier : recharger vm dernier fichier
with netCDF4.Dataset(files[-1]) as ds:
    vml = np.nanmean(rd(ds, "VEGET_MAX"), axis=0)
wt_last = np.where(np.isfinite(vml), vml, 0)[TREE][:, py, px]
cls_fm = np.clip(np.rint(wfm), 1, 5).astype(int)
for c in range(1, 6):
    fm_share[c - 1] = np.where(cls_fm == c, wt_last, 0).sum(axis=0)
fm_share = fm_share / np.maximum(fm_share.sum(axis=0), 1e-9)

# --- Quantification int16 par série
data = cube[:, :, py, px]                              # (NS,nyr,npix)
data = np.transpose(data, (0, 2, 1)).copy()            # (NS,npix,nyr)
# Écrêtage physique : les créneaux quasi vides portent des RDI/DA aberrants qui
# écraseraient la résolution int16 ; la fiche les masque de toute façon (AREA_MIN).
CLIP = {"rdi_c1": 3.0, "rdi_c2": 3.0, "rdi_c3": 3.0, "rdi_c4": 3.0,
        "da_fire": 0.5, "da_storm": 0.5, "da_pest": 0.5, "da_harv": 0.5,
        "ind_c1": 1.0, "ind_c2": 1.0, "ind_c3": 1.0, "ind_c4": 1.0}
meta_s = []
q = np.zeros(data.shape, dtype="<i2")
for s in range(NS):
    a = data[s]
    if SERIES[s] in CLIP:
        a = np.minimum(a, CLIP[SERIES[s]])
        data[s] = a
    fin = np.isfinite(a)
    lo = float(np.nanmin(a)) if fin.any() else 0.0
    hi = float(np.nanmax(a)) if fin.any() else 1.0
    if hi <= lo: hi = lo + 1.0
    sc = (hi - lo) / 32000.0
    v = np.where(fin, np.rint((a - lo) / sc), -32768).astype("<i2")
    q[s] = v
    meta_s.append({"name": SERIES[s], "off": lo, "scale": sc})

blob = gzip.compress(q.tobytes(), 9)
b64 = base64.b64encode(blob).decode()

# --- Blob PAR ESSENCE (int8) : essences >= 3 % du couvert forestier, plafond 4/pixel.
# Chaque couple (essence, classe) = un PFT unique -> séries brutes du PFT, sans pondération.
# Layout : [pixel][essence][var 0..4][classe 0..3][année], var = area,dia,hgt,age,rdi.
ESS_STEPS = [1 / 250, 0.8 / 250, 45 / 250, 300 / 250, 2.5 / 250, 0.5 / 250, 300 / 250]
pix_ess = []
for i in range(npix):
    shares = ess_share[:, i]
    ids = [int(g) for g in np.argsort(-shares) if shares[g] >= 0.03][:4]
    pix_ess.append(ids)
segs = []
for i, (j, ii) in enumerate(zip(py, px)):
    for g in pix_ess[i]:
        a, b = TREE_BLOCKS[g]
        for v in range(7):
            stp = ESS_STEPS[v]
            for k in range(4):
                arr = pftann[v, :, tpos[a + k], j, ii]
                qq = np.where(np.isfinite(arr),
                              np.clip(np.rint(arr / stp), 0, 250) - 124,
                              -128).astype("<i1")
                segs.append(qq)
eraw = np.concatenate(segs) if segs else np.zeros(0, dtype="<i1")
eblob = gzip.compress(eraw.tobytes(), 9)
eb64 = base64.b64encode(eblob).decode()
meta = {
    "years": years, "nyr": nyr, "npix": npix, "nseries": NS,
    "series": meta_s, "essences": ESSENCES,
    "lat": [round(float(lat[j]), 3) for j in py],
    "lon": [round(float(lon[i]), 3) for i in px],
    "iy": py.tolist(), "ix": px.tolist(),
    "grid": {"lat": lat.tolist(), "lon": lon.tolist()},
    "ess_share": np.round(ess_share, 4).T.tolist(),    # (npix,11)
    "fm_share": np.round(fm_share, 4).T.tolist(),      # (npix,5)
    "land_area_km2": [round(float(land_area[j, i]) / 1e6, 1) for j, i in zip(py, px)],
    "bare_iy": np.where(landm & ~pixm)[0].tolist(),
    "bare_ix": np.where(landm & ~pixm)[1].tolist(),
    "pix_ess": pix_ess,
    "ess_steps": ESS_STEPS,
}
with open(OUT + "euage_cube_meta.json", "w") as f:
    json.dump(meta, f)
with open(OUT + "euage_cube.b64", "w") as f:
    f.write(b64)
with open(OUT + "euage_ess.b64", "w") as f:
    f.write(eb64)
print(f"OK: cube {q.nbytes/1e6:.1f} Mo brut, {len(blob)/1e6:.2f} Mo gzip, "
      f"{len(b64)/1e6:.2f} Mo b64 ; essences {eraw.nbytes/1e6:.1f} Mo brut, "
      f"{len(eb64)/1e6:.2f} Mo b64", flush=True)
