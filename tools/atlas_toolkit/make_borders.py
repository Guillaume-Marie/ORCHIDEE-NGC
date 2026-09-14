#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Clip a basemap (Natural Earth 50 m coastlines + borders) to the simulation
domain, so that the atlas stays readable outside Europe.

Usage: python3 make_borders.py [atlas_meta.json] [borders_domain.json]

Without cartopy (or without its cached shapefiles), the script stops cleanly:
assemble_atlas.py then falls back to the shipped borders_eu.json.
"""
import json, math, os, sys

meta_path = sys.argv[1] if len(sys.argv) > 1 else "atlas_meta.json"
out_path = sys.argv[2] if len(sys.argv) > 2 else "borders_domain.json"

with open(meta_path) as f:
    m = json.load(f)
g = m["grid"]
lat0, lat1 = min(g["lat"]) - g["dlat"], max(g["lat"]) + g["dlat"]
lon0, lon1 = min(g["lon"]) - g["dlon"], max(g["lon"]) + g["dlon"]
# margin: the coastline must extend a little beyond the frame
mlat, mlon = 0.05 * (lat1 - lat0) + 1, 0.05 * (lon1 - lon0) + 1
lat0, lat1, lon0, lon1 = lat0 - mlat, lat1 + mlat, lon0 - mlon, lon1 + mlon

try:
    from cartopy.io import shapereader
    import shapely.geometry as sgeom
except ImportError:
    sys.exit("cartopy/shapely missing: European basemap kept "
             "(pip install cartopy to clip the basemap to the domain)")

# Below ~60° of extent the 50 m is sharp without being heavy; beyond that the
# 110 m is enough and divides the page weight by five.
res = "50m" if max(lat1 - lat0, lon1 - lon0) < 60 else "110m"
box = sgeom.box(lon0, lat0, lon1, lat1)
tol = max(lat1 - lat0, lon1 - lon0) / 1200.0     # simplification ~1 map px

segs = []


def add(geom):
    if geom.is_empty:
        return
    if geom.geom_type.startswith("Multi") or geom.geom_type == "GeometryCollection":
        for gg in geom.geoms:
            add(gg)
        return
    if geom.geom_type == "Polygon":
        add(geom.exterior)
        for r in geom.interiors:
            add(r)
        return
    xy = list(geom.simplify(tol).coords)
    if len(xy) > 1:
        segs.append([[round(x, 3), round(y, 3)] for x, y in xy])


for cat, name in (("physical", "coastline"), ("cultural", "admin_0_boundary_lines_land")):
    try:
        path = shapereader.natural_earth(resolution=res, category=cat, name=name)
    except Exception as e:
        print(f"  {name} unavailable ({type(e).__name__}) — skipped")
        continue
    for rec in shapereader.Reader(path).geometries():
        try:
            add(rec.intersection(box))
        except Exception:
            pass

if not segs:
    sys.exit("no coastline in the domain: basemap unchanged")
with open(out_path, "w") as f:
    json.dump(segs, f, separators=(",", ":"))
print(f"OK {out_path}: {len(segs)} segments, {sum(len(s) for s in segs)} points, "
      f"{os.path.getsize(out_path) / 1e3:.0f} kB (Natural Earth {res}, "
      f"lat {lat0:.1f}..{lat1:.1f}, lon {lon0:.1f}..{lon1:.1f})")
