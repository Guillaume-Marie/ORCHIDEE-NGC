#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Assemble atlas.html: template + blobs produced by extract_atlas_cube.py.

Usage: python3 assemble_atlas.py [blob_directory] [output.html]
Defaults: blobs in the current directory, output ./atlas.html.
The produced file is STANDALONE (pako embedded): it opens by double-click,
offline, in any recent browser.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
blobs = sys.argv[1] if len(sys.argv) > 1 else "."
out = sys.argv[2] if len(sys.argv) > 2 else "atlas.html"


def rd(*parts):
    with open(os.path.join(*parts)) as f:
        return f.read()


# basemap: the one clipped to the domain if it exists, otherwise the shipped one
borders = os.path.join(blobs, "borders_domain.json")
if not os.path.exists(borders):
    borders = os.path.join(HERE, "borders_eu.json")
    print(f"note: {borders} used (run make_borders.py to clip the basemap "
          "to the domain)")

meta = rd(blobs, "atlas_meta.json")
m = json.loads(meta)
html = (rd(HERE, "atlas_template.html")
        .replace("__META__", meta)
        .replace("__CUBE__", rd(blobs, "atlas_cube.b64"))
        .replace("__ESS__", rd(blobs, "atlas_ess.b64"))
        .replace("__BORDERS__", rd(borders))
        .replace("__PAKO__", rd(HERE, "pako.min.js")))
for ph in ("__META__", "__CUBE__", "__ESS__", "__BORDERS__", "__PAKO__"):
    assert ph not in html, f"remaining placeholder: {ph}"
with open(out, "w") as f:
    f.write(html)
print(f"OK {out}: {len(html) / 1e6:.2f} MB · {m['run']} · {m['years'][0]}-{m['years'][-1]} · "
      f"{m['npix']} pixels · {len(m['essences'])} species x {m['ncls']} classes")
