#!/usr/bin/env python3
"""Assemble atlas.html à partir du gabarit et des blobs produits par extract_atlas_cube.py.

Usage : python3 assemble_atlas.py [dossier_des_blobs] [sortie.html]
Défauts : blobs dans le dossier courant, sortie ./atlas.html.
Le fichier produit est AUTONOME (pako embarqué) : il s'ouvre par double-clic,
hors ligne, dans n'importe quel navigateur récent.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
blobs = sys.argv[1] if len(sys.argv) > 1 else "."
out = sys.argv[2] if len(sys.argv) > 2 else "atlas.html"

tpl = open(os.path.join(HERE, "atlas_template.html")).read()
meta = open(os.path.join(blobs, "euage_cube_meta.json")).read()
m = json.loads(meta)
html = (tpl.replace("__META__", meta)
           .replace("__CUBE__", open(os.path.join(blobs, "euage_cube.b64")).read())
           .replace("__ESS__", open(os.path.join(blobs, "euage_ess.b64")).read())
           .replace("__BORDERS__", open(os.path.join(HERE, "borders_eu.json")).read())
           .replace("__PAKO__", open(os.path.join(HERE, "pako.min.js")).read())
           .replace("__Y0__", str(m["years"][0]))
           .replace("__Y1__", str(m["years"][-1])))
for ph in ("__META__", "__CUBE__", "__ESS__", "__BORDERS__", "__PAKO__"):
    assert ph not in html, f"placeholder restant : {ph}"
open(out, "w").write(html)
print(f"OK {out} : {len(html)/1e6:.2f} Mo, années {m['years'][0]}-{m['years'][-1]}, "
      f"{m['npix']} pixels")
