# ORCHIDEE — NextGenCarbon development version (source code)

This repository provides the **source code** of the ORCHIDEE land surface model
version developed and used within the **NextGenCarbon** project, covering forest
age-class dynamics, edge effects, disturbance regimes (windthrow, bark beetle,
fire) and forest management developments.

- Base: ORCHIDEE trunk, SVN r9551 (`svn://forge.ipsl.fr/orchidee/trunk/ORCHIDEE`)
- Plus the age-class / edge-effect / maturity-transfer / stand-age developments
  carried out for NextGenCarbon.

## Scope

Source code only: Fortran sources (`src_*`), XIOS XML definitions (`src_xml/`),
the FCM build system (`makeorchidee_fcm`, `config.fcm`, `bld.cfg`, `arch/`) and
the canonical parameter reference (`orchidee.default`). Running the model
additionally requires the IPSL offline environment (modipsl, libIGCM, IOIPSL,
XIOS) and forcing datasets, which are not part of this repository.

## License

This software is governed by the **CeCILL license** (GPL-compatible French free
software license) — see `ORCHIDEE_CeCILL.LIC`.

## Official distribution

This repository is a project snapshot, **not** the official ORCHIDEE
distribution channel. For the reference model, documentation and released
versions, see https://orchidee.ipsl.fr and the IPSL forge
(https://forge.ipsl.fr/orchidee).

## Tools

The post-processing tools live in their own repositories:

- **orchidee_atlas** — https://github.com/Guillaume-Marie/orchidee_atlas — standalone
  interactive atlas generator: builds a self-contained `atlas.html` (clickable map,
  per-cell time-series plates with age classes, species filter, disturbances) from
  ORCHIDEE history files, driven by a single configuration file. The produced file
  opens offline in any browser.
- **orchidee_extract** — https://github.com/Guillaume-Marie/orchidee_extract —
  recipe-driven extraction and aggregation of ORCHIDEE outputs to zarr / netCDF /
  parquet, with an explicit semantics layer (area basis of each field, units).
