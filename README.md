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

- **`tools/atlas_toolkit/`** — standalone interactive atlas generator for ORCHIDEE
  outputs: builds a self-contained `atlas.html` (clickable European map, per-pixel
  panels with age classes, species filter, disturbances, forest edge) from monthly
  `stomate_history` files. No server or account needed — the produced file opens
  offline in any browser. See `tools/atlas_toolkit/README.md` (French).
