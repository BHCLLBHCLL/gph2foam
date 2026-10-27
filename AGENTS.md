# AGENTS.md

## Project Overview

**gph2foam** orchestrates GPH → OpenFOAM via sibling repos, then optionally
rewrites the case for air-only transient flow and Farassat 1A acoustics:

1. **gphdecoding** — `gph2cgns.parse_gph_mesh` / `write_cgns`
2. **repair** (ported from cradlecgns) — delete empty PointList BCs, etc.
3. **cgns2foam** — `src.convert.convert_file` / `--cht-direct`
4. **sph.py** — Cradle `.sph` → regions JSON (`flow_only` for transient)
5. **transient_heat.apply_transient_flow** — air-only `pimpleFoam`
6. **fwh.apply_fwh** — libAcoustics Farassat 1A (`--fwh`)

Package import name: `gph2foam`. CLI: `python -m gph2foam`.

Technical detail: [docs/TECHNICAL.md](docs/TECHNICAL.md).

## Agent instructions

- Do **not** reimplement GPH binary parsing or OpenFOAM polyMesh writers here;
  extend the sibling projects or call them through `gph2foam.deps`.
- Empty ZoneBC PointLists (no `" data"`) must stay deleted before foam convert —
  cgns2foam's `_read_bc` crashes on `None.reshape`.
- Prefer `--keep-bc-types` only when debugging; default BCWall rewrite matches
  cradlecgns and works with cgns2foam's `open*` name → `patch` override.
- Large meshes live under `tests/` / `cases/` and are gitignored; use
  `git add -f` for small test modules / JSON sidecars.
- Env overrides: `GPH2FOAM_GPHDECODING`, `GPH2FOAM_CGNS2FOAM`.
- `--transient` / `--fwh` need a sibling `.sph` (or `--sph`). They imply
  cht-direct for the fluid mesh, drop solids, and write `pimpleFoam`.
  There is **no** `--transient-heat` / energy equation on this path.
- `--fwh` implies `--transient`. Keep impeller walls as separate FWH
  surfaces; do not merge two fans into one sampledSurface.
- When updating behaviour of SPH / transient / FWH, refresh
  `docs/TECHNICAL.md` and the CLI help / README sections together.

## Quick checks

```bash
python -m gph2foam --help
python -c "from gph2foam.deps import find_gphdecoding_root, find_cgns2foam_root; print(find_gphdecoding_root()); print(find_cgns2foam_root())"
python -m unittest tests.test_repair tests.test_transient_heat -v
python -m gph2foam tests/laptop_thermal_steady_scaled_v3_fanonly.gph --cgns-only /tmp/fanonly.cgns
```
