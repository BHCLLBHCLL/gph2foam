# AGENTS.md

## Project Overview

**gph2foam** orchestrates GPH → OpenFOAM via sibling repos:

1. **gphdecoding** — `gph2cgns.parse_gph_mesh` / `write_cgns`
2. **repair** (ported from cradlecgns) — delete empty PointList BCs, etc.
3. **cgns2foam** — `src.convert.convert_file` / `--cht-direct`

Package import name: `gph2foam`. CLI: `python -m gph2foam`.

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

## Quick checks

```bash
python -m gph2foam --help
python -c "from gph2foam.deps import find_gphdecoding_root, find_cgns2foam_root; print(find_gphdecoding_root()); print(find_cgns2foam_root())"
python -m gph2foam tests/laptop_thermal_steady_scaled_v3_fanonly.gph --cgns-only /tmp/fanonly.cgns
```
