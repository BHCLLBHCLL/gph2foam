# gph2foam

Convert Software Cradle **GPH** (`CRDL-FLD`) meshes to **OpenFOAM** cases.

```
GPH  →  CGNS/HDF5  →  repair empty/invalid BCs  →  OpenFOAM
       gphdecoding      cradlecgns logic             cgns2foam
```

## Dependencies

Sibling toolchains (discovered automatically under `../`, or via env):

| Env / path | Role |
|------------|------|
| `GPH2FOAM_GPHDECODING` or `../gphdecoding` | GPH parse + CGNS write (`gph2cgns.py`) |
| `GPH2FOAM_CGNS2FOAM` or `../cgns2foam` | CGNS → OpenFOAM (`python -m src`) |

Python packages:

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Mono-block OpenFOAM case (prefer --single-zone to avoid duplicating FluidRegion+parts)
python -m gph2foam tests/laptop_thermal_steady_scaled_v3_fanonly.gph cases/fanonly --single-zone

# Keep intermediate CGNS; later reuse it without re-parsing GPH
python -m gph2foam mesh.gph out_case --cgns mesh.cgns
python -m gph2foam mesh.gph out_case --cgns mesh.cgns --reuse-cgns

# CGNS only (repaired)
python -m gph2foam mesh.gph --cgns-only mesh_fixed.cgns

# Multi-region CHT (needs regions JSON matching zone names from GPH)
python -m gph2foam mesh.gph out_cht --cht-direct --regions mesh.json

# Coupling scan
python -m gph2foam --scan mesh.gph --report couplings.json
```

### Why repair?

`gph2cgns` mirrors vendor FLDUTIL: surface regions missing from a zone are still
emitted as ZoneBC groups whose `PointList` has **no** `" data"` dataset.
Those empty BCs break the OpenFOAM converter and are deleted by default
(same as cradlecgns step 5). ZoneBC type `Null` is rewritten to `BCWall`
(openings still become `patch` when the name starts with `open`).

By default the full-mesh `FluidRegion` zone is **omitted** (`--keep-fluid-region`
to restore FLDUTIL parity). Prefer `--cht-direct` for multi-part thermal cases.

## Test mesh

Place large samples under `tests/` (gitignored). Example:

`tests/laptop_thermal_steady_scaled_v3_fanonly.gph`

CHT sidecar (tracked): `tests/laptop_thermal_steady_scaled_v3_fanonly.json`
