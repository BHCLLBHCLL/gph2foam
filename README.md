# gph2foam

Convert Software Cradle **GPH** (`CRDL-FLD`) meshes to **OpenFOAM** cases.

```
GPH  →  CGNS/HDF5  →  repair empty/invalid BCs  →  OpenFOAM
       gphdecoding      cradlecgns logic             cgns2foam
                                                  (+ SPH → --transient / --fwh)
```

Full pipeline, SPH mapping, transient rewrite, and Farassat 1A details:
**[docs/TECHNICAL.md](docs/TECHNICAL.md)**.

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

# Air-only transient flow from the GPH plus the sibling .sph (or --sph)
# pimpleFoam, solids dropped, sliding mesh when the SPH has ALE zones
python -m gph2foam mesh.gph out_transient --transient

# Same air-only case, plus libAcoustics Farassat 1A microphone noise
python -m gph2foam mesh.gph out_fwh --fwh

# Coupling scan
python -m gph2foam --scan mesh.gph --report couplings.json
```

### `--transient`

Reads the SPH beside the GPH (`<stem>.sph`, or `--sph`) and writes an air-only
transient flow case. The fluid mesh (including AMI fans) still comes from
`cgns2foam --cht-direct`, with `solid_regions` left empty. The switch then:

- moves `constant/air/polyMesh` to `constant/polyMesh` and deletes every solid
- turns fluid–solid `mappedWall` patches into walls
- sets `application pimpleFoam` with Euler time and PIMPLE (`correctPhi`)
- writes `constant/transportProperties` with `nu = mu/rho` and kinematic `p`
- replaces `MRFProperties` with `constant/dynamicMeshDict` when the SPH has ALE
  zones, and sets impeller walls to `movingWallVelocity`

There is no energy equation, no `T` field, and no solid conduction. `pimpleFoam`
has no `-region` option, so `Allrun` does not pass one. With fans, `endTime`
is ten revolutions of the fastest rotor. `GRAV` is not used: this solver is
incompressible and has no buoyancy term.

`--transient` cannot be combined with `--cgns-only` or `--scan`.
If `--regions` is also given, the SPH wins. The default output directory is
`<gph-directory>/<gph-stem>/`.

### `--fwh`

Implies `--transient`, then writes `system/fwhControl` and
`system/commonSettings` for libAcoustics Farassat 1A. The output directory
is the noise case itself (it is still deleted and rewritten). Pass a new
directory when the flow-only case must be kept.

- Each `wall` whose name contains `impeller` is its own FWH surface
  (`PartSurface_impeller1` → `impeller1`). libAcoustics orients faces
  from one face-centre average per surface, so the two fans are not merged
- Every other solid `wall` is one surface named `laptop_wall`. AMI patches
  and openings are not source surfaces
- Six microphones sit 0.5 m outside the `polyMesh/points` bounding box,
  in the plane of the fan origins when the SPH has ALE zones
- `rhoInf` is the SPH fluid density. `p` is kinematic, so `pInf` is 0
- `nonUniformSurfaceMotion` is on when a sliding zone exists. A second
  function object uses `fixedResponseDelay` of one fan period
- On Windows the function object loads `libAcoustics.dll`

`--fwh` cannot be combined with `--cgns-only` or `--scan`.

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

Unit tests:

```bash
python -m unittest tests.test_repair tests.test_transient_heat -v
```
