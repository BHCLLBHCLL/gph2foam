# gph2foam Technical Documentation

This document describes the conversion pipeline, CLI modes, SPH physics
mapping, transient air-only rewrite, and libAcoustics Farassat 1A (FWH)
setup. Implementation lives under `gph2foam/`; sibling mesh writers are
not reimplemented here.

## 1. Architecture

```
GPH (CRDL-FLD)  →  CGNS/HDF5  →  repair  →  OpenFOAM
   gphdecoding      (temp/keep)   cradlecgns     cgns2foam
                                                  │
                    optional SPH ─────────────────┤
                                                  ├─ --cht-direct   multi-region CHT
                                                  ├─ --transient    air-only pimpleFoam
                                                  └─ --fwh          + Farassat 1A mics
```

| Module | Role |
|--------|------|
| `deps.py` | Locate sibling `gphdecoding` / `cgns2foam` (`GPH2FOAM_*` env or `../`) |
| `convert.py` | Orchestration: GPH→CGNS→repair→foam; optional SPH / FWH steps |
| `repair.py` | HDF5 CGNS cleanup (empty PointList BC delete, BCWall, shapes) |
| `sph.py` | Parse Cradle `.sph` text → `SphPhysics` / regions JSON |
| `transient_heat.py` | `apply_transient_flow`: CHT air mesh → single-region `pimpleFoam` |
| `fwh.py` | `apply_fwh`: libAcoustics Farassat 1A dictionaries |
| `__main__.py` | CLI `python -m gph2foam` |

Package import name: `gph2foam`. Do not duplicate GPH binary parsing or
polyMesh writers; extend the sibling projects or call them through `deps`.

## 2. Core pipeline (all modes)

### 2.1 GPH → CGNS

`gph_to_cgns()` imports `gph2cgns.parse_gph_mesh` / `write_cgns`.

- Default: omit full-mesh `FluidRegion` (`skip_fluid_region=True`) so part
  zones are not duplicated in mono / CHT merges.
- `--keep-fluid-region` restores FLDUTIL parity.
- `--single-zone` emits one zone (implies keeping FluidRegion naming).

### 2.2 CGNS repair (default on)

Ported from `cradlecgns/cgns_repair_summary.py`. Always:

1. `CGNSLibraryVersion` → 4.2  
2. ZoneBC type → `BCWall` (GPH writes `Null`; `open*` names still map to
   `patch` in cgns2foam)  
3. PointList shape `(n,)` → `(n, 1)`  
5. **Delete** ZoneBC groups that have `PointList` but no `" data"`  
6. Zone_t data shape `(1,3)` → `(3,1)` when needed  
7. Empty `ElementConnectivity` datasets  

Optional: `-d` delete `@…` groups; `-i` int32→int64.

Empty PointList BCs must stay deleted: cgns2foam `_read_bc` crashes on
`None.reshape` when `" data"` is missing.

### 2.3 CGNS → OpenFOAM

Delegates to `cgns2foam` `convert_file` / `scan_file`.

| Mode | Flag | Output |
|------|------|--------|
| Mono-block | (default) | One `polyMesh` + scaffold `0/` |
| Coupling scan | `--scan` | Report only |
| Multi-region CHT | `--cht-direct` | `chtMultiRegionSimpleFoam` tree |
| Air-only transient | `--transient` | `pimpleFoam` (after rewrite) |
| Transient + noise | `--fwh` | implies `--transient` + FWH FO |

`--cgns PATH` / `--reuse-cgns` keep or skip re-parsing the intermediate
CGNS. `--regions PATH` is the CHT sidecar; with `--transient`/`--fwh` the
SPH-generated JSON wins.

## 3. SPH parser (`sph.py`)

Input: Cradle scFLOW `.sph` text (UTF-8). Default path: `<gph-stem>.sph`
or `--sph PATH`.

### 3.1 Commands read

| SPH block | Mapped to |
|-----------|-----------|
| `TM_CYCLE` / `TM_CYCLE_PROG` | cycle count / steady flag |
| `EQUA` `TURB` | laminar vs RANS |
| `GRAV` | gravity (used by CHT JSON; **ignored** for `pimpleFoam`) |
| `BASI` | `pbas`, `tmbas` (Celsius → Kelvin offset) |
| `PROPERTY` | materials (fluid/solid, ρ, μ, Cp, κ, buoyancy) |
| `PART_MATERIAL` | part → material → fluid/solid lists |
| `INIT` | default and per-part initial T (°C) |
| `BC_FLOWIO` | opening `TEMP` / `PRES` |
| `SOURCE` | heat sources (W) — CHT only |
| `REF_AXIS` + `ALE_MOVE` | rotation origin, axis, ω (rad/s) |

Command detection: unindented `UPPER_SNAKE` after `#` / `%` / column-0
`/` / SOF. Column-0 `Moving` under `ALE_MOVE` is a body token, not a new
command. Part names may contain digits (`rotation1`, `fin_1`); name vs
numeric rows are separated with `_has_alpha`.

### 3.2 Regions JSON variants

- `to_regions_json()` — full CHT sidecar (`fluid_regions`, `solid_regions`,
  materials, heat sources, optional `buoyancy`, `mrf_regions`).
- `to_flow_regions_json()` — **`--transient` / `--fwh`**: `solid_regions=[]`,
  `g=[0,0,0]`, fluid `mu` only; used so `--cht-direct` meshes air (+ AMI)
  without solids. Fluid–solid faces remain as air boundary patches and are
  later turned into walls.

Helpers: `kinematic_nu()` = μ/ρ; `opening_kinematic_pressure()` = PRES/ρ.

## 4. `--transient` (air-only `pimpleFoam`)

### 4.1 Sequence

1. Parse SPH → temp regions JSON (`flow_only=True`).  
2. Force `cht_direct=True`; build multi-region case (air only).  
3. `apply_transient_flow(out_dir, physics)`:

| Step | Detail |
|------|--------|
| Mesh promote | `constant/air/polyMesh` → `constant/polyMesh` |
| Strip | Delete all other region dirs under `constant/`, `0/`, `system/` |
| Couplings | `mappedWall` → `wall` |
| Transport | `constant/transportProperties` with `nu = mu/rho` |
| Pressure | Kinematic `p`; openings use `totalPressure` with `p0 = PRES/ρ` |
| Velocity | Impeller walls → `movingWallVelocity`; else `noSlip` / inletOutlet |
| Motion | If ALE: drop `MRFProperties`, write `constant/dynamicMeshDict` (`multiSolidBodyMotionSolver` / `rotatingMotion`) |
| Time | Euler + adjustable Δt; with fans `endTime = 10 × T_fastest` |
| Solver | `application pimpleFoam`; PIMPLE + `correctPhi` |
| Turbulence | Optional `k`/`omega`/`nut` when SPH `TURB ≠ 0` |
| Run script | `Allrun` without `-region` (`pimpleFoam` has none) |

### 4.2 What is not included

- No energy equation, no `T` field, no solid conduction.  
- Gravity / Boussinesq not applied (`pimpleFoam` incompressible).  
- Cannot combine with `--cgns-only` or `--scan`.

### 4.3 Case layout (after rewrite)

```
out_transient/
  0/          U, p  [, k, omega, nut]
  constant/   polyMesh/, transportProperties, turbulenceProperties
              [, dynamicMeshDict]
  system/     controlDict, fvSchemes, fvSolution, decomposeParDict
  Allrun
  setup_report.json   # mode: "transient-flow"
```

## 5. `--fwh` (Farassat 1A microphones)

Implies `--transient`, then `apply_fwh(case_dir, physics, summary)`.

### 5.1 Source surfaces

From `constant/polyMesh/boundary`:

- Every `wall` whose name contains `impeller` → its own sampled surface  
  (`PartSurface_impeller1` → surface name `impeller1`).  
  Two fans are **not** merged: libAcoustics orients all faces of one
  surface from a single face-centre average; merging would point normals
  into the gap between rotors.
- All other solid `wall` patches → one surface `laptop_wall`.  
- AMI / openings are **not** FWH sources.  
- At least one impeller wall is required.

### 5.2 Observers

Six mics, 0.5 m outside the AABB of `polyMesh/points` (`mesh_bounds`
supports binary LE float64 and ASCII). Reference point = mean of ALE
origins when present, else box centre. Names: `mic_px/mx/py/my/pz/mz`.

### 5.3 Acoustics settings

| Setting | Value |
|---------|--------|
| Formulation | `FfowcsWilliamsHawkings` / `Farassat1AFormulation` |
| Library | Windows: `libAcoustics.dll`; else `libAcoustics.so` |
| `pName` / `pInf` | `p`, **0** (kinematic pressure) |
| `rhoInf` | SPH fluid density |
| `c0` / `pRef` | 343 m/s, 2e-5 Pa |
| `probeFrequency` | 1 (sample every write / time step) |
| Motion | `nonUniformSurfaceMotion` when sliding zones exist |
| Dual FO | `FwhImpellers` and `FwhImpellersFixedDelay` (delay = one fan period) |

Writes `system/commonSettings`, `system/fwhControl`; injects
`#include "fwhControl"` into `controlDict` `functions`. `Allrun` notes
`acousticData/*.dat` on the master rank. `setup_report.json` mode →
`transient-fwh`.

Use a **separate** output directory if the flow-only case must be kept;
`--fwh` rewrites the same tree.

## 6. CLI reference (modes)

```bash
# Steady CHT (regions JSON)
python -m gph2foam mesh.gph out_cht --cht-direct --regions mesh.json

# Air-only transient flow
python -m gph2foam mesh.gph out_transient --transient
python -m gph2foam mesh.gph out_transient --transient --sph mesh.sph

# Transient + Farassat 1A
python -m gph2foam mesh.gph out_fwh --fwh

# Intermediate CGNS
python -m gph2foam mesh.gph --cgns-only mesh.cgns
python -m gph2foam mesh.gph out --cgns mesh.cgns --reuse-cgns
```

Environment:

| Variable | Meaning |
|----------|---------|
| `GPH2FOAM_GPHDECODING` | Root containing `gph2cgns.py` |
| `GPH2FOAM_CGNS2FOAM` | Root containing `src/convert.py` |

## 7. Tests

```bash
python -m unittest tests.test_repair tests.test_transient_heat -v
```

`tests/test_transient_heat.py` covers SPH fixtures (including real laptop
SPH when present on disk), `apply_transient_flow`, `apply_fwh` (binary
points bounds, impeller requirement, mic layout), and CLI mutual exclusion.

Large meshes under `tests/*.gph` / `cases/` are gitignored; force-add small
JSON / `test_*.py` fixtures.

## 8. Design decisions

1. **Compose, do not fork** mesh I/O — gphdecoding + cgns2foam remain source of truth.  
2. **`--transient` drops solids** — `pimpleFoam` cannot do multi-region CHT; fluid–solid interfaces become walls.  
3. **Sliding mesh instead of MRF** for ALE — needed when the acoustic / near-field unsteady signature of the fan matters; MRF stays available via plain `--cht-direct`.  
4. **FWH surfaces split by impeller** — matches libAcoustics normal orientation behaviour.  
5. **Kinematic pressure** — `pInf=0`, density only through `rhoInf`.  
6. **No `--transient-heat`** — earlier buoyantPimpleFoam / Boussinesq CHT path was replaced by air-only flow; use `--cht-direct` for steady conjugate heat transfer.
