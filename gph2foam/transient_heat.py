"""Turn a cht-direct case into a single-region ``pimpleFoam`` case.

``cgns2foam --cht-direct`` is still used to build the merged air mesh,
including AMI fans. This module then:

* drops every solid region and the heat-transfer fields
* moves ``constant/air/polyMesh`` to ``constant/polyMesh``
* turns fluid–solid ``mappedWall`` patches into ordinary walls
* replaces ``MRFProperties`` with ``constant/dynamicMeshDict``
* writes incompressible ``p`` / ``U`` and ``nu = mu/rho``

``pimpleFoam`` in OpenFOAM v2412 has no ``-region`` option and does not
solve an energy equation.
"""

from __future__ import annotations

import json
import math
import re
import shutil
from pathlib import Path
from typing import Any

from .sph import SphPhysics

_FOAM_BANNER = """\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2412                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
"""

_HEADER_END = (
    b"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //"
)

_PATCH_RE = re.compile(
    r"^(?P<indent>[ \t]+)(?P<name>\S+)\n(?P=indent)\{\n(?P<body>.*?)\n(?P=indent)\}",
    re.M | re.S,
)


def _header(obj: str, location: str, cls: str = "dictionary") -> str:
    return (
        _FOAM_BANNER
        + "FoamFile\n{\n"
        + "    version     2.0;\n"
        + "    format      ascii;\n"
        + f"    class       {cls};\n"
        + f'    location    "{location}";\n'
        + f"    object      {obj};\n"
        + "}\n"
        + "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
    )


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8").replace("\r\n", "\n")


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def _fmt(value: float) -> str:
    return f"{value:.16g}"


def _parse_mrf(text: str) -> list[dict[str, Any]]:
    blocks = re.split(r"\n(?=\S)", text)
    out: list[dict[str, Any]] = []
    for block in blocks:
        cz = re.search(r"cellZone\s+(\S+);", block)
        origin = re.search(r"origin\s+\(([^)]+)\);", block)
        axis = re.search(r"axis\s+\(([^)]+)\);", block)
        omega = re.search(r"omega\s+([^;]+);", block)
        if not (cz and origin and axis and omega):
            continue
        out.append(
            {
                "cellZone": cz.group(1),
                "origin": [float(x) for x in origin.group(1).split()],
                "axis": [float(x) for x in axis.group(1).split()],
                "omega": float(omega.group(1)),
            }
        )
    return out


def _dynamic_mesh_dict(entries: list[dict[str, Any]]) -> str:
    zones: list[str] = []
    for entry in entries:
        ox, oy, oz = entry["origin"]
        ax, ay, az = entry["axis"]
        zones.append(
            f"""    {entry['cellZone']}
    {{
        solidBodyMotionFunction rotatingMotion;
        rotatingMotionCoeffs
        {{
            origin      ({_fmt(ox)} {_fmt(oy)} {_fmt(oz)});
            axis        ({_fmt(ax)} {_fmt(ay)} {_fmt(az)});
            omega       {_fmt(entry['omega'])};
        }}
    }}"""
        )
    body = "\n\n".join(zones)
    return (
        _header("dynamicMeshDict", "constant")
        + f"""
dynamicFvMesh   dynamicMotionSolverFvMesh;

solver          multiSolidBodyMotionSolver;
motionSolver    multiSolidBodyMotionSolver;

multiSolidBodyMotionSolverCoeffs
{{
{body}
}}

// ************************************************************************* //
"""
    )


def _end_time(entries: list[dict[str, Any]]) -> tuple[float, float, float, float]:
    """Return ``(endTime, deltaT, writeInterval, maxDeltaT)``.

    With fans, ten revolutions of the fastest rotor. ``omega`` is rad/s.
    """
    omegas = [abs(float(e["omega"])) for e in entries if abs(float(e["omega"])) > 0]
    if not omegas:
        return (1.0, 1e-3, 0.1, 0.05)
    period = (2.0 * math.pi) / max(omegas)
    return (10.0 * period, 1e-5, max(2.0 * period, 1e-3), 2e-4)


def _control_dict(end_time: float, delta_t: float, write_interval: float, max_delta_t: float) -> str:
    return (
        _header("controlDict", "system")
        + f"""
application     pimpleFoam;

startFrom       startTime;
startTime       0;
stopAt          endTime;

// Ten fan revolutions when ALE rotation is present, else 1 s.
endTime         {_fmt(end_time)};

deltaT          {_fmt(delta_t)};

writeControl    adjustable;
writeInterval   {_fmt(write_interval)};
purgeWrite      0;

writeFormat     binary;
writePrecision  8;
writeCompression off;

timeFormat      general;
timePrecision   8;
runTimeModifiable true;

adjustTimeStep  yes;
maxCo           0.6;
maxDeltaT       {_fmt(max_delta_t)};

// ************************************************************************* //
"""
    )


def _allrun(n_procs: int, *, moving: bool) -> str:
    reconstruct = ""
    if moving:
        reconstruct = (
            "runApplication -o -s reconstructParMesh reconstructParMesh -constant\n"
        )
    return f"""#!/bin/sh
set -e
cd "${{0%/*}}" || exit
. ${{WM_PROJECT_DIR:?}}/bin/tools/RunFunctions
#------------------------------------------------------------------------------
# Single-region transient flow. Solids and energy are not in this case.
# pimpleFoam reads constant/polyMesh and, when fans are present,
# constant/dynamicMeshDict. It has no -region option.
runApplication -o -s decomposePar decomposePar -force
runParallel -o -np {n_procs} $(getApplication)
{reconstruct}runApplication -o -s reconstructPar reconstructPar
#------------------------------------------------------------------------------
"""


def _schemes(*, turbulent: bool) -> str:
    turb = ""
    if turbulent:
        turb = (
            "    div(phi,k)      Gauss upwind;\n"
            "    div(phi,omega)  Gauss upwind;\n"
        )
    return (
        _header("fvSchemes", "system")
        + f"""
ddtSchemes {{ default Euler; }}
gradSchemes {{ default Gauss linear; }}

divSchemes
{{
    div(phi,U)      Gauss upwind;
{turb}    div((nuEff*dev2(T(grad(U))))) Gauss linear;
    div(meshPhi,p)  Gauss linear;
}}
laplacianSchemes {{ default Gauss linear limited 0.333; }}
interpolationSchemes {{ default linear; }}
snGradSchemes {{ default limited 0.333; }}
wallDist {{ method meshWave; }}

// ************************************************************************* //
"""
    )


def _solution(*, turbulent: bool) -> str:
    if turbulent:
        velocity = """\
    "(U|k|omega)"
    {
        solver           PBiCGStab;
        preconditioner   DILU;
        tolerance        1e-7;
        relTol           0.1;
    }

    "(U|k|omega)Final"
    {
        solver           PBiCGStab;
        preconditioner   DILU;
        tolerance        1e-7;
        relTol           0;
    }"""
    else:
        velocity = """\
    U
    {
        solver           PBiCGStab;
        preconditioner   DILU;
        tolerance        1e-7;
        relTol           0.1;
    }

    UFinal
    {
        $U;
        relTol           0;
    }"""
    return (
        _header("fvSolution", "system")
        + f"""
solvers
{{
    p
    {{
        solver           GAMG;
        smoother         GaussSeidel;
        tolerance        1e-7;
        relTol           0.01;
    }}

    pFinal
    {{
        $p;
        relTol           0;
    }}

    "pcorr.*"
    {{
        solver           GAMG;
        smoother         GaussSeidel;
        tolerance        1e-2;
        relTol           0;
    }}

{velocity}
}}

PIMPLE
{{
    momentumPredictor yes;
    correctPhi        yes;
    nOuterCorrectors  2;
    nCorrectors       2;
    nNonOrthogonalCorrectors 1;
}}

relaxationFactors
{{
    equations
    {{
        "U.*"           1;
    }}
}}

// ************************************************************************* //
"""
    )


def _transport(nu: float) -> str:
    return (
        _header("transportProperties", "constant")
        + f"""
transportModel  Newtonian;

nu              [0 2 -1 0 0 0 0] {_fmt(nu)};

// ************************************************************************* //
"""
    )


def _turbulence(turbulent: bool) -> str:
    if not turbulent:
        body = "simulationType  laminar;\n"
    else:
        body = """\
simulationType  RAS;

RAS
{
    RASModel        kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}
"""
    return _header("turbulenceProperties", "constant") + "\n" + body + "\n// ************************************************************************* //\n"


def _is_opening(name: str, patch_type: str) -> bool:
    return patch_type == "patch" or name == "open" or name.startswith("open")


def _field_U(patches: list[tuple[str, str]]) -> str:
    blocks: list[str] = []
    for name, patch_type in patches:
        if patch_type == "cyclicAMI":
            body = "        type            cyclicAMI;"
        elif _is_opening(name, patch_type):
            body = (
                "        type            pressureInletOutletVelocity;\n"
                "        value           uniform (0 0 0);"
            )
        elif "impeller" in name.lower():
            body = (
                "        type            movingWallVelocity;\n"
                "        value           uniform (0 0 0);"
            )
        else:
            body = "        type            noSlip;"
        blocks.append(f"    {name}\n    {{\n{body}\n    }}")
    joined = "\n\n".join(blocks)
    return (
        _header("U", "0", "volVectorField")
        + f"""
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);

boundaryField
{{
    #includeEtc "caseDicts/setConstraintTypes"

{joined}
}}

// ************************************************************************* //
"""
    )


def _field_p(patches: list[tuple[str, str]], p_open: float) -> str:
    blocks: list[str] = []
    for name, patch_type in patches:
        if patch_type == "cyclicAMI":
            body = "        type            cyclicAMI;"
        elif _is_opening(name, patch_type):
            body = (
                "        type            totalPressure;\n"
                f"        p0              uniform {_fmt(p_open)};"
            )
        else:
            body = "        type            zeroGradient;"
        blocks.append(f"    {name}\n    {{\n{body}\n    }}")
    joined = "\n\n".join(blocks)
    return (
        _header("p", "0", "volScalarField")
        + f"""
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;

boundaryField
{{
    #includeEtc "caseDicts/setConstraintTypes"

{joined}
}}

// ************************************************************************* //
"""
    )


def _scalar_field(
    name: str,
    dimensions: str,
    internal: float,
    patches: list[tuple[str, str]],
    *,
    wall: str,
    opening: str,
) -> str:
    blocks: list[str] = []
    for patch_name, patch_type in patches:
        if patch_type == "cyclicAMI":
            body = "        type            cyclicAMI;"
        elif _is_opening(patch_name, patch_type):
            body = opening
        else:
            body = wall
        blocks.append(f"    {patch_name}\n    {{\n{body}\n    }}")
    joined = "\n\n".join(blocks)
    return (
        _header(name, "0", "volScalarField")
        + f"""
dimensions      {dimensions};
internalField   uniform {_fmt(internal)};

boundaryField
{{
    #includeEtc "caseDicts/setConstraintTypes"

{joined}
}}

// ************************************************************************* //
"""
    )


def _patches(text: str) -> list[tuple[str, str]]:
    found: list[tuple[str, str]] = []
    for match in _PATCH_RE.finditer(text):
        typ = re.search(r"type\s+(\S+);", match.group("body"))
        if typ is None:
            continue
        found.append((match.group("name"), typ.group(1)))
    if not found:
        raise ValueError("polyMesh/boundary has no patches")
    return found


def _walls_from_mapped(text: str) -> str:
    """Fluid–solid ``mappedWall`` faces become ordinary walls."""

    def _replace(match: re.Match[str]) -> str:
        body = match.group("body")
        typ = re.search(r"type\s+(\S+);", body)
        if typ is None or typ.group(1) != "mappedWall":
            return match.group(0)
        n_faces = re.search(r"nFaces\s+(\d+)\s*;", body)
        start = re.search(r"startFace\s+(\d+)\s*;", body)
        if n_faces is None or start is None:
            raise ValueError(
                f"mappedWall patch {match.group('name')} is missing nFaces/startFace"
            )
        indent = match.group("indent")
        inner = indent + ("\t" if "\t" in indent else "    ")
        name = match.group("name")
        return (
            f"{indent}{name}\n"
            f"{indent}{{\n"
            f"{inner}type wall;\n"
            f"{inner}inGroups 1(wall);\n"
            f"{inner}nFaces {n_faces.group(1)};\n"
            f"{inner}startFace {start.group(1)};\n"
            f"{indent}}}"
        )

    updated = _PATCH_RE.sub(_replace, text)
    if "mappedWall" in updated or "sampleRegion" in updated:
        raise ValueError("failed to remove mappedWall/sampleRegion from the air boundary")
    return updated


def _retarget_mesh_header(path: Path) -> None:
    data = path.read_bytes()
    idx = data.find(_HEADER_END)
    if idx < 0:
        idx = min(len(data), 4000)
    head, tail = data[:idx], data[idx:]
    replaced = head.replace(b"constant/air/polyMesh", b"constant/polyMesh")
    if replaced != head:
        path.write_bytes(replaced + tail)


def _n_procs(case_dir: Path) -> int:
    path = case_dir / "system" / "decomposeParDict"
    if not path.is_file():
        return 8
    match = re.search(r"numberOfSubdomains\s+(\d+)", _read(path))
    return int(match.group(1)) if match else 8


def _strip_regions(case_dir: Path) -> None:
    constant = case_dir / "constant"
    for child in list(constant.iterdir()):
        if child.name == "polyMesh":
            continue
        if child.is_dir():
            shutil.rmtree(child)
        else:
            child.unlink()
    zero = case_dir / "0"
    if zero.exists():
        shutil.rmtree(zero)
    system = case_dir / "system"
    if system.is_dir():
        for child in list(system.iterdir()):
            if child.is_dir():
                shutil.rmtree(child)


def apply_transient_flow(case_dir: str | Path, physics: SphPhysics) -> dict[str, Any]:
    """Rewrite ``case_dir`` into an air-only ``pimpleFoam`` case."""
    case_dir = Path(case_dir)
    air = case_dir / "constant" / "air"
    mesh_src = air / "polyMesh"
    if not mesh_src.is_dir():
        raise FileNotFoundError(
            f"{case_dir} has no constant/air/polyMesh; --transient needs "
            "the fluid mesh from --cht-direct"
        )

    nu = physics.kinematic_nu()
    p_open = physics.opening_kinematic_pressure()
    rotations: list[dict[str, Any]] = []
    mrf_path = air / "MRFProperties"
    if mrf_path.is_file():
        rotations = _parse_mrf(_read(mrf_path))
    if physics.rotations and not rotations:
        raise FileNotFoundError(
            f"{case_dir} has no usable constant/air/MRFProperties, but the SPH "
            "defines ALE rotation "
            + ", ".join(rot.part for rot in physics.rotations)
        )
    n_procs = _n_procs(case_dir)

    boundary_path = mesh_src / "boundary"
    boundary = _walls_from_mapped(_read(boundary_path))
    _write(boundary_path, boundary)
    for entry in mesh_src.iterdir():
        if entry.is_file():
            _retarget_mesh_header(entry)

    dest = case_dir / "constant" / "polyMesh"
    if dest.exists():
        shutil.rmtree(dest)
    shutil.move(str(mesh_src), str(dest))
    _strip_regions(case_dir)

    if rotations:
        _write(case_dir / "constant" / "dynamicMeshDict", _dynamic_mesh_dict(rotations))
    _write(case_dir / "constant" / "transportProperties", _transport(nu))
    _write(
        case_dir / "constant" / "turbulenceProperties",
        _turbulence(physics.turbulent),
    )

    patches = _patches(_read(dest / "boundary"))
    _write(case_dir / "0" / "U", _field_U(patches))
    _write(case_dir / "0" / "p", _field_p(patches, p_open))
    if physics.turbulent:
        _write(
            case_dir / "0" / "k",
            _scalar_field(
                "k",
                "[0 2 -2 0 0 0 0]",
                0.1,
                patches,
                wall="        type            kqRWallFunction;\n        value           uniform 0.1;",
                opening=(
                    "        type            inletOutlet;\n"
                    "        inletValue      uniform 0.1;\n"
                    "        value           uniform 0.1;"
                ),
            ),
        )
        _write(
            case_dir / "0" / "omega",
            _scalar_field(
                "omega",
                "[0 0 -1 0 0 0 0]",
                10.0,
                patches,
                wall="        type            omegaWallFunction;\n        value           uniform 10;",
                opening=(
                    "        type            inletOutlet;\n"
                    "        inletValue      uniform 10;\n"
                    "        value           uniform 10;"
                ),
            ),
        )
        _write(
            case_dir / "0" / "nut",
            _scalar_field(
                "nut",
                "[0 2 -1 0 0 0 0]",
                0.0,
                patches,
                wall="        type            nutkWallFunction;\n        value           uniform 0;",
                opening="        type            calculated;\n        value           uniform 0;",
            ),
        )

    end_time, delta_t, write_interval, max_delta_t = _end_time(rotations)
    _write(
        case_dir / "system" / "controlDict",
        _control_dict(end_time, delta_t, write_interval, max_delta_t),
    )
    _write(case_dir / "system" / "fvSchemes", _schemes(turbulent=physics.turbulent))
    _write(case_dir / "system" / "fvSolution", _solution(turbulent=physics.turbulent))
    _write(case_dir / "Allrun", _allrun(n_procs, moving=bool(rotations)))

    report_path = case_dir / "setup_report.json"
    if report_path.is_file():
        report = json.loads(report_path.read_text(encoding="utf-8"))
        report["mode"] = "transient-flow"
        report["solver"] = "pimpleFoam"
        report["dynamicMesh"] = rotations
        report.pop("mrf", None)
        report.pop("heat_sources", None)
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    return {
        "solver": "pimpleFoam",
        "endTime": end_time,
        "deltaT": delta_t,
        "nu": nu,
        "rotations": rotations,
    }
