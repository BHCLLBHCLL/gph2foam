"""Add libAcoustics Farassat 1A microphone noise to an air-only case.

``--fwh`` runs after :func:`gph2foam.transient_heat.apply_transient_flow`.
Each impeller wall is its own source surface. libAcoustics orients every
face of one sampledSurface from a single face-centre average, so two fans
on one surface point the normals away from the gap between them. Every
other solid wall shares one surface named ``laptop_wall``. AMI patches and
openings are left out. Observers sit 0.5 m outside the mesh bounding box.
``pimpleFoam`` pressure is kinematic, so ``pInf`` stays 0 and ``rhoInf``
is the SPH density.
"""

from __future__ import annotations

import json
import math
import re
import struct
import sys
from pathlib import Path
from typing import Any

from .sph import SphPhysics
from .transient_heat import _fmt, _patches, _read, _write

_MIC_OFFSET = 0.5
_C0 = 343.0
_PREF = 2e-5
_END_MARK = "// ************************************************************************* //"


def library_name() -> str:
    """Shared library token for this platform's OpenFOAM build.

    v2412 on Windows loads ``libAcoustics.dll``. A bare name is prefixed
    with ``lib``, so ``libAcoustics`` would look for ``liblibAcoustics.dll``.
    """
    if sys.platform == "win32":
        return "libAcoustics.dll"
    return "libAcoustics.so"


def mesh_bounds(path: str | Path) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Axis-aligned bounds of an OpenFOAM ``points`` file.

    Binary files are little-endian doubles, which is what cgns2foam writes
    for the v2412 Windows build. ASCII files list ``(x y z)`` triples.
    """
    data = Path(path).read_bytes()
    list_at = -1
    n_points = 0
    search = 0
    while True:
        at = data.find(b"\n(", search)
        if at < 0:
            break
        line_start = data.rfind(b"\n", 0, at)
        count_text = data[line_start + 1:at].strip()
        if count_text.isdigit() and int(count_text) > 0:
            list_at = at
            n_points = int(count_text)
            break
        search = at + 2
    if list_at < 0:
        raise ValueError(f"{path} has no point count")
    binary = bool(re.search(br"format\s+binary\s*;", data[:list_at]))
    payload = data[list_at + 2:]

    mins = [math.inf, math.inf, math.inf]
    maxs = [-math.inf, -math.inf, -math.inf]

    def _accumulate(x: float, y: float, z: float) -> None:
        mins[0] = min(mins[0], x)
        mins[1] = min(mins[1], y)
        mins[2] = min(mins[2], z)
        maxs[0] = max(maxs[0], x)
        maxs[1] = max(maxs[1], y)
        maxs[2] = max(maxs[2], z)

    if binary:
        if payload.startswith(b"\r\n"):
            payload = payload[2:]
        elif payload.startswith(b"\n"):
            payload = payload[1:]
        need = n_points * 24
        if len(payload) < need:
            raise ValueError(
                f"{path} binary point block is {len(payload)} bytes, "
                f"expected at least {need}"
            )
        for offset in range(0, need, 24):
            _accumulate(*struct.unpack_from("<3d", payload, offset))
    else:
        nums = [
            float(token)
            for token in re.findall(
                r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?",
                payload.decode("utf-8", errors="replace"),
            )
        ]
        if len(nums) < n_points * 3:
            raise ValueError(f"{path} ascii point list is short")
        for point in range(n_points):
            _accumulate(*nums[3 * point: 3 * point + 3])

    return (mins[0], mins[1], mins[2]), (maxs[0], maxs[1], maxs[2])


def _impeller_surface_name(patch: str, used: set[str]) -> str:
    """``PartSurface_impeller1`` becomes ``impeller1``."""
    name = patch[len("PartSurface_"):] if patch.startswith("PartSurface_") else patch
    if not name or name in used or name == "laptop_wall":
        name = patch
    used.add(name)
    return name


def _source_surfaces(boundary: str) -> list[tuple[str, list[str]]]:
    """One surface per impeller wall, then ``laptop_wall`` for the rest.

    libAcoustics ``fwhFormulation::initialize`` / ``update`` flip each face
    normal using one face-centre average per sampledSurface. Two impellers
    in that average sit outside both fans, so each fan is its own surface.
    """
    walls = [name for name, typ in _patches(boundary) if typ == "wall"]
    impellers = [name for name in walls if "impeller" in name.lower()]
    others = [name for name in walls if "impeller" not in name.lower()]
    if not impellers:
        raise ValueError(
            "--fwh needs at least one wall patch whose name contains "
            "'impeller' (the rotating FWH surface)"
        )
    used: set[str] = set()
    surfaces = [(_impeller_surface_name(name, used), [name]) for name in impellers]
    if others:
        surfaces.append(("laptop_wall", others))
    return surfaces


def _origins_and_period(rotations: list[dict[str, Any]]) -> tuple[list[list[float]], float]:
    origins = [list(entry["origin"]) for entry in rotations if "origin" in entry]
    omegas = [abs(float(entry["omega"])) for entry in rotations if abs(float(entry.get("omega", 0.0))) > 0]
    period = (2.0 * math.pi) / max(omegas) if omegas else 0.01
    return origins, period


def _observers(
    bounds_min: tuple[float, float, float],
    bounds_max: tuple[float, float, float],
    reference: tuple[float, float, float],
) -> list[tuple[str, tuple[float, float, float]]]:
    rx, ry, rz = reference
    xmin, ymin, zmin = bounds_min
    xmax, ymax, zmax = bounds_max
    return [
        ("mic_px", (xmax + _MIC_OFFSET, ry, rz)),
        ("mic_mx", (xmin - _MIC_OFFSET, ry, rz)),
        ("mic_py", (rx, ymax + _MIC_OFFSET, rz)),
        ("mic_my", (rx, ymin - _MIC_OFFSET, rz)),
        ("mic_pz", (rx, ry, zmax + _MIC_OFFSET)),
        ("mic_mz", (rx, ry, zmin - _MIC_OFFSET)),
    ]


def _common_settings(
    *,
    rho_inf: float,
    cofr: tuple[float, float, float],
    end_time: float,
    observers: list[tuple[str, tuple[float, float, float]]],
) -> str:
    obs = []
    for name, (x, y, z) in observers:
        obs.append(
            f"""            {name}
            {{
                position    ({_fmt(x)} {_fmt(y)} {_fmt(z)});
                pRef        {_fmt(_PREF)};
                fftFreq     1024;
            }}"""
        )
    cx, cy, cz = cofr
    return f"""        // Included inside each FfowcsWilliamsHawkings object.
        libs            ("{library_name()}");

        log             true;

        // v2412 calls functionObject::write() every time step, and
        // libAcoustics samples inside write(). 1 keeps every step.
        probeFrequency  1;

        timeStart       0;
        timeEnd         {_fmt(end_time)};

        writeFft        false;

        c0              {_fmt(_C0)};

        pName           p;

        // pimpleFoam pressure is kinematic. The library multiplies by
        // rhoInf when p is not in pascal, then subtracts pInf.
        pInf            0;

        rhoName         rho;
        rho             rhoInf;
        rhoInf          {_fmt(rho_inf)};

        // Moment centre only. Sliding-mesh face velocity comes from
        // nonUniformSurfaceMotion, not from one rigid rotation here.
        CofR            ({_fmt(cx)} {_fmt(cy)} {_fmt(cz)});

        observers
        {{
{chr(10).join(obs)}
        }}
"""


def _fwh_control(
    surfaces: list[tuple[str, list[str]]],
    *,
    moving: bool,
    response_delay: float,
) -> str:
    motion = "true" if moving else "false"
    surface_body = "\n\n".join(
        f"""        {name}
        {{
            type            patch;
            patches         ({" ".join(patches)});
            interpolate     false;
        }}"""
        for name, patches in surfaces
    )
    blocks = []
    for name, fixed in (
        ("FwhImpellers", "false"),
        ("FwhImpellersFixedDelay", "true"),
    ):
        blocks.append(
            f"""{name}
{{
    type            FfowcsWilliamsHawkings;
    formulationType Farassat1AFormulation;
    #include        "commonSettings"
    patches         ();
    interpolationScheme cellPoint;

    surfaces
    (
{surface_body}
    );

    cleanFreq       100;
    nonUniformSurfaceMotion {motion};
    Ufwh            (0 0 0);
    U0              (0 0 0);
    fixedResponseDelay {fixed};
    responseDelay   {_fmt(response_delay)};
}}"""
        )
    return "\n\n".join(blocks) + "\n"


def _insert_functions(control: str) -> str:
    include = '#include "fwhControl"'
    if include in control:
        return control
    block = "functions\n{\n    " + include + "\n}\n\n"
    idx = control.rfind(_END_MARK)
    if idx < 0:
        return control.rstrip() + "\n\n" + block
    return control[:idx] + block + control[idx:]


def _note_allrun(text: str) -> str:
    if "acousticData" in text:
        return text
    note = (
        "# libAcoustics writes acousticData/FwhImpellers-time.dat and\n"
        "# acousticData/FwhImpellersFixedDelay-time.dat on the master rank.\n"
    )
    marker = "#------------------------------------------------------------------------------\n"
    first = text.find(marker)
    if first < 0:
        return note + text
    insert_at = first + len(marker)
    return text[:insert_at] + note + text[insert_at:]


def apply_fwh(
    case_dir: str | Path,
    physics: SphPhysics,
    summary: dict[str, Any],
) -> dict[str, Any]:
    """Write Farassat 1A dictionaries into an air-only ``pimpleFoam`` case."""
    case_dir = Path(case_dir)
    boundary_path = case_dir / "constant" / "polyMesh" / "boundary"
    points_path = case_dir / "constant" / "polyMesh" / "points"
    control_path = case_dir / "system" / "controlDict"
    if not boundary_path.is_file() or not points_path.is_file():
        raise FileNotFoundError(
            f"{case_dir} has no constant/polyMesh; --fwh runs on the "
            "air-only mesh from --transient"
        )
    if not control_path.is_file():
        raise FileNotFoundError(f"{control_path} is missing")

    fluid = physics.fluid_material()
    if fluid is None or not fluid.rho:
        raise ValueError("--fwh needs the SPH fluid density (rhoInf)")

    surfaces = _source_surfaces(_read(boundary_path))
    bounds_min, bounds_max = mesh_bounds(points_path)
    rotations = list(summary.get("rotations") or [])
    if not rotations:
        rotations = [
            {"origin": list(rot.origin), "omega": rot.omega}
            for rot in physics.rotations
        ]
    origins, period = _origins_and_period(rotations)
    if origins:
        reference = (
            sum(origin[0] for origin in origins) / len(origins),
            sum(origin[1] for origin in origins) / len(origins),
            sum(origin[2] for origin in origins) / len(origins),
        )
    else:
        reference = tuple((lo + hi) / 2.0 for lo, hi in zip(bounds_min, bounds_max))
    observers = _observers(bounds_min, bounds_max, reference)
    end_time = float(summary.get("endTime", 0.1))
    moving = bool(rotations)

    _write(
        case_dir / "system" / "commonSettings",
        _common_settings(
            rho_inf=float(fluid.rho),
            cofr=reference,
            end_time=end_time,
            observers=observers,
        ),
    )
    _write(
        case_dir / "system" / "fwhControl",
        _fwh_control(surfaces, moving=moving, response_delay=period),
    )
    _write(control_path, _insert_functions(_read(control_path)))

    allrun = case_dir / "Allrun"
    if allrun.is_file():
        _write(allrun, _note_allrun(_read(allrun)))

    report_path = case_dir / "setup_report.json"
    if report_path.is_file():
        report = json.loads(report_path.read_text(encoding="utf-8"))
        report["mode"] = "transient-fwh"
        report["fwh"] = {
            "formulation": "Farassat1A",
            "surfaces": [
                {"name": name, "patches": patches} for name, patches in surfaces
            ],
            "rhoInf": float(fluid.rho),
            "pInf": 0.0,
            "responseDelay": period,
            "observers": [
                {"name": name, "position": list(pos)} for name, pos in observers
            ],
        }
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    return {
        "rhoInf": float(fluid.rho),
        "surfaces": [
            {"name": name, "patches": patches} for name, patches in surfaces
        ],
        "observers": [
            {"name": name, "position": list(pos)} for name, pos in observers
        ],
        "responseDelay": period,
        "lib": library_name(),
        "CofR": list(reference),
        "nonUniformSurfaceMotion": moving,
    }
