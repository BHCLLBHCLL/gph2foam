"""Read a Cradle scFLOW ``.sph`` text file into a regions JSON.

``cgns2foam --cht-direct`` does not read the SPH.  This parser covers the
blocks the transient heat case needs: materials, gravity, reference state,
part assignment, initial temperature, opening temperature, heat sources,
and ALE rotation axes.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

# OpenFOAM perfectGas: rho = p * W / (R_specific * T), W in g/mol,
# R_universal = 8.314462618 J/(mol·K).  molWeight is W.
_R_UNIVERSAL = 8.314462618

_COMMAND_RE = re.compile(r"^[A-Z][A-Z0-9_]*$")
_FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")
_MATERIAL_HEAD_RE = re.compile(r"^\s*(\d+)\s+general\s+(.+?)\s*$")


@dataclass
class Material:
    name: str
    kind: str  # "fluid" | "solid"
    rho: float | None = None
    mu: float | None = None
    cp: float | None = None
    kappa: float | None = None
    buoyancy: bool = False
    beta: float | None = None
    tref_c: float | None = None


@dataclass
class Rotation:
    part: str
    axis_name: str
    origin: tuple[float, float, float]
    axis: tuple[float, float, float]
    omega: float


@dataclass
class SphPhysics:
    """Physics taken from one SPH file."""

    path: Path
    cycle: int | None = None
    steady: bool = True
    turbulent: bool = False
    gravity: tuple[float, float, float] = (0.0, 0.0, 0.0)
    pbas: float = 101325.0
    tmbas: float = 273.15
    materials: dict[str, Material] = field(default_factory=dict)
    part_material: dict[str, str] = field(default_factory=dict)
    default_temp_c: float = 25.0
    part_temp_c: dict[str, float] = field(default_factory=dict)
    opening_temp_c: float | None = None
    opening_pres: float | None = None
    heat_sources_w: dict[str, float] = field(default_factory=dict)
    rotations: list[Rotation] = field(default_factory=list)

    def kelvin(self, celsius: float) -> float:
        return float(celsius) + self.tmbas

    @property
    def fluid_parts(self) -> list[str]:
        out: list[str] = []
        for part, mat_name in self.part_material.items():
            mat = self.materials.get(mat_name)
            if mat is not None and mat.kind == "fluid":
                out.append(part)
        return out

    @property
    def solid_parts(self) -> list[str]:
        out: list[str] = []
        for part, mat_name in self.part_material.items():
            mat = self.materials.get(mat_name)
            if mat is not None and mat.kind == "solid":
                out.append(part)
        return out

    def fluid_temperature_k(self) -> float:
        """One temperature for the merged ``air`` region.

        ``FluidRegion`` wins when the SPH sets it.  Otherwise the first
        fluid part that has its own initial temperature is used.
        """
        if "FluidRegion" in self.part_temp_c:
            return self.kelvin(self.part_temp_c["FluidRegion"])
        for part in self.fluid_parts:
            if part in self.part_temp_c:
                return self.kelvin(self.part_temp_c[part])
        return self.kelvin(self.default_temp_c)

    def solid_temperature_k(self) -> float:
        return self.kelvin(self.default_temp_c)

    def opening_temperature_k(self) -> float:
        if self.opening_temp_c is None:
            return self.solid_temperature_k()
        return self.kelvin(self.opening_temp_c)

    def fluid_material(self) -> Material | None:
        for part in self.fluid_parts:
            mat = self.materials.get(self.part_material[part])
            if mat is not None and mat.kind == "fluid":
                return mat
        for mat in self.materials.values():
            if mat.kind == "fluid":
                return mat
        return None

    def to_regions_json(self) -> dict[str, Any]:
        """Sidecar consumed by ``cgns2foam --cht-direct``.

        Part names are matched onto CGNS zone names by that loader.
        Extra keys (``buoyancy``, ``transient_heat``) are ignored there
        and kept so the SPH revision is visible next to the CGNS.
        """
        fluid = self.fluid_material()
        materials: dict[str, dict[str, float]] = {}
        buoyancy: dict[str, float] | None = None
        if fluid is not None:
            entry: dict[str, float] = {}
            if fluid.mu is not None:
                entry["mu"] = fluid.mu
            if fluid.cp is not None:
                entry["Cp"] = fluid.cp
            if fluid.mu is not None and fluid.cp is not None and fluid.kappa:
                entry["Pr"] = fluid.mu * fluid.cp / fluid.kappa
            rho_ref = fluid.rho
            t_ref = self.fluid_temperature_k()
            if fluid.buoyancy and fluid.tref_c is not None:
                t_ref = self.kelvin(fluid.tref_c)
            if rho_ref is not None and self.pbas > 0 and t_ref > 0:
                entry["molWeight"] = (
                    rho_ref * _R_UNIVERSAL * t_ref / self.pbas * 1000.0
                )
            if entry:
                materials["air"] = entry
            if fluid.buoyancy and fluid.rho is not None and fluid.beta is not None:
                buoyancy = {
                    "rho0": fluid.rho,
                    "beta": fluid.beta,
                    "T0": self.kelvin(fluid.tref_c if fluid.tref_c is not None else self.default_temp_c),
                }
        for part in self.solid_parts:
            mat = self.materials[self.part_material[part]]
            solid: dict[str, float] = {}
            if mat.rho is not None:
                solid["rho"] = mat.rho
            if mat.cp is not None:
                solid["Cp"] = mat.cp
            if mat.kappa is not None:
                solid["kappa"] = mat.kappa
            if solid:
                materials[part] = solid

        data: dict[str, Any] = {
            "transient_heat": True,
            "sph": str(self.path),
            "fluid_regions": list(self.fluid_parts),
            "solid_regions": list(self.solid_parts),
            "heat_sources": dict(self.heat_sources_w),
            "materials": materials,
            "g": list(self.gravity),
            "turbulence_model": "RANS" if self.turbulent else "laminar",
            "initial_conditions": {
                "T": self.solid_temperature_k(),
                "p": self.pbas,
            },
            "n_procs": 8,
            "endTime": 100,
            "writeInterval": 100,
        }
        if self.rotations:
            data["mrf_regions"] = [
                {
                    "cellZone": rot.part,
                    "origin": list(rot.origin),
                    "axis": list(rot.axis),
                    "omega": rot.omega,
                }
                for rot in self.rotations
            ]
        if buoyancy is not None:
            data["buoyancy"] = buoyancy
        return data

    def kinematic_nu(self) -> float:
        """``nu = mu/rho`` for incompressible ``pimpleFoam``."""
        fluid = self.fluid_material()
        if fluid is None or fluid.mu is None or not fluid.rho:
            raise ValueError(
                "pimpleFoam needs a fluid density and viscosity (nu = mu/rho)"
            )
        return fluid.mu / fluid.rho

    def opening_kinematic_pressure(self) -> float:
        """Gauge kinematic pressure ``PRES/rho`` on ``total_pressure`` openings."""
        fluid = self.fluid_material()
        if fluid is None or not fluid.rho:
            raise ValueError("pimpleFoam needs a fluid density")
        pres = 0.0 if self.opening_pres is None else self.opening_pres
        return pres / fluid.rho

    def to_flow_regions_json(self) -> dict[str, Any]:
        """Regions JSON that meshes the fluid parts only.

        Solids are omitted so ``--cht-direct`` does not write them. The
        fluid side of each fluid–solid face stays in the air boundary and
        is later turned into a wall for ``pimpleFoam``.
        """
        fluid = self.fluid_material()
        materials: dict[str, dict[str, float]] = {}
        if fluid is not None and fluid.mu is not None:
            materials["air"] = {"mu": fluid.mu}
        data: dict[str, Any] = {
            "transient_flow": True,
            "sph": str(self.path),
            "fluid_regions": list(self.fluid_parts),
            "solid_regions": [],
            "materials": materials,
            "g": [0.0, 0.0, 0.0],
            "turbulence_model": "RANS" if self.turbulent else "laminar",
            "initial_conditions": {"T": 300.0, "p": 0.0},
            "n_procs": 8,
            "endTime": 100,
            "writeInterval": 100,
        }
        if self.rotations:
            data["mrf_regions"] = [
                {
                    "cellZone": rot.part,
                    "origin": list(rot.origin),
                    "axis": list(rot.axis),
                    "omega": rot.omega,
                }
                for rot in self.rotations
            ]
        return data


def _commands(text: str) -> dict[str, list[str]]:
    """Split an SPH into ``{COMMAND: body lines}``.

    A command is an unindented ``UPPER_SNAKE`` token that follows ``#``,
    ``%``, a column-0 ``/``, or the start of the file.  Body lines such as
    ``Moving`` sit at column 0 too, but they are not commands.  A column-0
    ``/`` stays in the body so INIT / ALE blocks can close, and it also
    allows the next real command to start.
    """
    blocks: dict[str, list[str]] = {}
    current: str | None = None
    armed = True
    for raw in text.splitlines():
        line = raw.rstrip()
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("%") or stripped.startswith("#"):
            armed = True
            continue
        if stripped == "/" and line[:1] not in (" ", "\t"):
            armed = True
            if current is not None:
                blocks[current].append(line)
            continue
        if (
            armed
            and _COMMAND_RE.match(stripped)
            and line[:1] not in (" ", "\t")
        ):
            current = stripped
            blocks.setdefault(current, [])
            armed = False
            continue
        armed = False
        if current is not None:
            blocks[current].append(line)
    return blocks


def _floats(line: str) -> list[float]:
    return [float(tok) for tok in _FLOAT_RE.findall(line)]


def _has_alpha(line: str) -> bool:
    """True when *line* is a name, not a numeric row.

    Part names such as ``rotation1`` and ``fin_1`` contain digits, so a
    float scan alone cannot tell them from ``25 0``.
    """
    return any(ch.isalpha() for ch in line)


def _parse_property(lines: list[str]) -> dict[str, Material]:
    materials: dict[str, Material] = {}
    current: Material | None = None
    expect: str | None = None
    for line in lines:
        head = _MATERIAL_HEAD_RE.match(line)
        if head:
            current = Material(name=head.group(2).strip(), kind="solid")
            materials[current.name] = current
            expect = None
            continue
        if current is None:
            continue
        token = line.strip()
        if token in ("fluid", "solid"):
            current.kind = token
            continue
        if token in ("density", "viscosity", "specific_heat", "heat_conductivity"):
            expect = token
            continue
        if token in ("incompressible_buoyancy", "incompressible"):
            if token == "incompressible_buoyancy":
                current.buoyancy = True
            continue
        if token in ("value", "/"):
            continue
        nums = _floats(line)
        if expect == "density" and nums:
            current.rho = nums[0]
            if current.buoyancy and len(nums) >= 3:
                current.beta = nums[1]
                current.tref_c = nums[2]
            expect = None
        elif expect == "viscosity" and nums:
            current.mu = nums[0]
            expect = None
        elif expect == "specific_heat" and nums:
            current.cp = nums[0]
            expect = None
        elif expect == "heat_conductivity" and nums:
            current.kappa = nums[0]
            expect = None
    return materials


def _parse_part_material(lines: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped == "/":
            continue
        part, _, mat = stripped.partition(" ")
        if mat:
            out[part] = mat.strip()
    return out


def _parse_init(lines: list[str]) -> tuple[float, dict[str, float]]:
    """Return ``(default_celsius, {part: celsius})``.

    A ``TEMP`` block with no part names is the default.  ``/`` closes a
    block; it does not end the INIT command (several blocks share it).
    """
    default = 25.0
    per_part: dict[str, float] = {}
    value: float | None = None
    parts: list[str] = []
    want_value = False

    def _flush() -> None:
        nonlocal default, value, parts
        if value is None:
            parts = []
            return
        if parts:
            for part in parts:
                per_part[part] = value
        else:
            default = value
        value = None
        parts = []

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("%"):
            continue
        if stripped == "/":
            _flush()
            want_value = False
            continue
        if stripped == "TEMP":
            _flush()
            want_value = True
            continue
        if want_value and not _has_alpha(stripped):
            nums = _floats(stripped)
            if nums:
                value = nums[0]
                want_value = False
            continue
        if value is not None and _has_alpha(stripped):
            parts.append(stripped)
    _flush()
    return default, per_part


def _parse_bc_flowio(lines: list[str]) -> tuple[float | None, float | None]:
    temp: float | None = None
    pres: float | None = None
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("TEMP"):
            nums = _floats(stripped)
            if nums:
                temp = nums[0]
        elif stripped.startswith("PRES"):
            nums = _floats(stripped)
            if nums:
                pres = nums[0]
    return temp, pres


def _parse_source(lines: list[str]) -> dict[str, float]:
    out: dict[str, float] = {}
    pending: float | None = None
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped == "/":
            pending = None
            continue
        if stripped.startswith("TEMP"):
            nums = _floats(stripped)
            pending = nums[0] if nums else None
            continue
        if stripped.lower().startswith("source"):
            continue
        if pending is not None and _has_alpha(stripped):
            out[stripped] = pending
    return out


def _parse_axes(lines: list[str]) -> dict[str, tuple[tuple[float, float, float], tuple[float, float, float]]]:
    axes: dict[str, tuple[tuple[float, float, float], tuple[float, float, float]]] = {}
    name: str | None = None
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped == "/":
            name = None
            continue
        nums = _floats(stripped)
        if name and len(nums) >= 6 and not _has_alpha(stripped):
            origin = (nums[0], nums[1], nums[2])
            axis = (nums[3], nums[4], nums[5])
            axes[name] = (origin, axis)
            name = None
            continue
        if _has_alpha(stripped):
            name = stripped
    return axes


def _parse_ale(
    lines: list[str],
    axes: dict[str, tuple[tuple[float, float, float], tuple[float, float, float]]],
) -> list[Rotation]:
    """One ``Moving`` block: omega, axis name, part name."""
    rotations: list[Rotation] = []
    chunk: list[str] = []

    def _flush(body: list[str]) -> None:
        if not body:
            return
        omega: float | None = None
        names: list[str] = []
        for line in body:
            stripped = line.strip()
            if not stripped or stripped == "/":
                continue
            if _has_alpha(stripped):
                names.append(stripped)
                continue
            nums = _floats(stripped)
            # ``omega  0`` — the first value is rad/s.  Leading integers
            # (1, 2) are Cradle motion-mode flags and are ignored.
            if omega is None and len(nums) >= 2:
                omega = nums[0]
        if omega is None or len(names) < 2:
            return
        axis_name, part = names[-2], names[-1]
        origin_axis = axes.get(axis_name)
        if origin_axis is None:
            raise ValueError(
                f"ALE part {part!r} references {axis_name!r}, "
                "which is not in REF_AXIS"
            )
        origin, axis = origin_axis
        rotations.append(
            Rotation(
                part=part,
                axis_name=axis_name,
                origin=origin,
                axis=axis,
                omega=float(omega),
            )
        )

    for line in lines:
        token = line.strip()
        # ``MovingAxis_*`` is a name inside the block, not a new block.
        if token == "Moving" or token.startswith("Moving["):
            _flush(chunk)
            chunk = []
            continue
        chunk.append(line)
    _flush(chunk)
    return rotations


def parse_sph(path: str | Path) -> SphPhysics:
    """Parse ``path`` and return the physics the transient case uses."""
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="replace")
    blocks = _commands(text)
    physics = SphPhysics(path=path.resolve())

    if "TM_CYCLE" in blocks:
        nums = _floats(" ".join(blocks["TM_CYCLE"]))
        if nums:
            physics.cycle = int(nums[0])
    if "TM_CYCLE_PROG" in blocks:
        physics.steady = "steady" in " ".join(blocks["TM_CYCLE_PROG"]).lower()
    if "EQUA" in blocks:
        for line in blocks["EQUA"]:
            if line.strip().startswith("TURB"):
                nums = _floats(line)
                physics.turbulent = bool(nums and nums[0] != 0)
    if "GRAV" in blocks:
        nums = _floats(" ".join(blocks["GRAV"]))
        if len(nums) >= 3:
            physics.gravity = (nums[0], nums[1], nums[2])
    if "BASI" in blocks:
        nums = _floats(" ".join(blocks["BASI"]))
        if nums:
            physics.pbas = nums[0]
        if len(nums) >= 2:
            physics.tmbas = nums[1]
    if "PROPERTY" in blocks:
        physics.materials = _parse_property(blocks["PROPERTY"])
    if "PART_MATERIAL" in blocks:
        physics.part_material = _parse_part_material(blocks["PART_MATERIAL"])
    if "INIT" in blocks:
        physics.default_temp_c, physics.part_temp_c = _parse_init(blocks["INIT"])
    if "BC_FLOWIO" in blocks:
        physics.opening_temp_c, physics.opening_pres = _parse_bc_flowio(blocks["BC_FLOWIO"])
    if "SOURCE" in blocks:
        physics.heat_sources_w = _parse_source(blocks["SOURCE"])
    axes = _parse_axes(blocks.get("REF_AXIS", []))
    if "ALE_MOVE" in blocks:
        physics.rotations = _parse_ale(blocks["ALE_MOVE"], axes)

    missing = [
        name for name in physics.part_material.values() if name not in physics.materials
    ]
    if missing:
        raise ValueError(
            "SPH PART_MATERIAL names missing from PROPERTY: " + ", ".join(missing)
        )
    if not physics.fluid_parts:
        raise ValueError(f"SPH {path} has no fluid part in PART_MATERIAL")
    return physics


def write_regions_json(
    physics: SphPhysics,
    dest: str | Path,
    *,
    flow_only: bool = False,
) -> Path:
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    data = physics.to_flow_regions_json() if flow_only else physics.to_regions_json()
    dest.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return dest
