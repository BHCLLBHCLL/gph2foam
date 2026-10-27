"""SPH physics and the transient-flow dictionary rewrite.

The fixture follows the laptop thermal SPH: column-0 ``/`` inside INIT,
column-0 ``Moving`` under ALE_MOVE, and part names that contain digits.
"""

from __future__ import annotations

import contextlib
import io
import json
import struct
import tempfile
import unittest
from pathlib import Path

from gph2foam.__main__ import main
from gph2foam.fwh import apply_fwh, library_name, mesh_bounds
from gph2foam.sph import parse_sph
from gph2foam.transient_heat import apply_transient_flow

_SPH = """\
#
TM_CYCLE
    1000
#
TM_CYCLE_PROG
    steady
#
EQUA
 TURB   0
/
#
GRAV
    0 -9.81 0
#
BASI
    101325 273.15 0
#
PROPERTY
    1 general air(incompressible/20C)
      fluid
        density
          incompressible_buoyancy
            1.206 0.003411 25
        viscosity
          value
            1.83e-05
        specific_heat
          value
            1007
        heat_conductivity
          value
            0.0256
      /
    2 general quartz_glass(300K)
      solid
        density
          incompressible
            2190
        specific_heat
          value
            740
        heat_conductivity
          value
            1.38
      /
    3 general copper(Cu)(300K)
      solid
        density
          incompressible
            8880
        specific_heat
          value
            386
        heat_conductivity
          value
            398
      /
/
#
PART_MATERIAL
  air_domain air(incompressible/20C)
  CPU quartz_glass(300K)
  fin_1 copper(Cu)(300K)
  rotation1 air(incompressible/20C)
  rotation2 air(incompressible/20C)
/
#
INIT
% @InitialDefaultTemperature
  TEMP
    25 0
% Initial
  TEMP
    20 -1
      FluidRegion
/
% Initial[2]
  TEMP
    20 -1
      rotation1
      rotation2
/
/
#
BC_FLOWIO
  Flux[4]
    total_pressure
    PRES 0
    TEMP 25
/
    open
/
/
#
SOURCE
  Source[2]
    TEMP   15   0   11
         CPU
/
/
#
ALE_MOVE
Moving
  1
  2
  628.318530717959 0
  MovingAxis_1
  rotation1
/
Moving[2]
  1
  2
  628.318530717959 0
  MovingAxis_2
  rotation2
/
/
#
REF_AXIS
  MovingAxis_1
    -0.067799999999992 -0.004198957848 0.081000000000001 0 1 0
  MovingAxis_2
    0.0801999999999901 -0.001798957847803 0.0810000000000001 0 -1 0
/
#
"""

_SPH_STILL = """\
#
EQUA
 TURB   0
/
#
GRAV
    0 0 0
#
BASI
    101325 273.15 0
#
PROPERTY
    1 general air(incompressible/20C)
      fluid
        density
          incompressible
            1.2
        viscosity
          value
            1.8e-05
        specific_heat
          value
            1000
        heat_conductivity
          value
            0.026
      /
    2 general quartz_glass(300K)
      solid
        density
          incompressible
            2190
        specific_heat
          value
            740
        heat_conductivity
          value
            1.38
      /
/
#
PART_MATERIAL
  air_domain air(incompressible/20C)
  CPU quartz_glass(300K)
/
#
INIT
  TEMP
    25 0
/
#
"""

_BOUNDARY = """\
FoamFile
{
    version     2.0;
    format      ascii;
    class       polyBoundaryMesh;
    location    "constant/air/polyMesh";
    object      boundary;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

4
(
	ami_rotation1_air_domain
	{
		type cyclicAMI;
		neighbourPatch ami_air_domain_rotation1;
		startFace 0;
		nFaces 10;
	}
	air_to_CPU
	{
		type mappedWall;
		sampleMode nearestPatchFace;
		sampleRegion CPU;
		samplePatch CPU_to_air;
		startFace 10;
		nFaces 4;
	}
	open
	{
		type patch;
		startFace 14;
		nFaces 2;
	}
	PartSurface_impeller1
	{
		type wall;
		startFace 16;
		nFaces 3;
	}
)
"""

_MRF = """\
MRF1
{
    cellZone            FPHPARTS_rotation1;
    active              yes;
    nonRotatingPatches ();
    origin              (-0.067799999999992 -0.004198957848 0.081000000000001);
    axis                (0 1 0);
    omega               628.318530717959;
}
"""


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


class ParseSphTest(unittest.TestCase):
    def test_laptop_shaped_sph(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "case.sph"
            path.write_text(_SPH, encoding="utf-8")
            physics = parse_sph(path)

        self.assertEqual(physics.gravity, (0.0, -9.81, 0.0))
        self.assertEqual(physics.pbas, 101325.0)
        self.assertEqual(physics.tmbas, 273.15)
        self.assertFalse(physics.turbulent)
        self.assertEqual(physics.default_temp_c, 25.0)
        self.assertEqual(physics.part_temp_c["FluidRegion"], 20.0)
        self.assertEqual(physics.part_temp_c["rotation1"], 20.0)
        self.assertEqual(physics.part_temp_c["rotation2"], 20.0)
        self.assertAlmostEqual(physics.fluid_temperature_k(), 293.15)
        self.assertAlmostEqual(physics.solid_temperature_k(), 298.15)
        self.assertAlmostEqual(physics.opening_temperature_k(), 298.15)
        self.assertEqual(physics.heat_sources_w, {"CPU": 15.0})
        self.assertEqual(physics.fluid_parts, ["air_domain", "rotation1", "rotation2"])
        self.assertEqual(physics.solid_parts, ["CPU", "fin_1"])
        self.assertEqual(len(physics.rotations), 2)
        self.assertEqual(physics.rotations[0].part, "rotation1")
        self.assertEqual(physics.rotations[0].axis, (0.0, 1.0, 0.0))
        self.assertAlmostEqual(physics.rotations[0].omega, 628.318530717959)
        self.assertEqual(physics.rotations[1].axis, (0.0, -1.0, 0.0))
        self.assertAlmostEqual(physics.rotations[0].origin[0], -0.067799999999992)

        data = physics.to_regions_json()
        self.assertEqual(data["g"], [0.0, -9.81, 0.0])
        self.assertEqual(data["heat_sources"], {"CPU": 15.0})
        self.assertEqual(data["initial_conditions"]["T"], 298.15)
        self.assertEqual(data["mrf_regions"][0]["cellZone"], "rotation1")
        self.assertEqual(data["materials"]["CPU"]["kappa"], 1.38)
        self.assertEqual(data["materials"]["fin_1"]["rho"], 8880)
        air = data["materials"]["air"]
        self.assertAlmostEqual(air["Pr"], 1.83e-5 * 1007 / 0.0256)
        self.assertAlmostEqual(
            air["molWeight"],
            1.206 * 8.314462618 * 298.15 / 101325.0 * 1000.0,
        )
        self.assertAlmostEqual(data["buoyancy"]["rho0"], 1.206)
        self.assertAlmostEqual(data["buoyancy"]["beta"], 0.003411)
        self.assertAlmostEqual(data["buoyancy"]["T0"], 298.15)
        flow = physics.to_flow_regions_json()
        self.assertEqual(flow["solid_regions"], [])
        self.assertNotIn("heat_sources", flow)
        self.assertEqual(flow["g"], [0.0, 0.0, 0.0])
        self.assertIn("rotation1", flow["fluid_regions"])
        self.assertNotIn("CPU", flow["fluid_regions"])
        self.assertAlmostEqual(physics.kinematic_nu(), 1.83e-5 / 1.206)
        self.assertAlmostEqual(physics.opening_kinematic_pressure(), 0.0)

    def test_real_laptop_sph(self) -> None:
        path = Path(
            r"D:\training\cradle\laptop\laptop_thermal_steady_scaled_v3_transient"
            r"\laptop_thermal_steady_scaled_v3_BCs.sph"
        )
        if not path.is_file():
            self.skipTest("laptop SPH is not on this machine")
        physics = parse_sph(path)
        self.assertEqual(physics.gravity, (0.0, -9.81, 0.0))
        self.assertEqual(physics.heat_sources_w["CPU"], 15.0)
        self.assertEqual(len(physics.rotations), 2)
        self.assertTrue(all(abs(r.omega - 628.318530717959) < 1e-9 for r in physics.rotations))
        self.assertAlmostEqual(physics.fluid_temperature_k(), 293.15)
        self.assertAlmostEqual(physics.opening_temperature_k(), 298.15)
        data = physics.to_regions_json()
        self.assertAlmostEqual(data["buoyancy"]["rho0"], 1.206)
        self.assertAlmostEqual(data["buoyancy"]["beta"], 0.003411)
        self.assertIn("CPU", data["solid_regions"])
        self.assertIn("rotation1", data["fluid_regions"])
        self.assertIn("Cu_block", data["materials"])
        flow = physics.to_flow_regions_json()
        self.assertEqual(flow["solid_regions"], [])
        self.assertNotIn("CPU", flow["fluid_regions"])
        self.assertAlmostEqual(physics.kinematic_nu(), 1.83e-5 / 1.206)


class ApplyTransientFlowTest(unittest.TestCase):
    def _case(self, root: Path, *, with_mrf: bool) -> None:
        _write(root / "constant" / "air" / "polyMesh" / "boundary", _BOUNDARY)
        _write(root / "constant" / "air" / "polyMesh" / "points", "points\n")
        _write(root / "constant" / "g", "value (0 -9.81 0);\n")
        _write(root / "constant" / "regionProperties", "regions ( fluid (air) solid (CPU) );\n")
        _write(root / "constant" / "CPU" / "polyMesh" / "points", "solid\n")
        _write(root / "0" / "air" / "T", "internalField uniform 298.15;\n")
        _write(root / "0" / "CPU" / "T", "internalField uniform 298.15;\n")
        _write(root / "system" / "CPU" / "fvSchemes", "ddtSchemes { default steadyState; }\n")
        _write(root / "system" / "air" / "fvSolution", "SIMPLE { }\n")
        _write(root / "system" / "decomposeParDict", "numberOfSubdomains 4;\n")
        _write(
            root / "setup_report.json",
            json.dumps(
                {
                    "mode": "cht-direct",
                    "solver": "chtMultiRegionSimpleFoam",
                    "mrf": [1],
                    "heat_sources": {"CPU": 15},
                }
            ),
        )
        if with_mrf:
            _write(root / "constant" / "air" / "MRFProperties", _MRF)

    def test_air_only_pimplefoam(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sph = Path(tmp) / "case.sph"
            sph.write_text(_SPH, encoding="utf-8")
            physics = parse_sph(sph)
            case = Path(tmp) / "openfoam"
            self._case(case, with_mrf=True)
            summary = apply_transient_flow(case, physics)

            self.assertFalse((case / "constant" / "air").exists())
            self.assertFalse((case / "constant" / "CPU").exists())
            self.assertFalse((case / "constant" / "g").exists())
            self.assertFalse((case / "constant" / "regionProperties").exists())
            self.assertFalse((case / "0" / "air").exists())
            self.assertFalse((case / "0" / "T").exists())
            self.assertFalse((case / "system" / "CPU").exists())
            self.assertFalse((case / "system" / "air").exists())

            boundary = (case / "constant" / "polyMesh" / "boundary").read_text(encoding="utf-8")
            self.assertIn('location    "constant/polyMesh"', boundary)
            self.assertIn("type wall;", boundary)
            self.assertNotIn("mappedWall", boundary)
            self.assertNotIn("sampleRegion", boundary)
            self.assertIn("cyclicAMI", boundary)

            mesh = (case / "constant" / "dynamicMeshDict").read_text(encoding="utf-8")
            self.assertIn('location    "constant"', mesh)
            self.assertIn("multiSolidBodyMotionSolver", mesh)
            self.assertIn("FPHPARTS_rotation1", mesh)
            self.assertFalse((case / "constant" / "air" / "dynamicMeshDict").exists())

            control = (case / "system" / "controlDict").read_text(encoding="utf-8")
            self.assertIn("application     pimpleFoam;", control)
            self.assertNotIn("buoyantPimpleFoam", control)
            self.assertAlmostEqual(summary["endTime"], 0.1, places=6)
            self.assertAlmostEqual(summary["nu"], 1.83e-5 / 1.206)

            velocity = (case / "0" / "U").read_text(encoding="utf-8")
            self.assertIn("movingWallVelocity", velocity)
            self.assertIn("pressureInletOutletVelocity", velocity)
            self.assertIn("noSlip", velocity)
            pressure = (case / "0" / "p").read_text(encoding="utf-8")
            self.assertIn("dimensions      [0 2 -2 0 0 0 0];", pressure)
            self.assertIn("internalField   uniform 0;", pressure)
            self.assertIn("totalPressure", pressure)
            self.assertNotIn("p_rgh", pressure)

            transport = (case / "constant" / "transportProperties").read_text(encoding="utf-8")
            self.assertIn("Newtonian", transport)
            self.assertIn(f"{summary['nu']:.16g}", transport)
            schemes = (case / "system" / "fvSchemes").read_text(encoding="utf-8")
            self.assertIn("Euler", schemes)
            self.assertNotIn("div(phi,h)", schemes)
            solution = (case / "system" / "fvSolution").read_text(encoding="utf-8")
            self.assertIn("correctPhi", solution)
            self.assertNotIn("p_rgh", solution)
            self.assertNotIn("SIMPLE", solution)

            allrun = (case / "Allrun").read_text(encoding="utf-8")
            self.assertIn("decomposePar -force", allrun)
            commands = [line for line in allrun.splitlines() if not line.startswith("#")]
            self.assertFalse(any("-region" in line for line in commands))
            self.assertNotIn("-copyZero", allrun)
            self.assertIn("reconstructParMesh -constant", allrun)
            self.assertIn("-np 4 $(getApplication)", allrun)

            report = json.loads((case / "setup_report.json").read_text(encoding="utf-8"))
            self.assertEqual(report["solver"], "pimpleFoam")
            self.assertEqual(report["mode"], "transient-flow")
            self.assertNotIn("mrf", report)
            self.assertNotIn("heat_sources", report)

    def test_no_rotation_stays_static(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sph = Path(tmp) / "still.sph"
            sph.write_text(_SPH_STILL, encoding="utf-8")
            physics = parse_sph(sph)
            case = Path(tmp) / "openfoam"
            self._case(case, with_mrf=False)
            summary = apply_transient_flow(case, physics)
            self.assertFalse((case / "constant" / "dynamicMeshDict").exists())
            self.assertEqual(summary["endTime"], 1.0)
            self.assertAlmostEqual(summary["nu"], 1.8e-5 / 1.2)
            allrun = (case / "Allrun").read_text(encoding="utf-8")
            self.assertNotIn("reconstructParMesh", allrun)
            commands = [line for line in allrun.splitlines() if not line.startswith("#")]
            self.assertFalse(any("-region" in line for line in commands))


class CliTest(unittest.TestCase):
    def test_help_lists_the_switch(self) -> None:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            with self.assertRaises(SystemExit) as ctx:
                main(["--help"])
        self.assertEqual(ctx.exception.code, 0)
        text = buf.getvalue()
        self.assertIn("--transient", text)
        self.assertIn("--fwh", text)
        self.assertNotIn("--transient-heat", text)
        self.assertIn("pimpleFoam", text)
        self.assertIn("Farassat", text)
        self.assertIn("--sph", text)

    def test_scan_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            gph = Path(tmp) / "mesh.gph"
            gph.write_text("not a mesh", encoding="utf-8")
            code = main(["--transient", "--scan", str(gph)])
            fwh_code = main(["--fwh", "--scan", str(gph)])
        self.assertEqual(code, 2)
        self.assertEqual(fwh_code, 2)


_POINTS = """\
FoamFile
{
    version     2.0;
    format      ascii;
    class       vectorField;
    location    "constant/air/polyMesh";
    object      points;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

2
(
(-0.2 -0.1 -0.3)
(0.3 0.2 0.4)
)

"""


class ApplyFwhTest(unittest.TestCase):
    def test_writes_farassat_mics(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sph = root / "case.sph"
            sph.write_text(_SPH, encoding="utf-8")
            physics = parse_sph(sph)
            case = root / "openfoam"
            ApplyTransientFlowTest()._case(case, with_mrf=True)
            boundary = case / "constant" / "air" / "polyMesh" / "boundary"
            boundary.write_text(
                boundary.read_text(encoding="utf-8").replace(
                    "\tnFaces 3;\n\t}\n)\n",
                    "\tnFaces 3;\n\t}\n"
                    "\tPartSurface_impeller2\n\t{\n"
                    "\t\ttype wall;\n"
                    "\t\tstartFace 19;\n"
                    "\t\tnFaces 3;\n"
                    "\t}\n)\n",
                ),
                encoding="utf-8",
            )
            _write(case / "constant" / "air" / "polyMesh" / "points", _POINTS)
            summary = apply_transient_flow(case, physics)
            acoustics = apply_fwh(case, physics, summary)

            self.assertEqual(
                acoustics["surfaces"],
                [
                    {"name": "impeller1", "patches": ["PartSurface_impeller1"]},
                    {"name": "impeller2", "patches": ["PartSurface_impeller2"]},
                    {"name": "laptop_wall", "patches": ["air_to_CPU"]},
                ],
            )
            self.assertAlmostEqual(acoustics["rhoInf"], 1.206)
            self.assertAlmostEqual(acoustics["responseDelay"], 0.01, places=6)
            self.assertEqual(acoustics["lib"], library_name())
            self.assertTrue(acoustics["nonUniformSurfaceMotion"])
            self.assertAlmostEqual(acoustics["CofR"][0], -0.067799999999992)
            self.assertAlmostEqual(acoustics["CofR"][2], 0.081)

            by_name = {item["name"]: item["position"] for item in acoustics["observers"]}
            self.assertAlmostEqual(by_name["mic_px"][0], 0.8)
            self.assertAlmostEqual(by_name["mic_mx"][0], -0.7)
            self.assertAlmostEqual(by_name["mic_py"][1], 0.7)
            self.assertAlmostEqual(by_name["mic_my"][1], -0.6)
            self.assertAlmostEqual(by_name["mic_pz"][2], 0.9)
            self.assertAlmostEqual(by_name["mic_mz"][2], -0.8)
            self.assertLess(by_name["mic_px"][0], 0.3 + 0.5 + 1e-9)
            self.assertGreater(by_name["mic_px"][0], 0.3)

            common = (case / "system" / "commonSettings").read_text(encoding="utf-8")
            head, _, _observers = common.partition("observers")
            self.assertIn("pInf            0;", common)
            self.assertIn("rhoInf          1.206;", common)
            self.assertIn(f'libs            ("{library_name()}");', common)
            self.assertNotIn("pRef", head)
            self.assertIn("pRef", _observers)

            control_surfaces = (case / "system" / "fwhControl").read_text(encoding="utf-8")
            self.assertIn("Farassat1AFormulation", control_surfaces)
            self.assertIn("patches         (PartSurface_impeller1);", control_surfaces)
            self.assertIn("patches         (PartSurface_impeller2);", control_surfaces)
            self.assertIn("patches         (air_to_CPU);", control_surfaces)
            self.assertNotIn(
                "PartSurface_impeller1 PartSurface_impeller2",
                control_surfaces,
            )
            self.assertNotIn("cyclicAMI", control_surfaces)
            self.assertIn("fixedResponseDelay false;", control_surfaces)
            self.assertIn("fixedResponseDelay true;", control_surfaces)
            self.assertIn("nonUniformSurfaceMotion true;", control_surfaces)

            control = (case / "system" / "controlDict").read_text(encoding="utf-8")
            self.assertIn('#include "fwhControl"', control)
            self.assertIn("application     pimpleFoam;", control)
            report = json.loads((case / "setup_report.json").read_text(encoding="utf-8"))
            self.assertEqual(report["mode"], "transient-fwh")
            allrun = (case / "Allrun").read_text(encoding="utf-8")
            self.assertIn("acousticData", allrun)
            commands = [line for line in allrun.splitlines() if not line.startswith("#")]
            self.assertFalse(any("-region" in line for line in commands))

    def test_rejects_a_case_without_an_impeller(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sph = root / "case.sph"
            sph.write_text(_SPH, encoding="utf-8")
            physics = parse_sph(sph)
            case = root / "openfoam"
            ApplyTransientFlowTest()._case(case, with_mrf=True)
            boundary = case / "constant" / "air" / "polyMesh" / "boundary"
            boundary.write_text(
                boundary.read_text(encoding="utf-8").replace(
                    "PartSurface_impeller1", "air_to_cover"
                ),
                encoding="utf-8",
            )
            _write(case / "constant" / "air" / "polyMesh" / "points", _POINTS)
            summary = apply_transient_flow(case, physics)
            with self.assertRaises(ValueError):
                apply_fwh(case, physics, summary)

    def test_binary_point_bounds(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "points"
            header = (
                "FoamFile\n{\n    format      binary;\n    class       vectorField;\n}\n"
                "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
                "2\n("
            ).encode("ascii")
            blob = struct.pack("<3d", -1.5, 0.25, 3.0) + struct.pack("<3d", 4.0, -2.0, 0.5)
            path.write_bytes(header + blob + b")\n")
            lo, hi = mesh_bounds(path)
        self.assertEqual(lo, (-1.5, -2.0, 0.5))
        self.assertEqual(hi, (4.0, 0.25, 3.0))


if __name__ == "__main__":
    unittest.main()
