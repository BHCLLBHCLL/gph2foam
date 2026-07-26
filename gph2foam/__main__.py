"""CLI: ``python -m gph2foam <input.gph> [output_dir]``."""

from __future__ import annotations

import argparse
import os
import sys

from .convert import convert_gph


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="gph2foam",
        description=(
            "Convert a Software Cradle GPH (CRDL-FLD) mesh to an OpenFOAM "
            "case. Pipeline: GPH -> CGNS (gphdecoding) -> repair empty/invalid "
            "ZoneBCs -> OpenFOAM (cgns2foam)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  %(prog)s mesh.gph out_case
  %(prog)s mesh.gph out_case --keep-cgns
  %(prog)s --scan mesh.gph
  %(prog)s --cht-direct mesh.gph out_cht --regions mesh.json
  %(prog)s mesh.gph --cgns-only mesh.cgns

environment:
  GPH2FOAM_GPHDECODING   path to gphdecoding root (gph2cgns.py)
  GPH2FOAM_CGNS2FOAM     path to cgns2foam root (src/convert.py)
""".rstrip(),
    )
    p.add_argument("gph_file", help="Path to the input .gph file")
    p.add_argument(
        "output_dir",
        nargs="?",
        default=None,
        help=(
            "OpenFOAM case directory (default: "
            "<dirname-of-input>/<basename>/)."
        ),
    )
    p.add_argument(
        "-q", "--quiet", action="store_true", help="Suppress progress messages."
    )
    p.add_argument(
        "--keep-cgns",
        action="store_true",
        help="Keep the intermediate .cgns next to the output case (or --cgns).",
    )
    p.add_argument(
        "--cgns",
        metavar="PATH",
        default=None,
        help="Write intermediate CGNS at PATH (implies keep).",
    )
    p.add_argument(
        "--reuse-cgns",
        action="store_true",
        help="If --cgns PATH already exists, skip GPH parse and repair.",
    )
    p.add_argument(
        "--cgns-only",
        metavar="PATH",
        default=None,
        help="Only write repaired CGNS to PATH; skip OpenFOAM conversion.",
    )
    p.add_argument(
        "--skip-repair",
        action="store_true",
        help="Skip CGNS repair (empty PointList BC deletion, etc.).",
    )
    p.add_argument(
        "--keep-bc-types",
        action="store_true",
        help="Do not force ZoneBC types to BCWall (keep GPH Null/etc.).",
    )
    p.add_argument(
        "-d",
        "--delete-at-groups",
        action="store_true",
        help="Also delete HDF5 groups whose names start with '@'.",
    )
    p.add_argument(
        "-i",
        "--int64",
        action="store_true",
        help="Upgrade Element*/Zone int32 datasets to int64.",
    )
    p.add_argument(
        "--single-zone",
        action="store_true",
        help="Emit a single FluidRegion zone instead of GPH partitions.",
    )
    p.add_argument(
        "--keep-fluid-region",
        action="store_true",
        help=(
            "Keep the full-mesh FluidRegion zone (default: omit it so "
            "OpenFOAM cases are not duplicated)."
        ),
    )
    p.add_argument(
        "--zone",
        default=None,
        help="Zone name when --single-zone is set (default: FluidRegion).",
    )
    p.add_argument(
        "--openfoam-native",
        action="store_true",
        help="Write OpenFOAM-native binary polyMesh (default: ANSA-compatible).",
    )
    p.add_argument(
        "--scan",
        action="store_true",
        help="Scan fluid/solid couplings only (no polyMesh write).",
    )
    p.add_argument(
        "--cht-direct",
        action="store_true",
        help=(
            "Write a multi-region chtMultiRegionSimpleFoam case. "
            "Requires a regions JSON (--regions or <gph>.json)."
        ),
    )
    p.add_argument(
        "--regions",
        metavar="PATH",
        default=None,
        help="Regions/physics sidecar JSON for --cht-direct (or coupling scan).",
    )
    p.add_argument(
        "--report",
        metavar="PATH",
        default=None,
        help="Write coupling scan JSON to PATH.",
    )
    p.add_argument(
        "--solid-pattern",
        action="append",
        default=None,
        metavar="REGEX",
        help="Regex for solid zone names (repeatable).",
    )
    p.add_argument(
        "--fluid-pattern",
        action="append",
        default=None,
        metavar="REGEX",
        help="Regex for fluid zone names (repeatable).",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    gph_path = args.gph_file
    if not os.path.isfile(gph_path):
        print(f"error: not a file: {gph_path}", file=sys.stderr)
        return 2

    verbose = not args.quiet

    # CGNS-only path: GPH → CGNS → repair, then exit.
    if args.cgns_only:
        from .convert import gph_to_cgns
        from .repair import repair_cgns

        cgns_out = args.cgns_only
        gph_to_cgns(
            gph_path,
            cgns_out,
            single_zone=args.single_zone,
            zone_name=args.zone,
            skip_fluid_region=not args.keep_fluid_region,
            verbose=verbose,
        )
        if not args.skip_repair:
            repair_cgns(
                cgns_out,
                force_bcwall=not args.keep_bc_types,
                delete_at_groups=args.delete_at_groups,
                upgrade_int64=args.int64,
                verbose=verbose,
            )
        return 0

    cgns_path = args.cgns
    keep_cgns = args.keep_cgns or (cgns_path is not None)

    convert_gph(
        gph_path,
        args.output_dir,
        cgns_path=cgns_path,
        keep_cgns=keep_cgns,
        reuse_cgns=args.reuse_cgns,
        skip_repair=args.skip_repair,
        force_bcwall=not args.keep_bc_types,
        delete_at_groups=args.delete_at_groups,
        upgrade_int64=args.int64,
        single_zone=args.single_zone,
        zone_name=args.zone,
        skip_fluid_region=not args.keep_fluid_region,
        cht_direct=args.cht_direct,
        regions_json=args.regions,
        openfoam_native=args.openfoam_native,
        scan_only=args.scan,
        report_path=args.report,
        solid_patterns=args.solid_pattern,
        fluid_patterns=args.fluid_pattern,
        verbose=verbose,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
