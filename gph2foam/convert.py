"""GPH -> CGNS -> repair -> OpenFOAM orchestration."""

from __future__ import annotations

import shutil
import tempfile
import time
from pathlib import Path
from typing import Any

from .deps import import_cgns2foam_convert, import_cgns2foam_writer, import_gph2cgns
from .repair import repair_cgns


def _resolve_regions_json(
    gph_path: Path,
    cgns_path: Path,
    regions_json: str | Path | None,
) -> Path | None:
    """Place a sidecar ``<cgns>.json`` next to the intermediate CGNS if needed."""
    if regions_json is not None:
        src = Path(regions_json)
        if not src.is_file():
            raise FileNotFoundError(f"regions JSON not found: {src}")
    else:
        # Prefer JSON beside the GPH, then beside a pre-existing CGNS stem.
        candidates = [
            gph_path.with_suffix(".json"),
            gph_path.parent / f"{gph_path.stem}_fix.json",
        ]
        src = next((c for c in candidates if c.is_file()), None)
        if src is None:
            return None

    dest = cgns_path.with_suffix(".json")
    if src.resolve() != dest.resolve():
        shutil.copy2(src, dest)
    return dest


def gph_to_cgns(
    gph_path: str | Path,
    cgns_path: str | Path,
    *,
    single_zone: bool = False,
    zone_name: str | None = None,
    skip_fluid_region: bool = True,
    verbose: bool = True,
) -> Path:
    """Parse GPH and write a multi-zone (or single-zone) CGNS/HDF5 file.

    By default ``FluidRegion`` (full-mesh duplicate of all parts) is omitted
    so OpenFOAM mono/CHT cases are not inflated.  Pass
    ``skip_fluid_region=False`` to keep vendor FLDUTIL parity.
    """
    import h5py

    gph2cgns = import_gph2cgns()
    gph_path = Path(gph_path)
    cgns_path = Path(cgns_path)
    cgns_path.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()
    if verbose:
        print(f"[gph2foam] reading GPH: {gph_path}")
    mesh = gph2cgns.parse_gph_mesh(str(gph_path))
    if mesh.get("vertices") is None or mesh.get("link_data") is None:
        raise RuntimeError(f"failed to extract mesh from {gph_path}")

    ld = mesh["link_data"]
    if verbose:
        print(
            f"[gph2foam] mesh: {mesh['n_vertices']} verts, "
            f"{ld['n_faces']} faces, {ld['n_cells']} cells, "
            f"{len(ld['boundary_faces'])} boundary faces "
            f"[{time.perf_counter() - t0:.1f}s]"
        )

    if single_zone:
        zone_names = (zone_name or "FluidRegion",)
        skip_fluid_region = False
    else:
        zone_names = None

    plan = gph2cgns._build_zone_plan(mesh, override_zone_names=zone_names)
    if skip_fluid_region:
        plan = [(n, m) for n, m in plan if n != "FluidRegion"]
    if not plan:
        raise RuntimeError("no zones left to write after FluidRegion filter")

    if verbose:
        print(f"[gph2foam] zones: {len(plan)}")
        for zname, mask in plan:
            print(f"           - {zname}  ({int(mask.sum())} cells)")

    t1 = time.perf_counter()
    if verbose:
        print(f"[gph2foam] writing CGNS: {cgns_path}")

    # gph2cgns.write_cgns has no skip_fluid_region flag; wrap zone plan.
    _orig_plan = gph2cgns._build_zone_plan

    def _filtered_plan(mesh_dict, override_zone_names=None):
        built = _orig_plan(mesh_dict, override_zone_names=override_zone_names)
        if skip_fluid_region:
            built = [(n, m) for n, m in built if n != "FluidRegion"]
        return built

    gph2cgns._build_zone_plan = _filtered_plan  # type: ignore[assignment]
    try:
        gph2cgns.write_cgns(mesh, str(cgns_path), zone_names=zone_names)
    finally:
        gph2cgns._build_zone_plan = _orig_plan  # type: ignore[assignment]

    if verbose:
        print(f"[gph2foam] CGNS done [{time.perf_counter() - t1:.1f}s]")
        if skip_fluid_region:
            # Sanity: FluidRegion should not be present.
            with h5py.File(cgns_path, "r") as f:
                base = f.get("Base")
                if base is not None and "FluidRegion" in base:
                    print("[gph2foam] warning: FluidRegion still present in CGNS")
    return cgns_path


def convert_gph(
    gph_path: str | Path,
    out_dir: str | Path | None = None,
    *,
    cgns_path: str | Path | None = None,
    keep_cgns: bool = False,
    reuse_cgns: bool = False,
    skip_repair: bool = False,
    force_bcwall: bool = True,
    delete_at_groups: bool = False,
    upgrade_int64: bool = False,
    single_zone: bool = False,
    zone_name: str | None = None,
    skip_fluid_region: bool = True,
    cht_direct: bool = False,
    regions_json: str | Path | None = None,
    openfoam_native: bool = False,
    scan_only: bool = False,
    report_path: str | Path | None = None,
    solid_patterns: list[str] | None = None,
    fluid_patterns: list[str] | None = None,
    verbose: bool = True,
) -> Any:
    """Convert a ``.gph`` file to an OpenFOAM case directory.

    Steps:
      1. GPH -> temporary (or ``cgns_path``) CGNS via gphdecoding/gph2cgns
      2. Repair empty/invalid ZoneBC PointLists (cradlecgns logic)
      3. CGNS -> OpenFOAM via cgns2foam (mono or ``--cht-direct``)

    Returns the cgns2foam result (``Mesh`` or ``CouplingReport``).
    """
    gph_path = Path(gph_path).resolve()
    if not gph_path.is_file():
        raise FileNotFoundError(gph_path)

    stem = gph_path.stem
    if out_dir is None:
        out_dir = gph_path.parent / stem
    out_dir = Path(out_dir)

    tmp_dir: tempfile.TemporaryDirectory[str] | None = None

    try:
        if cgns_path is None:
            if keep_cgns:
                cgns_path = out_dir.parent / f"{stem}.cgns"
            else:
                tmp_dir = tempfile.TemporaryDirectory(prefix="gph2foam_")
                cgns_path = Path(tmp_dir.name) / f"{stem}.cgns"
        else:
            keep_cgns = True
        cgns_path = Path(cgns_path)

        reused = False
        if reuse_cgns and cgns_path.is_file():
            reused = True
            if verbose:
                print(f"[gph2foam] reusing CGNS: {cgns_path}")
        else:
            gph_to_cgns(
                gph_path,
                cgns_path,
                single_zone=single_zone,
                zone_name=zone_name,
                skip_fluid_region=skip_fluid_region,
                verbose=verbose,
            )

        if not skip_repair and not reused:
            if verbose:
                print(f"[gph2foam] repairing CGNS: {cgns_path}")
            repair_cgns(
                cgns_path,
                force_bcwall=force_bcwall,
                delete_at_groups=delete_at_groups,
                upgrade_int64=upgrade_int64,
                verbose=verbose,
            )
        elif verbose and reused:
            print("[gph2foam] repair skipped (reused CGNS)")
        elif verbose:
            print("[gph2foam] repair skipped")

        if cht_direct or regions_json is not None:
            placed = _resolve_regions_json(gph_path, cgns_path, regions_json)
            if verbose and placed is not None:
                print(f"[gph2foam] regions sidecar: {placed}")
            elif cht_direct and placed is None:
                raise FileNotFoundError(
                    "CHT mode needs a regions JSON beside the GPH "
                    f"({gph_path.with_suffix('.json')}) or --regions PATH"
                )

        foam_convert = import_cgns2foam_convert()
        writer = import_cgns2foam_writer()
        write_opts = (
            writer.WriteOptions.openfoam_native() if openfoam_native else None
        )

        if scan_only and not cht_direct:
            return foam_convert.scan_file(
                str(cgns_path),
                report_path=str(report_path) if report_path else None,
                verbose=verbose,
                solid_patterns=solid_patterns,
                fluid_patterns=fluid_patterns,
            )

        if verbose:
            mode = "cht-direct" if cht_direct else "mono-block"
            print(f"[gph2foam] OpenFOAM ({mode}) -> {out_dir}")

        result = foam_convert.convert_file(
            str(cgns_path),
            str(out_dir),
            verbose=verbose,
            write_options=write_opts,
            cht_direct=cht_direct,
            solid_patterns=solid_patterns,
            fluid_patterns=fluid_patterns,
        )

        if cht_direct and report_path:
            src = Path(out_dir) / "coupling_scan.json"
            if src.is_file():
                Path(report_path).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, report_path)
                if verbose:
                    print(f"[gph2foam] report copied to {report_path}")

        if keep_cgns and verbose:
            print(f"[gph2foam] CGNS kept at {cgns_path}")

        return result
    finally:
        if tmp_dir is not None:
            tmp_dir.cleanup()
