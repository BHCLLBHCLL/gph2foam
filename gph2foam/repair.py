"""HDF5 CGNS repair steps for Cradle/GPH exports (from cradlecgns).

GPH→CGNS writes ZoneBC entries whose PointList group has no ``" data"``
dataset for surface regions absent from a zone.  Those empty BCs crash
cgns2foam's reader and must be deleted before OpenFOAM conversion.

Default (foam-safe) repair:
  1. CGNSLibraryVersion → 4.2
  2. ZoneBC type → BCWall  (GPH emits ``Null``; name rules still override openings)
  3. PointList shape (n,) → (n, 1)
  5. Delete BCs with PointList but no ``" data"``
  6. Zone_t data shape (1, 3) → (3, 1)
  7. Create empty ElementConnectivity ``" data"`` when missing

Optional:
  4. Delete ``@…`` groups (``delete_at_groups=True``)
  8. int32 → int64 for Element*/Zone data (``upgrade_int64=True``)
"""

from __future__ import annotations

import shutil
from pathlib import Path

import h5py
import numpy as np

DATA_DSET = " data"
TARGET_VERSION = 4.2
BC_TYPE_WALL = "BCWall"


def _find_zonebc_groups(f: h5py.File, root: str = "/") -> list[h5py.Group]:
    zonebcs: list[h5py.Group] = []
    node = f if root == "/" else f[root]
    if not isinstance(node, h5py.Group):
        return zonebcs
    for key in node.keys():
        if key == "ZoneBC":
            zonebcs.append(node[key])
        elif isinstance(node[key], h5py.Group):
            sub = (root if root != "/" else "") + "/" + key
            zonebcs.extend(_find_zonebc_groups(f, sub))
    return zonebcs


def step_set_version(f: h5py.File, dry_run: bool, log) -> None:
    name = "CGNSLibraryVersion"
    version_data = np.array([TARGET_VERSION], dtype=np.float32)
    if name not in f:
        log(f"[1 version] create {name} = {TARGET_VERSION}")
        if not dry_run:
            grp = f.create_group(name)
            grp.create_dataset(DATA_DSET, data=version_data)
        return
    grp = f[name]
    if DATA_DSET not in grp:
        log(f"[1 version] write {name} = {TARGET_VERSION}")
        if not dry_run:
            grp.create_dataset(DATA_DSET, data=version_data)
        return
    current = float(grp[DATA_DSET][()].flat[0])
    if abs(current - TARGET_VERSION) < 1e-5:
        log(f"[1 version] already {TARGET_VERSION}")
        return
    log(f"[1 version] {current} -> {TARGET_VERSION}")
    if not dry_run:
        del grp[DATA_DSET]
        grp.create_dataset(DATA_DSET, data=version_data)


def step_bc_to_bcwall(f: h5py.File, dry_run: bool, log) -> None:
    to_fix: list[tuple[h5py.Dataset, str]] = []

    def visit(name, obj):
        if not isinstance(obj, h5py.Dataset):
            return
        if "ZoneBC" not in name or not name.endswith("/" + DATA_DSET):
            return
        parent = obj.parent
        if "PointList" in parent.name or "GridLocation" in parent.name:
            return
        to_fix.append((obj, name))

    f.visititems(visit)
    if not to_fix:
        log("[2 BCWall] nothing to update")
        return
    val = BC_TYPE_WALL.encode("ascii")[:6]
    arr = np.frombuffer(val, dtype=np.int8)
    for ds, path in to_fix:
        log(f"[2 BCWall] {path}")
        if not dry_run:
            parent = ds.parent
            key = ds.name.split("/")[-1]
            del parent[key]
            parent.create_dataset(key, data=arr)


def step_fix_pointlist_shape(f: h5py.File, dry_run: bool, log) -> None:
    to_fix: list[tuple[h5py.Group, str]] = []

    def visit(name, obj):
        if not isinstance(obj, h5py.Group):
            return
        if not name.endswith("/PointList"):
            return
        if DATA_DSET not in obj:
            return
        ds = obj[DATA_DSET]
        if len(ds.shape) == 1:
            to_fix.append((obj, name))

    f.visititems(visit)
    if not to_fix:
        log("[3 PointList] nothing to reshape")
        return
    for pl_grp, path in to_fix:
        data = np.asarray(pl_grp[DATA_DSET][()])
        log(f"[3 PointList] {path}: {data.shape} -> ({data.size}, 1)")
        if not dry_run:
            del pl_grp[DATA_DSET]
            pl_grp.create_dataset(
                DATA_DSET, data=data.reshape(-1, 1), dtype=data.dtype
            )


def step_delete_at_groups(f: h5py.File, dry_run: bool, log) -> None:
    to_delete: list[str] = []

    def visit(name, obj):
        if isinstance(obj, h5py.Group) and name.rsplit("/", 1)[-1].startswith("@"):
            to_delete.append(name)

    f.visititems(visit)
    to_delete.sort(key=lambda p: p.count("/"), reverse=True)
    if not to_delete:
        log("[4 @group] none found")
        return
    for path in to_delete:
        log(f"[4 @group] delete {path}")
        if not dry_run:
            del f[path]


def step_delete_empty_pointlist_bc(f: h5py.File, dry_run: bool, log) -> None:
    """Delete ZoneBC children that have PointList but no ``\" data\"`` dataset."""
    to_delete: list[tuple[h5py.Group, str]] = []
    for zonebc in _find_zonebc_groups(f):
        for bc_key in list(zonebc.keys()):
            bc_grp = zonebc[bc_key]
            if not isinstance(bc_grp, h5py.Group):
                continue
            if "PointList" not in bc_grp:
                continue
            pointlist = bc_grp["PointList"]
            if not isinstance(pointlist, h5py.Group):
                continue
            if DATA_DSET not in pointlist:
                to_delete.append((zonebc, bc_key))
    if not to_delete:
        log("[5 empty-BC] none found")
        return
    for zonebc, bc_key in to_delete:
        log(f"[5 empty-BC] delete {zonebc.name}/{bc_key}")
        if not dry_run:
            del zonebc[bc_key]


def step_fix_zone_data_shape(f: h5py.File, dry_run: bool, log) -> None:
    to_fix: list[h5py.Group] = []

    def visit(name, obj):
        if not isinstance(obj, h5py.Group):
            return
        if obj.attrs.get("label") != b"Zone_t":
            return
        if DATA_DSET not in obj:
            return
        if tuple(obj[DATA_DSET].shape) == (1, 3):
            to_fix.append(obj)

    f.visititems(visit)
    if not to_fix:
        log("[6 Zone] nothing to reshape")
        return
    for zone_grp in to_fix:
        flat = np.asarray(zone_grp[DATA_DSET][()], dtype=np.int32).reshape(-1)
        if flat.size != 3:
            continue
        log(f"[6 Zone] {zone_grp.name}: (1,3) -> (3,1)")
        if not dry_run:
            del zone_grp[DATA_DSET]
            zone_grp.create_dataset(
                DATA_DSET, data=flat.reshape(3, 1), dtype=np.int32
            )


def step_fix_empty_element_connectivity(f: h5py.File, dry_run: bool, log) -> None:
    to_fix: list[h5py.Group] = []

    def visit(name, obj):
        if not isinstance(obj, h5py.Group):
            return
        if obj.attrs.get("label") != b"DataArray_t":
            return
        if not name.endswith("/ElementConnectivity"):
            return
        if DATA_DSET not in obj:
            to_fix.append(obj)

    f.visititems(visit)
    if not to_fix:
        log("[7 EmptyConn] none found")
        return
    for conn_grp in to_fix:
        log(f"[7 EmptyConn] create empty data under {conn_grp.name}")
        if not dry_run:
            conn_grp.create_dataset(DATA_DSET, data=np.array([], dtype=np.int32))


def step_upgrade_int64(f: h5py.File, dry_run: bool, log) -> None:
    to_convert: list[tuple[h5py.Dataset, str]] = []

    def visit(name, obj):
        if not isinstance(obj, h5py.Dataset):
            return
        if obj.dtype not in (np.int32, np.dtype("int32")):
            return
        if not name.endswith("/" + DATA_DSET):
            return
        if "PointList" in name:
            return
        if any(
            k in name
            for k in ("ElementConnectivity", "ElementRange", "ElementStartOffset")
        ):
            to_convert.append((obj, name))
            return
        parent = obj.parent
        if isinstance(parent, h5py.Group) and parent.attrs.get("label") == b"Zone_t":
            to_convert.append((obj, name))

    f.visititems(visit)
    seen: set[str] = set()
    unique: list[tuple[h5py.Dataset, str]] = []
    for ds, path in to_convert:
        if path in seen:
            continue
        seen.add(path)
        unique.append((ds, path))
    if not unique:
        log("[8 int64] nothing to upgrade")
        return
    for ds, path in unique:
        log(f"[8 int64] {path}")
        if not dry_run:
            data = np.asarray(ds[()], dtype=np.int64)
            parent = ds.parent
            key = ds.name.split("/")[-1]
            del parent[key]
            parent.create_dataset(key, data=data, dtype=np.int64)
            if parent.attrs.get("type", b"") == b"I4":
                del parent.attrs["type"]
                parent.attrs.create(
                    "type", np.bytes_("I8"), dtype=h5py.string_dtype(length=3)
                )


def repair_cgns(
    input_path: str | Path,
    output_path: str | Path | None = None,
    *,
    dry_run: bool = False,
    force_bcwall: bool = True,
    delete_at_groups: bool = False,
    upgrade_int64: bool = False,
    verbose: bool = True,
) -> Path:
    """Repair a CGNS/HDF5 file in place or write a copy.

    Returns the path of the repaired file (or the input path on dry-run).
    """
    input_path = Path(input_path)
    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    target = input_path
    if output_path is not None and not dry_run:
        target = Path(output_path)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(input_path, target)

    def log(msg: str) -> None:
        if verbose:
            print(f"[gph2foam:repair] {msg}")

    mode = "r" if dry_run else "r+"
    with h5py.File(target, mode) as f:
        step_set_version(f, dry_run, log)
        if force_bcwall:
            step_bc_to_bcwall(f, dry_run, log)
        else:
            log("[2 BCWall] skipped (--keep-bc-types)")
        step_fix_pointlist_shape(f, dry_run, log)
        if delete_at_groups:
            step_delete_at_groups(f, dry_run, log)
        else:
            log("[4 @group] skipped")
        step_delete_empty_pointlist_bc(f, dry_run, log)
        step_fix_zone_data_shape(f, dry_run, log)
        step_fix_empty_element_connectivity(f, dry_run, log)
        if upgrade_int64:
            step_upgrade_int64(f, dry_run, log)
        else:
            log("[8 int64] skipped")

    return target
