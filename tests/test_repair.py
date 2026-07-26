"""Unit tests for CGNS repair (empty PointList BC deletion)."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

from gph2foam.repair import DATA_DSET, repair_cgns


def _make_minimal_cgns(path: Path) -> None:
    with h5py.File(path, "w") as f:
        f.create_group("CGNSLibraryVersion").create_dataset(
            DATA_DSET, data=np.array([3.21], dtype=np.float32)
        )
        base = f.create_group("Base")
        zone = base.create_group("FluidRegion")
        zone.attrs.create("label", np.bytes_("Zone_t"))
        zone.create_dataset(
            DATA_DSET, data=np.array([[10], [2], [0]], dtype=np.int32)
        )
        zbc = zone.create_group("ZoneBC")
        # Valid BC
        good = zbc.create_group("wall")
        good.create_dataset(
            DATA_DSET, data=np.frombuffer(b"Null", dtype=np.int8)
        )
        pl_good = good.create_group("PointList")
        pl_good.create_dataset(
            DATA_DSET, data=np.array([1, 2, 3], dtype=np.int32)
        )
        # Invalid BC: PointList group without " data"
        bad = zbc.create_group("ghost_open")
        bad.create_dataset(
            DATA_DSET, data=np.frombuffer(b"Null", dtype=np.int8)
        )
        bad.create_group("PointList")


class TestRepair(unittest.TestCase):
    def test_deletes_empty_pointlist_bc(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            src = Path(td) / "in.cgns"
            _make_minimal_cgns(src)
            repair_cgns(src, verbose=False)
            with h5py.File(src, "r") as f:
                keys = list(f["Base/FluidRegion/ZoneBC"].keys())
                self.assertIn("wall", keys)
                self.assertNotIn("ghost_open", keys)
                # PointList reshaped to (n, 1)
                shape = f["Base/FluidRegion/ZoneBC/wall/PointList"][DATA_DSET].shape
                self.assertEqual(shape, (3, 1))
                # BC type forced to BCWall
                raw = f["Base/FluidRegion/ZoneBC/wall"][DATA_DSET][()]
                self.assertEqual(raw.tobytes().decode("ascii"), "BCWall")


if __name__ == "__main__":
    unittest.main()
