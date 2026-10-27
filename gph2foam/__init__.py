"""gph2foam — convert Software Cradle GPH (CRDL-FLD) meshes to OpenFOAM cases.

Pipeline: GPH → CGNS (gphdecoding) → optional HDF5 repair → OpenFOAM (cgns2foam).
"""

from .convert import convert_gph
from .sph import parse_sph

__all__ = ["convert_gph", "parse_sph"]
__version__ = "0.1.0"
