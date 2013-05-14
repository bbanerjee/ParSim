import sys
sys.path.append('../')

import numpy as np
from src.datawarehouse import DataWarehouse as Dw
from src.patch import Patch
from src.material import Material
from src.material import JacobianError
from src.boundcond import BoundaryCondition as Bc
from src.mpmutils import readableTime as readTime
from src.shape2 import GIMP as Shape    
from src import geomutils
import src.mpm2d as mpm

try:
    from src.shape2_c import GIMP as Shape_c
except Exception:
    Shape_c = Shape