import sys
sys.path.append('../')


import numpy as np
from src.datawarehouse import DataWarehouse as Dw
from src.patch import Patch
from src.material import Material
from src.contact import VelocityContact
from src.contact import FrictionlessContact
from src.contact import FrictionContact
from src.contact import FreeContact
from src.material import JacobianError
from src.boundcond import BoundaryCondition as Bc
from src.mpmutils import readableTime as readTime
from src.shape2 import Quad
from src.shape2 import GIMP
from src.shape2 import Linear
from src import geomutils as gu
from src import plotutil
import src.mpm2d as mpm2d