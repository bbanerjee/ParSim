import sys
sys.path.append('../')


import numpy as np
from src.datawarehouse import DataWarehouse as Dw
from src.patch import Patch
from src.material import Material
from src.simplecontact import VelocityContact
from src.simplecontact import FrictionlessContact
from src.simplecontact import FrictionContact
from src.simplecontact import FreeContact
from src.material import JacobianError
from src.boundcond import BoundaryCondition as Bc
from src.mpmutils import readableTime as readTime
from src.shape2 import Quad
from src.shape2 import GIMP
from src.shape2 import Linear
from src import geomutils as gu
from src import plotutil
import src.mpm2d as mpm2d