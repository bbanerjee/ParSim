#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.
#

import numpy as np
import Material
import cProfile
from itertools import izip, count
from inspect import isfunction
import mpmutils as util

try:
    import materialmodel2d_c as mmodel_c
    import mpmutils_c as util_c
except Exception:
    mmodel_c = mmodel
    util_c = util


#===============================================================================
class Material_Stationary(Material):
    # Stationary Rigid Material - 
    def __init__(self, props, model, dwi, shape, useCython=True):
        self.props = props
        self.dwi = dwi
        self.shape = shape
                       
        if useCython:
            self.util = util_c
            self.mmodel = mmodel_c
        else:
            self.util = util
            self.mmodel = mmodel
    
    def getdwi( self ):
        return self.dwi


    def updateContributions( self, dw, patch ):
        dw.zeroGrid( self.dwi )
        self.shape.updateContribList( dw, patch, self.dwi )                
    
    
    def setVelocity( self, dw, v ):
        pass
            
    def setExternalLoad( self, dw, fe ):
        pass
                
                
    def setExternalAcceleration( self, dw, acc ):
        pass


    def applyExternalLoads( self, dw, patch ):
        # Apply external loads to each material
        pass

            
    def interpolateParticlesToGrid( self, dw, patch ):
        # Interpolate particle mass and momentum to the grid
        cIdx,cW = dw.getMult( ['cIdx','cW'], self.dwi )     

        pp = dw.get( 'pm', self.dwi )                          # Mass
        gg = dw.get( 'gm', self.dwi)
        self.util.integrate( cIdx, cW, pp, gg )

     
    def computeStressTensor( self, dw, patch ):    
        pass 
                        
                        
    def computeInternalForce( self, dw, patch ):  
        # Compute internal body forces - integrate divergence of stress to grid
        pass


    def computeAndIntegrateAcceleration( self, dw, patch, tol ):
        # Integrate grid acceleration
        pass
            
            
    def interpolateToParticlesAndUpdate( self, dw, patch ):         
        pass
    