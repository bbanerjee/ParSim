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


#===============================================================================
class BoundaryCondition:
    def __init__(self, bc_type, bc_val, bc_var, dwis, fun ):
        # Set boundary condition - bc_type = 'X' or 'Y'
        # bc_val = value of x or y where condition is applied
        # bc_var is nodal variable to set
        # fun is function - takes a point as input
        self.bc_type = bc_type
        self.bc_val = bc_val
        self.bc_var = bc_var
        self.fun = fun
        self.dwis = dwis
        
    def setBoundCond( self, dw, patch, tol ):
        if( self.bc_type == 'X' ):
            self.bcX( dw, patch, tol )
        else:
            self.bcY( dw, patch, tol )
        
    def bcX( self, dw, patch, tol ):
        #  Set boundary condition on line x=val
        for dwi in self.dwis:
            gg = dw.get( self.bc_var, dwi )
            gx = dw.get( 'gx', dwi )
            for ii in range(len(gx)):
                if( np.abs(gx[ii][0]-self.bc_val) < tol ):
                    gg[ii] = self.fun( gx[ii] )
                
    def bcY( self, dw, patch, tol ):
        #  Set boundary condition on line y=val
        for dwi in self.dwis:
            gg = dw.get( self.bc_var, dwi )
            gx = dw.get( 'gx', dwi )
            for ii in range(len(gx)):
                if( np.abs(gx[ii][1]-self.bc_val) < tol ):
                    gg[ii] = self.fun( gx[ii] )                