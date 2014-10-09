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