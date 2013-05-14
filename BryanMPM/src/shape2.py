import numpy as np
import gimp2
try:
    import gimp2_c
except Exception:
    gimp2_c = gimp

#===============================================================================
class Shape:
    #  Shape functions - compute nodal contributions to particle values
    def __init__(self):
	self.dim = 2;
	self.S = np.zeros([self.dim,1])    # Value of Shape function
	self.G = np.zeros([self.dim,1])    # Value of Shape function derivative	


#===============================================================================
class GIMP(Shape):
    def __init__(self, useCython=True):
	self.nSupport = 9
	self.nGhost = 2
	Shape.__init__(self)
	
	if useCython:
	    self.gimp = gimp2_c
	else:
	    self.gimp = gimp2

    def updateContribList( self, dw, patch, mIdx ):
	self.gimp.updateContribList( dw, patch, mIdx )