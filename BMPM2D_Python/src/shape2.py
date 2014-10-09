import numpy as np

import shapes.gimp2 as gimp2
try: import shapes.gimp2_c as gimp2_c
except Exception: gimp2_c = gimp2

import shapes.quad2 as quad2
try: import shapes.quad2_c as quad2_c
except Exception: quad2_c = quad2

import shapes.linear2 as linear2
try: import shapes.linear2_c as linear2_c
except Exception: linear2_c = linear2

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

    def updateContribList( self, dw, patch, dwi ):
	self.gimp.updateContribList( dw, patch, dwi )


#===============================================================================
class Quad(Shape):
    def __init__(self, useCython=True):
	self.nSupport = 9
	self.nGhost = 2
	Shape.__init__(self)
	
	if useCython:
	    self.quad = quad2_c
	else:
	    self.quad = quad2

    def updateContribList( self, dw, patch, dwi ):
	self.quad.updateContribList( dw, patch, dwi )
	
#===============================================================================
class Linear(Shape):
    def __init__(self, useCython=True):
	self.nSupport = 4
	self.nGhost = 1
	Shape.__init__(self)
	
	if useCython:
	    self.linear = linear2_c
	else:
	    self.linear = linear2

    def updateContribList( self, dw, patch, dwi ):
	self.linear.updateContribList( dw, patch, dwi )
	
#===============================================================================
class Cubic(Shape):
    def __init__(self, useCython=True):
	self.nSupport = 12
	self.nGhost = 2
	Shape.__init__(self)
	
	if useCython:
	    self.cubic = cubic2_c
	else:
	    self.cubic = cubic2

    def updateContribList( self, dw, patch, dwi ):
	self.cubic.updateContribList( dw, patch, dwi )