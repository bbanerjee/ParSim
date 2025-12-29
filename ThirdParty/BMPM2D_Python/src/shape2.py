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

import shapes.gimp2 as gimp2
import shapes.quad2 as quad2
import shapes.linear2 as linear2

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
        self.gimp = gimp2

    def updateContribList( self, dw, patch, dwi ):
        self.gimp.updateContribList( dw, patch, dwi )


#===============================================================================
class Quad(Shape):
    def __init__(self, useCython=True):
        self.nSupport = 9
        self.nGhost = 2
        Shape.__init__(self)
        self.quad = quad2

    def updateContribList( self, dw, patch, dwi ):
        self.quad.updateContribList( dw, patch, dwi )
        
#===============================================================================
class Linear(Shape):
    def __init__(self, useCython=True):
        self.nSupport = 4
        self.nGhost = 1
        Shape.__init__(self)
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