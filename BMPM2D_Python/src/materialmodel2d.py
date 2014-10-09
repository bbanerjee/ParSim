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
class MaterialModel:
    # Defines material models - accessed using getStress 
    #  - actual computation done in static methods   
    # Returns stress tensor and jacobian of deformation
    def __init__(self, modelName, props):
        self.modelName = modelName               # Selects Material Model
        self.props = props
        
    def getStress( self, F ):
        model = getattr( self, self.modelName )
        S,Ja = model(self.props, F);    
        return (S,Ja)
    
    def changeProps( self, props ):
        self.props = props


    @staticmethod
    def planeStrainNeoHookean( props, F ):
        # Props - poisson, E
        I2 = F*0.
        I2[0,0] = I2[1,1] = 1.
        v = props['poisson']
        E = props['modulus']
        l = E * v / ((1.+v)*(1.-2.*v))
        m = 0.5 * E / (1.+v)
        Ja = F[0,0]*F[1,1] - F[1,0]*F[0,1]
        S = I2*l*np.log(Ja)/Ja + m/Ja * (np.dot(F, F.T) - I2)
        
        return (S,Ja)
    
    @staticmethod
    def planeStrainNeoHookeanMaxStress( props, F ):
        # Props - poisson, E, maxStress
        I2 = np.eye(2)
        v = props['poisson']
        E = props['modulus']
        sMax = props['maxStress']
        l = E * v / ((1.+v)*(1.-2.*v))
        m = 0.5 * E / (1.+v)
        Ja = np.linalg.det(F)
        S = I2*l*np.log(Ja)/Ja + m/Ja * (np.dot(F, F.T) - I2)
        if vonMises(S) > sMax: 
            S = I2*0.
            Ja = 1.
        
        return (S,Ja)
    
    @staticmethod
    def vonMises( S ):
        return np.sqrt( S[0,0]*S[0,0] - S[0,0]*S[1,1] + S[1,1]*S[1,1] +
            3.*S[1,0]*S[0,1] )