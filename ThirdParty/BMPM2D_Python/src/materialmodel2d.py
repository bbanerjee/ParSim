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
from numba import njit

def makeArray( props, modelName ):
    if( modelName == 'planeStrainNeoHookean' ):
        arr = np.zeros(3)
        arr[0] = props['modulus']
        arr[1] = props['poisson']
        arr[2] = props['density']
    elif( modelName == 'planeStrainNeoHookeanMaxStress' ):        
        arr = np.zeros(4)
        arr[0] = props['modulus']
        arr[1] = props['poisson']
        arr[2] = props['density']
        arr[3] = props['maxStress']
    else:
        arr = np.array([])
    return arr

@njit(cache=True)
def planeStrainNeoHookean( props, F ):
    S = np.zeros((2,2))
    v = props[1]
    E = props[0]
    l = E * v / ((1.+v)*(1.-2.*v))
    m = 0.5 * E / (1.+v)
    Ja = F[0,0]*F[1,1] - F[1,0]*F[0,1]
    
    FFT00 = F[0,0]*F[0,0] + F[0,1]*F[0,1] - 1.0
    FFT01 = F[0,0]*F[1,0] + F[0,1]*F[1,1]
    FFT10 = F[1,0]*F[0,0] + F[1,1]*F[0,1]
    FFT11 = F[1,0]*F[1,0] + F[1,1]*F[1,1] - 1.0
    
    logJa = np.log(Ja)
    
    S[0,0] = l*logJa/Ja + m/Ja * FFT00
    S[0,1] = m/Ja * FFT01
    S[1,0] = m/Ja * FFT10
    S[1,1] = l*logJa/Ja + m/Ja * FFT11
    
    return S, Ja

@njit(cache=True)
def planeStrainNeoHookeanMaxStress( props, F ):
    S = np.zeros((2,2))
    v = props[1]
    E = props[0]
    sMax = props[3]
    l = E * v / ((1.+v)*(1.-2.*v))
    m = 0.5 * E / (1.+v)
    Ja = F[0,0]*F[1,1] - F[1,0]*F[0,1]
    
    FFT00 = F[0,0]*F[0,0] + F[0,1]*F[0,1] - 1.0
    FFT01 = F[0,0]*F[1,0] + F[0,1]*F[1,1]
    FFT10 = F[1,0]*F[0,0] + F[1,1]*F[0,1]
    FFT11 = F[1,0]*F[1,0] + F[1,1]*F[1,1] - 1.0
    
    logJa = np.log(Ja)
    
    S[0,0] = l*logJa/Ja + m/Ja * FFT00
    S[0,1] = m/Ja * FFT01
    S[1,0] = m/Ja * FFT10
    S[1,1] = l*logJa/Ja + m/Ja * FFT11
    
    vm = np.sqrt( S[0,0]*S[0,0] - S[0,0]*S[1,1] + S[1,1]*S[1,1] + 3.0*S[1,0]*S[0,1] )     
    if vm > sMax:
        S[0,0] = 0.
        S[0,1] = 0.
        S[1,0] = 0.
        S[1,1] = 0.
        Ja = 1.
        
    return S, Ja

class MaterialModel:
    def __init__(self, modelName, props):
        self.modelName = modelName
        self.props = makeArray(props, modelName)
        
    def getStress( self, F ):
        if self.modelName == 'planeStrainNeoHookean':
            return planeStrainNeoHookean(self.props, F)
        elif self.modelName == 'planeStrainNeoHookeanMaxStress':
            return planeStrainNeoHookeanMaxStress(self.props, F)
        else:
             # Just in case there are other models not covered here but present in legacy code
            raise Exception("Model not supported in optimized version: " + self.modelName)

    def changeProps( self, props ):
        self.props = makeArray(props, self.modelName)
