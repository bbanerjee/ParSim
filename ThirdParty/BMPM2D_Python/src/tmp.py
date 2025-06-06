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

def getX():
    x = (np.array(range(401)) - 200.)/100.
    return x

def getquad(X):
    s = []
    g = []
    for x in X:
        r = np.abs(x)
        h = 1.0
        sgn = np.copysign(1.,x)
        
        if ( r < 0.5*h ):
            S = -r*r/(h*h) + 3./4.
            G = -2.*x/(h*h)
        elif ( r < 1.5*h ): 
            S = r*r/(2.*h*h) - 3.*r/(2.*h) + 9./8.
            G = x/(h*h) - sgn*3./(2.*h)        
        else: 
            S = 0.
            G = 0.	
        s.append(S)
        g.append(G)
        
    S = np.array(s)
    G = np.array(g)
    return (S,G)

def getgimp(X):
    s = []
    g = []
    for x in X:

        r = np.abs(x)
        l = .125
        h = 1.
        sgn = np.copysign(1.,x)		    
        if (r < l):
            S = 1. - (r*r+l*l)/(2.*h*l)
            G = -x/(h*l)
        elif (r < h-l):
            S = 1. - r/h
            G = -sgn/h
        elif (r < h+l):
            S = (h+l-r)*(h+l-r) / (4.*h*l)
            G = (h+l-r) / (-2.*sgn*h*l)
        else:
            S = 0.
            G = 0.
        s.append(S)
        g.append(G)
            
    S = np.array(s)
    G = np.array(g)
    return (S,G)    