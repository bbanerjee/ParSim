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