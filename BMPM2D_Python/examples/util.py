import numpy as np

def getMomentum( data, dwi ):
    t = data.keys()
    t.sort()
    dt = t[1]-t[0]
    pw = []
    pf = []
    for t0 in t:
        pw0 = data[t0].dw['pw',dwi]
        pm0 = data[t0].dw['pm',dwi]
        pa0 = data[t0].dw['pvI',dwi]
        pf0 = pa0*pm0*dt
        mpw0 = np.sqrt((pw0*pw0).sum(axis=1)[:,np.newaxis])
        mpf0 = np.sqrt((pf0*pf0).sum(axis=1)[:,np.newaxis])
        pw.append(sum(pw0))
        pf.append(sum(mpf0))
        
    t = np.array(t)
    pw = np.array(pw)
    pf = np.array(pf)
    return t,pw,pf