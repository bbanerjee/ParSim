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
import matplotlib.pyplot as plt
from copy import deepcopy as copy
from matplotlib.widgets import Slider, Button, RadioButtons

def vnorm( x ):
    xx = np.copysign(np.sqrt((x*x).sum(axis=1)[:,np.newaxis]), x[:,0])
    return np.reshape( xx, xx.size )

def pPosition( data, t ):
    px = copy( data[t].get('px',1) )
    px = np.append(px,data[t].get('px',2), axis=0 )
    return px

def pVelocity( data, t ):
    pv = copy( data[t].get('pxI',1) )
    pv = np.append(pv,data[t].get('pxI',2), axis=0 )
    return vnorm(pv)

def pVonMises( data, t ):
    pv = copy( data[t].get('pxI',1) )
    pv = np.append(pv,data[t].get('pxI',2), axis=0 )
    return vnorm(pv)


def plotInteractive( data ):
    ax = plt.subplot(111)
    plt.subplots_adjust(bottom=0.25)
    t = sorted(data.keys())
    
    pl = plt.scatter([0],[0],s=10,linewidths=(0,0,0))
    plt.axis([0, 2, 0, 2])        
    
    axTime = plt.axes([0.25, 0.1, 0.65, 0.03]) 
    sTime = Slider(axTime, 'Time', min(t), max(t), valinit=min(t))
    
    def updateTime(val):
        t0 = sTime.val
        dt = abs(np.array(t)-t0)
        t_in = t[dt.argmin()]
        px = getpx(data,t_in)
        pv = getpv(data,t_in)
        pl.set_offsets(px)  
        pl.set_array(pv)
        plt.draw()
        
    updateTime(min(t))
    sTime.on_changed(updateTime)
        
    plt.show(block=False)