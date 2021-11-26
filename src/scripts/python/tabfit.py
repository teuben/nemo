#! /usr/bin/env python
#


import sys
import numpy as np

from astropy.io import ascii
import pandas as pd

try:
    from hyperfit.linfit import LinFit
    from hyperfit.data import ExampleData
    Qfit = True
    print("hyperfit available")
except:
    Qfit = False

Qdebug = False    
Qplot = True

table = sys.argv[1]
xcol = 0
ycol = 1
deg  = 2

data = np.loadtxt(table).T
if Qdebug:
    print('TABLE:',data.shape)

x = data[xcol]
y = data[ycol]

rcond = None;
full  = False;
w     = None;
cov   = False;

if full==False and cov==True:
    (p,v) = np.polyfit(x,y,deg, rcond, full, w, cov)
    print('V:',v)
else:
    p     = np.polyfit(x,y,deg, rcond, full, w, cov)
print('P:',p)


if Qplot:
    import matplotlib.pyplot as plt
    xp = np.linspace(x.min(), x.max(), 2*len(x))
    z = np.poly1d(p)
    plt.figure()
    plt.plot(x,y,'.',xp,z(xp),'-')
    plt.show()
