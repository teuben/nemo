#! /usr/bin/env python
#
#  tabplot in python

import os
import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from nemopy import getparam

#import pdb
#pdb.set_trace()


keyval = [  # list of key values
    "in=???\n       input fits cube",
    "xcol=1\n       X column(s) (1=first column)",
    "ycol=2\n       Y column(s)",
    "xmin=\n        X-min value if not autoscaled",
    "xmax=\n        X-max value if not autoscaled",
    "ymin=\n        Y-min value if not autoscaled []",
    "ymax=\n        Y-max value if not autoscaled []",
    "line=\n        TBD",
    "color=\n       TBD",
    "point=\n       TBD",
    "VERSION=0.1\n  3-apr-2023 PJT",
]

usage = """  
plot a table like tabplot
"""

p = getparam.Param(keyval,usage)

# get CLI parameters
infile = p.get("in")
xinput = p.get("xcol")
yinput = p.get("ycol")

# split string to get ints for xcol and ycol
xcols = [int(x) for x in xinput.split(',')]
ycols = [int(y) for y in yinput.split(',')]

# if no xcol or ycol specified, assume col 1 is x and col 2 is y
if len(xcols) == 0:
    xcols = [1]
if len(ycols) == 0:
    ycols = [2]

# gather data
data = np.loadtxt(infile).T

xdata = [0] * len(xcols)
ydata = [0] * len(ycols)

for i in range(len(xcols)): # Read x data
    xdata[i] = data[xcols[i]-1]

for i in range(len(ycols)): # Read y data
    ydata[i] = data[ycols[i]-1]

print(xdata[0])
print(ydata[0])
# plot!

plt.figure()
#plt.plot(xdata, ydata)

if len(xdata) == 1: # Case: only 1 xcol, plot each y against the only x
    for i in range (len(ydata)):
        plt.plot(xdata[0], ydata[i])
else:               # Case: more than 1 xcol, plot each ycol with its xcol at the same index
    for i in range (len(xdata)):
        try:
            plt.plot(xdata[i], ydata[i])
        except IndexError:
            print("ERROR: xcols and ycols don't match\n")
            sys.exit()

plt.show()
