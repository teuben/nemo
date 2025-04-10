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


keyval = [
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
xcol = int(p.get("xcol"))
ycol = int(p.get("ycol"))

# gather data
data = np.loadtxt(infile).T

xdata = data[xcol-1]
ydata = data[ycol-1]

# plot!

plt.figure()
plt.plot(xdata, ydata)
plt.show()
