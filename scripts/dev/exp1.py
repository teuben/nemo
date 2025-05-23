#! /usr/bin/env python
#
#  Example usage:   ./exp1.py  run5/p41.cgs/snap.mrad run5/p42.cgs/snap.mrad
#
#  the file label will have the ".cgs/snap.mrad" removed
#
import os
import sys
# import docopt
import numpy as np
import matplotlib.pyplot as plt


data = []
tabs = sys.argv[1:]
irad = [8,9,10]         # 8=70% 9=%80 10=90%

for t in tabs:
    data.append(np.loadtxt(t).T)


#plt.figure(dpi=300,figsize=(20/2.54,20/2.54))
plt.figure()

for (d,f) in zip(data,tabs):
    fs = f[:f.index(".cgs")]
    for i in irad:
        plt.plot(d[0],d[i],label=fs+"[%d]"%i)

plt.xlabel('Time')
plt.ylabel('Lagrangian Radius')
plt.title('Mass Radii')
plt.legend()
plt.savefig('exp1.png')
plt.show()
