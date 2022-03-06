#! /usr/bin/env python
#
#  helper plotting function for the 3 (plummer) orbits
#  


import os
import sys
import numpy as np
import matplotlib.pyplot as plt

#  helper
if len(sys.argv) == 1:
    print("Usage: %s parameter_file [key=val]" % sys.argv[0])
    print("  parameter_file    set of (python) variable settings for the plotting script")
    print("  key=val           optionally override settings from the parameter_file")
    sys.exit(0)
    

#  read parameter file
run = sys.argv[1]
lines = open('%s' % run).readlines()
for line in lines:
    exec(line)

#  optionally override parameters
for arg in sys.argv[2:]:
    exec(arg)

plname = "%s_%d_%g_%d.%s" % (run, nbody, eps, seed, fig)
print('PLOT:',plname)

o = [0,0,0]

o[0] = np.loadtxt('%s.tab1' % run).T
o[1] = np.loadtxt('%s.tab2' % run).T
o[2] = np.loadtxt('%s.tab3' % run).T

def plot(ax, x, y, size=box):
    ax.plot(x,y)
    ax.scatter(x[0],y[0],color='red')
    ax.set_xlim(-size,size)
    ax.set_ylim(-size,size)

if fig == "x":
    interactive = True
else:
    interactive = False

    
if interactive:
    plt.ion()

fig, axs = plt.subplots(3,3,figsize=(8,8))
fig.suptitle("%s  nbody=%d eps=%g seed=%d"  % (run,nbody,eps,seed))

#   tab2, tab3, tab1

#   tab2: analytical orbit
plot(axs[0,0], o[1][0], o[1][1])    
plot(axs[0,1], o[1][0], o[1][2])
plot(axs[0,2], o[1][1], o[1][2])

#   tab3: frozen nbody orbit
plot(axs[1,0],o[2][0],o[2][1])
plot(axs[1,1],o[2][0],o[2][2])
plot(axs[1,2],o[2][1],o[2][2])

#   tab1: nbody orbit
plot(axs[2,0],o[0][0],o[0][1])
plot(axs[2,1],o[0][0],o[0][2])
plot(axs[2,2],o[0][1],o[0][2])

if interactive:
    plt.ioff()
    plt.show()
else:
    plt.savefig(plname)
