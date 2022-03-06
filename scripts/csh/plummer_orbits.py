#! /usr/bin/env python
#
#  helper plotting function for the 3 (plummer) orbits
#  


import sys
import numpy as np
import matplotlib.pyplot as plt

run = sys.argv[1]

lines = open('%s.rc' % run).readlines()
for line in lines:
    exec(line)

print("seed=",seed)

o = [0,0,0]

o[0] = np.loadtxt('%s.tab1' % run).T
o[1] = np.loadtxt('%s.tab2' % run).T
o[2] = np.loadtxt('%s.tab3' % run).T

def plot(ax, x, y, size=box):
    ax.plot(x,y)
    ax.scatter(x[0],y[0],color='red')
    ax.set_xlim(-size,size)
    ax.set_ylim(-size,size)
    

fig, axs = plt.subplots(3,3,figsize=(8,8))
fig.suptitle("%s  nbody=%d eps=%g seed=%d"  % (run,nbody,eps,seed))

#   tab2, tab3, tab1

plot(axs[0,0], o[1][0], o[1][1])
plot(axs[0,1], o[1][0], o[1][2])
plot(axs[0,2], o[1][1], o[1][2])

plot(axs[1,0],o[2][0],o[2][1])
plot(axs[1,1],o[2][0],o[2][2])
plot(axs[1,2],o[2][1],o[2][2])


plot(axs[2,0],o[0][0],o[0][1])
plot(axs[2,1],o[0][0],o[0][2])
plot(axs[2,2],o[0][1],o[0][2])


plt.show()
