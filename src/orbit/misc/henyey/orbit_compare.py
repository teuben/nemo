#! /usr/bin/env python
#
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# perorb
tab1 = 'tab4'
data1 = np.loadtxt(tab1).T
print(data1.shape)

# henyey
tab2 = 'tab1x'
data2 = np.loadtxt(tab2).T
print(data2.shape)

tab3 = 'tab1y'
data3 = np.loadtxt(tab3).T
print(data3.shape)

#plt.figure(dpi=300,figsize=(20/2.54,20/2.54))
plt.figure()
# X
plt.plot(data1[0],data1[1],label='perorb X')
plt.plot(data2[2],data2[5],label='henyey X')
# Y
plt.plot(data1[2],-data1[3],label='perorb Y')
plt.plot(data3[3],-data3[4],label='henyey Y')

plt.xlabel('X,Y')
plt.ylabel('VY,VX')
plt.title('comparing rotation curves of periodic orbits')
plt.legend()
#plt.savefig('pyplot.png')
plt.show()
