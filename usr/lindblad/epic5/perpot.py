#! /usr/bin/env python
#
#   Example of processing the annoying fortran wrapped "perpot" table
#
#   @todo    convert from/to nemo snapfour style output, some plotting

import numpy as np

#
#           NPER DRP MMAX NTH NC
#           AMPC(M), (M=1,10)
#           AMPS(M), (M=1,10)
#           for M=1,MMAX:
#              CP(I,M) (I=1,NPER)
#              SP(I,M) (I=1,NPER)

#   c[0][i],s[0][i],  c[1][i], s[1][i]    cs[mmax][2][nper]
#

help = """
     500    0.0505000          10          64          10
       0.5       0.5       0.5       0.5       0.5
       0.5       0.5       0.5       0.5       0.5
       -0.5       -0.5       -0.5       -0.5       -0.5
       -0.5       -0.5       -0.5       -0.5       -0.5
      -0.0000000      -0.0000000      -0.0000000      -0.0000000      -0.0000000
      -0.0000000      -0.0000000      -0.0000000      -0.0000000      -0.0000000
      -0.0000000      -0.0000000      -0.0000000      -0.0000000      -0.0000000
      -0.0000000      -0.0000000      -0.0000000      -0.0000000      -0.0000000
"""

def read_perpot(file):
    data = np.loadtxt(file).flatten()
    print(data.shape)
    nper = int(data[0])
    drp  = data[1]
    mmax = int(data[2])
    nth  = int(data[3])
    nc   = int(data[4])
    ampc = data[5:15]
    amps = data[15:25]
    print(nper,drp,mmax,nth,nc)
    print(ampc)
    print(amps)
    d3 = data[25:]
    print(d3.shape)
    d3 = d3.reshape(mmax,2,nper)
    print(d3.shape)
    cp = d3[:,0,:]
    sp = d3[:,1,:]
    print(cp.shape)
    for i in range(10):
        print('c',i,cp[i][:4],cp[i][:].max())
        print('s',i,sp[i][:4],sp[i][:].max())


read_perpot('perpot')
