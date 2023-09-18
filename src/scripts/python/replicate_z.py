#! /usr/bin/env python
#
#  Replicate a FITS cube in the 3rd dimensions
#
#  Peter Teuben - 16-nov-2022 - Created

import os
import sys
import numpy as np
from astropy.io import fits
from nemopy import getparam


keyval = [
    "in=\n          input fits cube",
    "out=\n         output fits cube",
    "n=2\n          replication factor",
    "VERSION=0.1\n  16-nov-2022 PJT",
]

usage = """
replicate a cube in the 3rd dimension n times.


"""

p = getparam.Param(keyval,usage)


fin  = p.get("in")
fout = p.get("out")
n = int(p.get("n"))


hdu = fits.open(fin)
naxis3 = hdu[0].header['NAXIS3']
data = hdu[0].data
shape = data.shape

print('old shape: ',shape)

shape2 = (1,naxis3*n,shape[2],shape[3])
print('new shape:', shape2)
data2 = np.zeros(shape2)
print(data2.shape)
zoff = 0
for z in range(n):
    data2[0,zoff:zoff+naxis3,:,:] = data
    zoff = zoff + naxis3

hdu[0].header['NAXIS3'] = n*naxis3
hdu[0].data = data2
hdu.writeto(fout)
