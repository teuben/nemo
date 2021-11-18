#! /usr/bin/env python
#
#

import os
import sys
import aplpy

#   dims for an [ra,dec,vel] cube
dims = [0,1]   # ra-dec
#dims = [1,2]   # dec-vel
#dims = [0,2]   # ra-vel

fitsfile = sys.argv[1]
if len(sys.argv) > 2:
    plane = int(sys.argv[2])
else:
    plane = -1
    
if plane < 0:
    f = aplpy.FITSFigure(fitsfile)
else:
    f = aplpy.FITSFigure(fitsfile, slices=[plane], dimensions=dims)
    
f.show_grayscale()
f.show_colorscale(cmap='gist_heat')
f.add_colorbar()

try:
    f.add_beam()
except:
    pass

# f.show_contour(fitsfile, levels=10)
f.add_grid()
f.save(fitsfile + ".pdf")
