#! /usr/bin/env python
#
#

import os
import sys
import aplpy
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


help_main = ["Simple color plot of a FITS image",
             "with options to pick a plane (if a cube) or slicing method",
             "colormaps and plot file extension can also be changed",
             "plot file name is derived from input FITS file",
             ]

help_color =  ['Popular colors: viridis, gist_heat gist_ncar (default)',
               '                rainbow, jet, nipy_spectral', 
               'https://matplotlib.org/stable/tutorials/colors/colormaps.html']


parser = argparse.ArgumentParser(description="\n".join(help_main),
                                 formatter_class=argparse.RawTextHelpFormatter)
                                 # formatter_class=argparse.ArgumentDefaultsHelpFormatter)

               

parser.add_argument('fitsfile',    help="input FITS file",        default=None)
parser.add_argument('--plane',     help="plane (if cube) [-1]",   default=-1,            type=int)
parser.add_argument('--pvar',      help="plane var (x,y,[z])",    default='z')
parser.add_argument('--color',     help="\n".join(help_color),    default='gist_ncar')
parser.add_argument('--ext',       help="plot type ([png],pdf)",  default='png')
parser.add_argument('--hist',      help="add histogram",          action="store_true")
parser.add_argument('--size',      help="plot size (inch)",       default=8,             type=float)

args  = parser.parse_args()


fitsfile = args.fitsfile
plane    = args.plane
color    = args.color
ext      = args.ext
pvar     = args.pvar
size     = args.size

if pvar == 'z':
    dims = [0,1]   # ra-dec
elif pvar == 'y':
    dims = [0,2]   # ra-vel
elif pvar == 'x':
    dims = [1,2]   # dec-vel
else:
    dims = [0,1]

hdu = fits.open(fitsfile)
if plane < 0:
    data = hdu[0].data
else:
    data = hdu[0].data[plane]

data = data[~np.isnan(data)]
data = data[data != 0]

dmin = np.min(data)
dmax = np.max(data)
dmean = np.mean(data)
dstd  = np.std(data)
print("Data min/max/mean/sig: %g %g %g %g" % (dmin,dmax,dmean,dstd))
dmin = dmean - 3*dstd;
dmax = dmean + 3*dstd;
print("Data min/max: %g %g" % (dmin,dmax))
bins = np.linspace(dmin, dmax, 32)
if args.hist:
    print("BINS: ",bins)

if args.hist:
    # side by side
    f=0.7
    #box1 = [0.1,0.1,0.8,0.8]   # full size, image
    box2 = [0.05,0.1,0.5/f,0.5/f]   # left side, image
    box3 = [0.55/f,0.1,0.2,f]  # right side, histo
else:
    # just a square image
    box1 = [0.1,0.1,0.8,0.8]   # full size, image
    #box2 = [0.1,0.1,0.5,0.5]   # left side, image
    #box3 = [0.7,0.15,0.2,0.4]  # right side, histo

try:
    if args.hist:
        fig = plt.figure(figsize=(size, f*size))
    else:
        fig = plt.figure(figsize=(size, size))
    
    if plane < 0:
        if args.hist:
            f1 = aplpy.FITSFigure(fitsfile, figure=fig, subplot=box2)
            ax_hist = fig.add_axes(box3)
            ax_hist.hist(data, bins=bins, orientation='horizontal', facecolor='blue',log=True)
        else:
            f1 = aplpy.FITSFigure(fitsfile, figure=fig, subplot=box1)
    else:
        if args.hist:
            f1 = aplpy.FITSFigure(fitsfile, slices=[plane], dimensions=dims, figure=fig, subplot=box2)
            ax_hist = fig.add_axes(box3)
            ax_hist.hist(data, bins=bins, orientation='horizontal', facecolor='blue',log=True)
        else:
            f1 = aplpy.FITSFigure(fitsfile, slices=[plane], dimensions=dims, figure=fig, subplot=box1)
except:
    print("problem processing %s in %s" % (fitsfile,os.getcwd()))
    sys.exit(0)
    
f1.show_grayscale()
f1.show_colorscale(cmap=color)
f1.add_colorbar()

try:
    f1.add_beam()
except:
    pass

# f.show_contour(fitsfile, levels=10)
f1.add_grid()
fig.canvas.draw()

idx = fitsfile.rfind('.fits')
if plane < 0:
    pfile = fitsfile[:idx] + ".%s" % ext
else:
    pfile = fitsfile[:idx] + ".%04d.%s" % (plane,ext)
# fig.subplots_adjust(right=0.15)   
fig.savefig(pfile)
print("Writing ",pfile)
# plt.show()
