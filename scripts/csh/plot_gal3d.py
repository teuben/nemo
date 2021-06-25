#! /usr/bin/env python
#
#  Plot the data created with mk_con3
#
#

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import aplpy

basename = sys.argv[1]

slicex = 64
slicey = 64
slicev = 32

fig = plt.figure(figsize=(11, 8))


# moments

gc1 = aplpy.FITSFigure(basename + ".den.fits", figure=fig, subplot=(3,3,1))
gc1.show_colorscale(cmap='gist_heat')
gc1.set_title("mom0")
gc1.axis_labels.hide()
gc1.tick_labels.hide()

gc2 = aplpy.FITSFigure(basename + ".vel.fits", figure=fig, subplot=(3,3,2))
gc2.show_colorscale(cmap='rainbow')
gc2.set_title("mom1")
gc2.axis_labels.hide()
gc2.tick_labels.hide()

gc3 = aplpy.FITSFigure(basename + ".sig.fits", figure=fig, subplot=(3,3,3))
gc3.show_colorscale(cmap='gist_heat')
gc3.show_colorscale()
gc3.set_title("mom2")
gc3.axis_labels.hide()
gc3.tick_labels.hide()

# cube slices

gc4 = aplpy.FITSFigure(basename + ".fits", figure=fig, subplot=(3,3,4), slices=[slicev])
gc4.show_colorscale(cmap='gist_heat')
gc4.show_colorscale()
gc4.set_title("XY slice")
gc4.axis_labels.hide()
gc4.tick_labels.hide()

gc5 = aplpy.FITSFigure(basename + ".fits", figure=fig, subplot=(3,3,5), slices=[slicey], dimensions=[0,2])
gc5.show_colorscale(cmap='gist_heat')
gc5.show_colorscale()
gc5.set_title("XV slice")
gc5.axis_labels.hide()
gc5.tick_labels.hide()

gc6 = aplpy.FITSFigure(basename + ".fits", figure=fig, subplot=(3,3,6), slices=[slicex], dimensions=[1,2])
gc6.show_colorscale(cmap='gist_heat')
gc6.show_colorscale()
gc6.set_title("YV slice")
gc6.axis_labels.hide()
gc6.tick_labels.hide()

# R-V

gc7 = aplpy.FITSFigure(basename + ".rvma.fits", figure=fig, subplot=(3,3,7))
gc7.show_colorscale(cmap='gist_heat')
gc7.show_colorscale()
gc7.set_title("R-V Major")
gc7.axis_labels.hide()
gc7.tick_labels.hide()

gc8 = aplpy.FITSFigure(basename + ".rvmi.fits", figure=fig, subplot=(3,3,8))
gc8.show_colorscale(cmap='gist_heat')
gc8.show_colorscale()
gc8.set_title("R-V Minor")
gc8.axis_labels.hide()
gc8.tick_labels.hide()


#gc.add_beam()

gc8.save('plot_gal3d.png')
gc8.save('plot_gal3d.pdf')

plt.show()
