#! /usr/bin/env python
#

import numpy

t = numpy.linspace(0,50,201) # grid in time
w = 0.6 # rotation frequency ω = V /R
x,y,z = numpy.cos(w*t), numpy.sin(w*t), t*0 # orbit coordinates
numpy.savetxt("center.txt", numpy.column_stack((t, x, y, z)))
ax,ay,az = -w**2*x, -w**2*y, 0*z # components of centrifugal acceleration
numpy.savetxt("accel.txt", numpy.column_stack((t, ax, ay, az)))
