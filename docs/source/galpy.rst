.. _galpy:

GALPY
-----

.. note::
   See https://docs.galpy.org/en/stable/getting_started.html


With some examples we show how things done in **galpy** can be done in **NEMO**


Rotation Curves
~~~~~~~~~~~~~~~

The following code example shows how to initialize a Miyamoto-Nagai disk potential and plot its rotation curve

.. tab:: galpy

   galpy prefers a normalization, which definese the value of the rotation curve at radius = 1

   .. code:: python

	     from galpy.potential import MiyamotoNagaiPotential
             #
	     mp= MiyamotoNagaiPotential(a=0.5,b=0.0375,normalize=1.)
	     mp.plotRotcurve(Rrange=[0.01,10.],grid=1001)


.. tab:: NEMO

   NEMO normally defines the mass as the first non-pattern speed parameter
   that normalizes the rotation curve

   .. code:: bash

      rotcurves miyamoto 0,1,0.5,0.0375 radii=0.01:10:0.01

and an example of a combined potential

.. tab:: galpy

   .. code:: python

	     from galpy.potential import NFWPotential, HernquistPotential
             #
	     mp= MiyamotoNagaiPotential(a=0.5,b=0.0375,normalize=.6)
	     np= NFWPotential(a=4.5,normalize=.35)
	     hp= HernquistPotential(a=0.6/8,normalize=0.05)
	     from galpy.potential import plotRotcurve
	     plotRotcurve(hp+mp+np,Rrange=[0.01,10.],grid=1001,yrange=[0.,1.2])
             #
	     mp.plotRotcurve(Rrange=[0.01,10.],grid=1001,overplot=True)
	     hp.plotRotcurve(Rrange=[0.01,10.],grid=1001,overplot=True)
	     np.plotRotcurve(Rrange=[0.01,10.],grid=1001,overplot=True)

.. tab:: NEMO

         galpy's rescale option is more complicated in NEMO. However, since speed scales with sqrt(mass),
         we can at least use the sqrt of galpy's normalization, and at least get the relative weights correct.

   .. code:: bash   

         rotcurves name1=miyamoto  pars1="0,sqrt(0.6),0.5,0.0375" \
	           name2=hernquist pars2="0,sqrt(0.06),0.6/8"     \
	           name3=nfw       pars3="0,1,1"                  \
	           radii=0.05:10:0.05


Orbit integration
~~~~~~~~~~~~~~~~~

.. tab:: galpy

   .. code:: python
      
	     from galpy.potential import MiyamotoNagaiPotential
	     from galpy.orbit import Orbit
	     #
	     mp= MiyamotoNagaiPotential(a=0.5,b=0.0375,amp=1.,normalize=1.)
	     o= Orbit([1.,0.1,1.1,0.,0.1])
             #
	     import numpy
	     ts= numpy.linspace(0,100,10000)
	     o.integrate(ts,mp,method='odeint')
	     #
	     o.plot()

.. tab:: NEMO

   .. code:: bash

	     # [R,vR,vT,z,vz] - not implemented
	     mkorbit r=1. vr=0.1 vt=1.1 z=0 vz=0.1   potname=miyamoto potpars=0,1,0.5,0.0375
	     orbint
	     orbplot



Surface of Section
~~~~~~~~~~~~~~~~~~

https://galaxiesbook.org/chapters/III-02.-Orbits-in-Triaxial-Mass-Distributions.html

.. tab:: galpy

   .. code:: python
      
	     from galpy.potential import LogarithmicHaloPotential
	     from galpy.orbit import Orbit
	     ts= numpy.linspace(0.,30,601)
	     lp= LogarithmicHaloPotential(normalize=True,b=0.9,core=0.2)
	     o= Orbit([0.1,0.,lp.vcirc(0.1,phi=0.),0.])
	     o.integrate(ts,lp)
	     o.animate(staticPlot=True); # remove the ; to display the animation

.. tab:: NEMO

   .. code:: bash

             mkorbit  x=0.1 y=0 z=0 vx=0 vy=0       potname=log potpars=0,1,0.2,0.9


      phi = v^2/2 log(x^2 + y^2/b^2 + a^2)
	     v^2/2 = m/a
