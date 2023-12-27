.. _galpy:

GALPY
=====

.. note::
   See https://docs.galpy.org/en/stable/getting_started.html


With some examples we show how things done in **galpy** can be done in **NEMO**


Rotation Curves
---------------

The following code example shows how to initialize a Miyamoto-Nagai disk potential and plot its rotation curve

.. tab:: galpy

	 we will

   .. code:: python

	     from galpy.potential import MiyamotoNagaiPotential
             #
	     mp= MiyamotoNagaiPotential(a=0.5,b=0.0375,normalize=1.)
	     mp.plotRotcurve(Rrange=[0.01,10.],grid=1001)

.. tab:: NEMO

	 we will

   .. code:: bash

      rotcurves miyamoto 0,a,b,m
      cp a b


and now a combined rotation curve

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

         we will
