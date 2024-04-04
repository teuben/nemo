Images
------


Initializing Images
~~~~~~~~~~~~~~~~~~~

There are a few programs with which images can be initialized:

- `ccdgen <https://teuben.github.io/nemo/man_html/ccdgen.1.html>`_
  is the simplest program to generate images from scratch, much like 
  `tabgen <https://teuben.github.io/nemo/man_html/tabgen.1.html>`_
  for tables. Here is an example creating and visualizing an image

.. code-block:: bash

  % ccdgen bar.ccd  object=bar spar=1,10,0.5,30 size=128,128
  % ccdplot bar.ccd
  % ccdfits bar.ccd  bar.fits
  % ds9 bar.fits


- **ccdmath** is the most straightforward program.  Here is an example
  of creating an image from scratch:


.. code-block:: bash

  % ccdmath out=ccd1 fie=%x+%y size=2,4
  Generating a map from scratch

  % tsf ccd1
  set Image
    set Parameters
      int Nx 2
      int Ny 4
      int Nz 1
      double Xmin 0.00000
      double Ymin 0.00000
      double Zmin 0.00000
      double Dx 1.00000
      double Dy 1.00000
      double Dz 1.00000
      double Xrefpix 0.00000 
      double Yrefpix 0.00000 
      double Zrefpix 0.00000 
      double MapMin -4.00000
      double MapMax 0.00000
      int BeamType 0
      double Beamx 0.00000
      double Beamy 0.00000
      double Beamz 0.00000
      double Restfreq 0.00000 
      double VLSR 0.00000 
      double Time 0.00000
      char Storage[5] "CDef"
      int Axis 1
    tes
    set Map
      double MapValues[2][4] 0.00000 1.00000 2.00000 3.00000 1.00000 
         2.00000 3.00000 4.00000 
    tes
  tes

  % ccdprint ccd1 x= y= label=x,y
   Y\X 0 1
 
  3  3 4 
  2  2 3 
  1  1 2 
  0  0 1 

- **snapgrid** converts a snapshot to an image.

- **fitsccd** converts a FITS file to an image. The inverse of this,
  **ccdfits** also exists. Storage model of FITS vs. NEMO's image
  format is currently not optimal, and thus for very large images
  a CPU penalty is incurred.


.. code-block:: bash

  nx,ny	->    data[nx][ny]

  e.g.	ccdmath out=ccd1   nx=10 ny=5
  	gives       double MapValues[10][5]  


  ccdmath "" - %x 3,2 | tsf - margin=100 | grep MapVal
	MapValues[3][2] -2.00000 -2.00000 -1.00000 -1.00000 0.00000 0.00000
  ccdmath "" - %y 3,2 | tsf - margin=100 | grep MapVal
	MapValues[3][2] -2.00000 -1.00000 -2.00000 -1.00000 -2.00000 -1.00000

	

Galactic Velocity Fields
~~~~~~~~~~~~~~~~~~~~~~~~

As an example, a
special section is devoted here to the analysis of 
galactic
velocity fields.\footnote{In this example
shell variables such as ``r=$(nemoinp 0:60)`` have been
replaced with the more portable macro files like
``@tmp.r``. Although the example uses ``0:60`` and works
fine in the shell the example was used under, increasing the
number to 256 would fail because of overflowing the maximum
characters allowed on the commandline}

The following programs are available:

.. code-block:: bash

	ccdvel          create a model velocity field, from scratch
	rotcur          tilted ring model velocity field fitting
	rotcurshape     annulus rotation curve shape fitting to a velocity field
	ccdmath         perform math on images, or use math to create images
	ccdplot         plot (contour/greyscale) an image
	ccdprint	print out pixel values in an imamge


    % nemoinp 0:60 > tmp.r
    % tabmath tmp.r - "100*%1/(20+%1)" all > tmp.v
    % ccdvel out=map1.vel rad=@tmp.r vrot=@tmp.v pa=30 inc=60
    % rotcurshape in=map1.vel radii=0,60 pa=30 inc=60 vsys=0 units=arcsec,1 \
                  rotcur1=core1,100,20,1,1 tab=-
 

    % ccdmath out=map0.vel fie=0 size=128,128
    % rotcurshape map0.vel 0,40 30 45 0 blank=-999 resid=map2.vel \
                  rotcur1=plummer,200,10,0,0 fixed=all units=arcsec,1



Since rotcurshape computes a residual velocity field, one can easily
create nice model velocity fields from any selected shape by 
*fitting* a rotation curve shape to a velocity field of all 0s
and keeping all parameters fixed to the requested values:

.. code-block:: bash

   % ccdmath out=map0.vel fie=0 size=128,128
   % rotcurshape map0.vel 0,40 30 45 0 blank=-999 resid=map.vel \
              rotcur1=plummer,200,10,0,0 fixed=all units=arcsec,1
   % ccdplot map.vel -100:100:10 blankval=0 cmode=1


..  rcshape1.ps


