.. _crosswalk:

Crosswalk
---------

Comparing functionality in NEMO with like-wise programs in other packages can be
helpful. Here are a few examples:  

.. list-table:: Comparing NEMO programs with those in other packages
   :header-rows: 1
      
   * - NEMO
     - python
     - CASA
     - MIRIAD
   * - tabtrend
     - numpy.diff()
     - N/A
     - N/A
   * - tabplot
     - plt.plot()
     - N/A
     - N/A
   * - fitsccd
     - hdu=astropy.io.fits.open()
     - importfits()
     - fits op=xyin
   * - ccdfits
     - hdu=astropy.io.fits.open()
     - exportfits()
     - fits op=xyout
   * - ccdprint
     - print(hdu[0].data)
     - N/A
     - imlist
   * - ccdmom
     - np.sum()
     - immoments()
     - moment
   * - ccdmath
     - numpy
     - immath()
     - maths
   * - ccdstat
     - hdu[0].data.std()
     - imstat()
     - imstat

An expanded example how to read a FITS (image) into the respective systems is as follows,
plus some sample commands one can in the respective package:


	
.. code-block::

  # NEMO
  $ fitsccd ngc0001.fits ngc0001.ccd
  $ ccdhead ngc0001.ccd
  $ ccdplot ngc0001.ccd yapp=/xs

  # PYTHON
  In [1]: hdu = astropy.io.fits.open('ngc0001.fits')

  # CASA
  CASA <1>: importfits('ngc0001.fits', 'ngc0001.im')
  CASA <2>: imhead('ngc0001.im')

  # MIRIAD
  $ fits in=ngc0001.fits out=ngc0001.mir op=xyin
  $ imhead in=ngc0001.mir
  $ implot in=ngc0001.mir device=/xs

   
  
