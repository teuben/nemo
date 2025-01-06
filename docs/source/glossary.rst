.. _glossary:

Glossary
--------

.. glossary::


  acceleration
      Dataformat used in the :term:`falcon` package. Related to the
      :term:`potential` in NEMO, except one can use multiple
      acceleration descriptors, normally separated by a semi-colon.

  bodytrans
      Dataformat (normally a C function) that is used to
      perform arbitrary operations on expression
      variables used in snapshot's. Example:   xvar=x

  ccd
      Synonymous for image; most programs in NEMO which handle
      images start or end with ``ccd``, e.g.
      ``ccdfits``, ``fitsccd``, ``ccdmath``. Most programs handle
      1D, 2D and 3D "ccd" images.

  ECSV
      (Enhanced Character Separated Values) a popular
      self-describing ASCII table format popularized by astropy.
      https://docs.astropy.org/en/stable/io/ascii/ecsv.html

  falcon
      A subpackage in NEMO that hosts the gyrfalcON code
      and many related tools.
      See ``$NEMO/usr/dehnen/falcON``.
  
  fie
      Most expressions that you give to
      :term:`program keyword` s are
      parsed by *nemofie* and eventually ``fie``. (Nomenclature
      borrowed from :term:`GIPSY`).  Example:   a='sqrt(pi)*G'

  FITS
      "Flexible Interchange Transport System", a
      standard dataformat used to interchange data between
      machines. Commonly used for images, but has since its inception
      be modified to store tables as well. See also :term:`SDFITS`

  FWHM
      (Full Width Half Max): the effective resolution of the
      beam if normally given in **FITS** keywords BMAJ,BMIN,BPA.  The
      term **resolution** is used interchangeably.  The FWHM of
      a gaussian beam is 2.355 :math:`\sigma`.
  

  GIPSY
      The Groningen Image Processing System. A number of concepts in NEMO, and
      a few routines, have been taken liberally from :term:`GIPSY`. See also
      https://www.astro.rug.nl/~gipsy/
  

  history
      Each NEMO dataset normally contains a data processing
      history in the form of of history items at the beginnging of the
      dataset. They are normally kept
      track of by the data processing programs, and can be displayed
      with the program ``hisf``.

  image
      Dataformat in NEMO, used to represent 2- and
      3-D images/data cubes. See also :term:`ccd`.

  MIRIAD
      Another astronomical data reduction package,
      from which we have borrowed some user interfaces
      which are plug-compatible with
      our command-line syntax.

  orbit
      Dataformat in NEMO used to represent a
      stellar orbit; most programs in NEMO which handle orbits start
      or end with *orb*, for example ``orbint``, an orbit integrator.

  pixel
      PIXEL/voxel: an area in 2D or 3D representing selected
      phase space coordinates.

  potential
      Dataformat in NEMO used to represent a
      potential; most programs in NEMO which handle potentials start
      or end with *pot*, for example ``potlist``.
      Related to the :term:`acceleration`

  program keyword
      Keywords that are defined by the
      program only. They can be seen by using the **help=** keyword
      (in itself being a :term:`system keyword`).

  review
      A small user interface that pops up when a
      program is interrupted. Type ``quit`` to exit it, or ``?``
      for help. This feature of the user interface may not be
      installed in your version.

  SDFITS
      A FITS extension that allows efficient storage of single dish
      spectra. See also :term:`FITS`

  set
      Compound hierarchical data-structure of a
      structured file. They are the equivalent of a C *struct*..
      A *set* is alwayd closed with a *tes*.

  snapshot
      Dataformat used in NEMO to represent an
      N-body system. Many programs that handle {\it snapshot}'s in
      NEMO start or end with *snap*, for example ``snapplot``.

  structured file
      The binary data NEMO writes is in a
      hierarchical structured format. Programs like
      `tsf  <https://teuben.github.io/nemo/man_html/tsf.1.html>`_
      `rsf  <https://teuben.github.io/nemo/man_html/rsf.1.html>`_,
      and 
      `csf  <https://teuben.github.io/nemo/man_html/csf.1.html>`_
      perform general and basic I/O functions on
      such files. They are hierarchical structured sets, much like
      how binary XML files would look.

  system keyword
      Global keyword that every NEMO
      program knows about, and are not listed in the :term:`program keyword` s
      that can be seen by issuing e.g. **help=** (in itself
      being a system keyword). This concept originated in :term:`GIPSY` where
      hidden keywords are also used.

  table
      A table consists of rows and columns of values, numbers or text.
      Most commonly stored in ASCII. Less well defined, it is one of the
      four data types commonly used in NEMO. Most programs that handle
      tables start or end with *tab*, for example ``tabplot``.
      

  voxel
      A three dimensional pixel. See also pixel.

  yapp
      "Yet Another Plotting Package", the library
      definition API that is used by all programs that produce graphics
      output. It is deliberately kept very simple to promote portability
      to lower level graphics packages. The **yapp=** system keyword
      controls the graphics device definitions/capabilities.  Examples
      are the ``ps`` and ``pgplot`` implementations.


