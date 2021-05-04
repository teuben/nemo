Glossary
--------

.. glossary::


  accelleration
      Dataformat 

  bodytrans
      Dataformat that is used to
      perform arbitrary operations on expression
      variables used in snapshot's.

  ccd
      Synonymous for image; most programs in NEMO which handle
      images start or end with ``ccd``.

  ECSV
      (Enhanced Character Separated Values) a popular
      self-describing ASCII table format popularized by astropy

  
  fie
      Most expressions that you give to program
      keywords are
      parsed by *nemofie* and eventually ``fie``. (Nomenclature
      borrowed from\GIPSY)

  FITS
      ``Flexible Interchange Transport System'', a
      standard dataformat used to interchange data between
      machines. Commonly used for images.

  FWHM
      (Full Width Half Max): the effective resolution of the
      beam if normally given in **FITS** keywords BMAJ,BMIN,BPA.  The
      term **resolution** is used interchangeably.
  

  history
      Each NEMO dataset normally contains a data
      history in the form of of history items at the beginnging of the
      dataset. They are normally kept
      track of by the data processing programs, and can be displayed
      with the program ``hisf``.

  image
      Dataformat in NEMO, used to represent 2- and
      3-D images/data cubes. See also :term:`ccd`.

  miriad
      Another astronomical data reduction package,
      from which we have borrowed some user interfaces
      which are plug-compatible with
      our command-line syntax.

  orbit
      Dataformat in NEMO used to represent a
      stellar orbit; most programs in NEMO which handle orbits start
      or end with *orb*.

  pixel
      PIXEL/voxel: an area in 2D or 3D representing

  potential
      Dataformat in NEMO used to represent a
      potential; most programs in NEMO which handle potentials start
      or end with *pot*. Related are the :term:`accelleration`

  program keyword
      Keywords that are defined by the
      program only. They can be seen by using the **help=** keyword
      (in itself being a system keyword).

  review
      A small user interface that pops up when a
      program is interrupted. Type ``quit`` to exit it, or ``?``
      for help. This feature of the user interface may not be
      installed in your version.

  set
      Compound hierarchical data-structure of a
      structured file. They are the equivalent of a C structure.

  snapshot
      Dataformat used in NEMO to represent an
      N-body system. Many programs that handle {\it snapshot}'s in
      NEMO start or end with *snap*.

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
      program knows about, and are not listed in the (program)
      keywords that can be seen by issuing e.g. **help=** (in itself
      being a system keyword).

  table
      A table consists of rows and columns of values, numbers or text.
      Most commonly stored in ASCII. Less well defined, it is one of the
      four data types in NEMO.

  yapp
      ``Yet Another Plotting Package'', the library
      definition that is used by all programs that produce graphics
      output. It is kept very simple. The **yapp=** system keyword
      controls the graphics device definitions/capabilities.


