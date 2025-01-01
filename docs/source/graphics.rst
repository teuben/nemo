.. _graphics:

Graphics and Image Display
==========================

.. note::
   Most NEMO graphics programs select their graphics output with
   the ``yapp=`` system keyword.

NEMO programs also need to display their data of course.
Here we will make a distinction between *graphics* and *image* data.
A simple but flexible *graphics* interface has been defined in NEMO and is used
extensively in programs that need such output.  These are controlled by
the user via the ``yapp=`` :ref:`system_keywords`.

To display *image* data we rely mostly (but see
`ccdplot(1NEMO) <https://teuben.github.io/nemo/man_html/ccdplot.1.html>`_
for an exception)
on external software.
Often images would need to be copied to a FITS file for this
(but see ``nds9`` for an example that can use NEMO's image format), and
use 3rd party programs such as *ds9*.

The YAPP graphics interface
---------------------------

The programs in NEMO which use *graphics* are rather simple and generally allow no
interactive processing, except perhaps for a simple 'hit-the-Enter-key'
or 'click-a-mouse-button' between successive plots or actions.  A very
simple interface (API) was defined (**yapp**, Yet Another Plotting Package)
with basic plot functions.  
There are currently a few yapp implementations
available, with a postscript-only device, and pgplot the most common ones.
If your output device is not supported by the ones available
in the current yapp directory
(``$NEMO/src/kernel/yapp``), you may have to write a new one!
A reasonably experienced programmer writes a functional yapp-interface in
an hour or so.

Although this method results in a flexible graphics interface, a
program can currently only be linked with one yapp-interface at a time, which
is selected at NEMO's installation time via the **configure** script.
This might
result in the existence of more than one version of the same
program, each for another graphics output device.  We use the 
convention that the ones for a
postscript printer have a ``_ps`` appended to their original name: the 
program which has the original name is the one whose display is the current
screen,
Hence we may see program names such as {\tt snapplot} (the default),
``snapplot_ps`` (postscript), or
``snapplot_cg`` (color Sun screen when this was popular) .
Again: actual names may differ on your system, and it may never be used.

If programs are linked with the multiplexing libraries
*yapp_pgplot* interface, several device drivers are transparently present through
pgplot, and the system keyword ``yapp=`` is then used to select
a device (a default can be set by using the **YAPP** environment
variable). 
See also the
:ref:`iface`
chapter or the
`yapp(5NEMO) <https://teuben.github.io/nemo/man_html/yapp.5.html>`_
manual page.

However, despite these grim sounding words, we currently
almost exclusively use the PGPLOT implementation of yapp, which is very flexible.
You will need to have pgplot installed on your system, which is another story
in itself.

.. note::
   The command ``mknemo pgplot`` may work if you cannot find a system package installer.
   This will install PGPLOT system files in ``$NEMOLIB``.  



pyplot=: python matplotlib
--------------------------

An experimental ``pyplot=`` keyword has been added to
``tabplot`` and ``tabhist``, which creates a simple python script
that reproduces (to some degree) what those program intended to do.
The intent is that these can be edited and made more functionality.
There are also pure python versions of tabplot and tabhist under
development.

General Graphics Display
------------------------

Another convenient way to present data in graphical form is by using
the table format. We have already encountered the *tables* created by
many NEMO programs. These tables can be used by NEMO programs
such as
`tabplot(1NEMO) <https://teuben.github.io/nemo/man_html/tabplot.1.html>`_
and 
`tabhist(1NEMO) <https://teuben.github.io/nemo/man_html/tabhist.1.html>`_, 
and other packages
such as
*gnuplot*,
*xgobi*,
*xmgrace*, 
*xgraphic*, 
*glueviz*, and
*topcat*.


Image Display Interface
-----------------------

Data in 
`image(5NEMO) <https://teuben.github.io/nemo/man_html/image.5.html>`_
format can be transferred in
`fits(5NEMO) <https://teuben.github.io/nemo/man_html/fits.5.html>`_
format and subsequently displayed and analyzed within
almost any astronomical image processing system.  They are generally much
better equipped to display and manipulate data of this kind of format. 
A number of standalone display programs can also understand FITS
format.  An excellent example of this is 
*ds9*, although it understands FITS files, can be used in
a client-server setting and NEMO image files can be directly sent
to the display server (a temporary fits file is created, which
can have drawbacks if they are large):

.. code-block::

    % ds9 &
    % nds9 map.ccd


Other programs you can consider are:  *carta*, *qfitsview*, *ginga*, and *fv*.
See also the list on https://blends.debian.org/astro/tasks/viewers
