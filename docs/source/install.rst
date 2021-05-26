.. _install:

Installation
============

This is covered at length elsewhere, but to summarize here is an example:

.. code-block:: bash

   git clone https://github.com/teuben/nemo
   cd nemo
   ./configure --with-yapp=pgplot
   make build check bench 
   source nemo_start.sh


install_nemo.sh
---------------

The script ``$NEMO/docs/install.sh`` is an example installing NEMO
with a number of extra libraries, and python with a number of
modules. 


   

Advanced Installation
---------------------

It's a fact of life that you will not be satisified with the compiler
or libraries that your system provides. Add to this that you don't
have admin privilages, and you might be in for a rude awakening.

No worries, NEMO has you covered (to some degree).  We provide an
environment (a poor man's container) where most open source libraries
can be installed with a supported ``$NEMO/opt`` prefix. This means you
can configure packages using

.. code-block::

      --with-prefix=$NEMO/opt


of for *cmake* based packages

.. code-block::

      -DCMAKE_INSTALL_PREFIX=$NEMO/opt

For some packages this has been automated using the ``mknemo`` command, described in
the next section.

mknemo
------

Using the ``mknemo`` command a number of supported libraries that
are often used can be compiled within ``$NEMO/opt``, as described
in the previous section. The supporting scripts are generally
located ``$NEMO/src/scripts/mknemo.d``.

Examples:

.. code-block::

   mknemo cfitsio fftw gsl hdf4 hdf5 hypre netcdf4 wcslib

