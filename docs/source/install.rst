.. _install:

Installation
============

.. note::
   NEMO is available on https://github.com/teuben/nemo

Installation from github
------------------------

Installation is normally done via github. Here is a simple example, just 3 lines in
your (bash) shell, using a configurable helper script:

.. code-block:: bash

   wget https://teuben.github.io/nemo/install_nemo.sh
   bash install_nemo.sh  nemo=$HOME/opt/nemo yapp=pgplot bench5=1 python=1
   source  $HOME/opt/nemo/nemo_start.sh

where the arguments to the
`install_nemo.sh <https://github.com/teuben/nemo/blob/master/docs/install_nemo.sh>`_
script are optional, but a few are
given to show some often use non-defaults. See that script for more details.

A more manual install, bypassing this script, can be:

.. code-block:: bash

   git clone https://github.com/teuben/nemo
   cd nemo
   ./configure --with-yapp=pgplot
   make build check bench bench5 python
   source nemo_start.sh


On a Mac with their new
`SIP protection <https://macpaw.com/how-to/disable-enable-system-integrity-protection>`_,
the ``--disable-shared`` flag needs to be added

.. code-block:: bash

   git clone https://github.com/teuben/nemo
   cd nemo
   ./configure --with-yapp=pgplot --disable-shared
   make build check bench bench5
   source nemo_start.sh

Disabling SIP is not recommended.		


Rebuilding
----------

If you have an existing installation, but many things have change, this is probably the preferred method:

.. code-block:: bash
   
   cd $NEMO
   git pull
   make rebuild

this will also preserve the possibly peculiar options for configure that you passed the first time it was installed.

Advanced Installation
---------------------

It's a fact of life that you will not be satisified with the compiler
or libraries that your system provides. Add to this that if you don't
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

as NEMO generally adds the $NEMO/opt tree search for include and library files, as
well as adding its binaries to the search path.

For some packages this has been automated using the ``mknemo`` command, described in
the next section.

mknemo
------

Although the ``mknemo`` script was intended to quickly compile a NEMO program
(from any directory), and without the need to know where the source code lives.
It is now also used to aid the installation
of a number of supported libraries that
can be used by NEMO. They are compiled within ``$NEMO/local``, and will be installed
in ``$NEMO/opt``, as described
in the previous section. The supporting scripts are generally
located ``$NEMO/src/scripts/mknemo.d`` for you to examine.

Examples:

.. code-block::

   mknemo cfitsio fftw gsl hdf4 hdf5 hypre netcdf4 wcslib


The :ref:`progr` will give some advanced examples how to
deal with other libraries, and writing your own programs
or one of the plugins.

python
------

With so many useful python packages around, and so many different methods
(anaconda, cond, venv etc.), we will not recommend a method, as this will
likely depend on your own situation. The installation examples below
should give you enough information how to adapt it for your python
installation.  It goes without saying (this is 2021) we only support
python3.

However, if you install python from within NEMO, there will be a
``$NEMO/anaconda3`` directory, that gets automatically activated once
NEMO is loaded. Here is how you can install that version:

.. code-block::

      cd $NEMO
      make python

This will install a few python modules we often wind up using:
**amuse-framework**,
**amuse-galactics**,
**amuse-gadget2**,
**amuse-bhtree**,
**astromartini**,
**gala**,
**galpy**,
**pynbody**,
**python-unsio**,
**python-unsiotools**,
and
**yt**

For a number of these we have small test scripts to see if they are functional:

.. code-block::

      cd $NEMO/src/scripts/python
      make tests
   

For the cases where you want some control and be in developer mode, we
suggest the recommended practice of placing the code in ``$NEMO/local``,
as can be seen in the example below


.. code-block::

      cd $NEMO/local
      git clone https://github.com/webbjj/clustertools
      pip install -e clustertools


For a few packages, we have a few existing examples in the ``$NEMO/usr`` tree
(e.g. amuse, martini, unsio and uns_projects)
