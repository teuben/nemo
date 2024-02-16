.. _install:

Installation
============

.. note::
   NEMO is available on https://github.com/teuben/nemo

Installation from github
------------------------

Installation is normally done by getting the source code via github.


The One Liners
~~~~~~~~~~~~~~

Here is an example, just 3 lines in
your (bash) shell, using a configurable helper script:

.. code-block:: bash

   wget https://teuben.github.io/nemo/install_nemo.sh
		
or

.. code-block:: bash
   
   curl -LO https://teuben.github.io/nemo/install_nemo.sh

after which installation can be done with something like

.. code-block:: bash
   
   bash install_nemo.sh nemo=$HOME/opt/nemo yapp=pgplot bench5=1 python=1
   source  $HOME/opt/nemo/nemo_start.sh

where the arguments to the
`install_nemo.sh <https://github.com/teuben/nemo/blob/master/docs/install_nemo.sh>`_
script are optional, but a few are
given to show some often use non-defaults. See that script for more details,
or use the `-h` flag

A more full example
~~~~~~~~~~~~~~~~~~~

A more manual install, bypassing this script, can be:

.. code-block:: bash

   git clone https://github.com/teuben/nemo
   cd nemo
   ./configure --with-yapp=pgplot
   make build check bench5 python
   source nemo_start.sh


On a Mac with 
`SIP protection <https://macpaw.com/how-to/disable-enable-system-integrity-protection>`_,
enabled, the ``--disable-shared`` flag may need to be added.

.. code-block:: bash

   git clone https://github.com/teuben/nemo
   cd nemo
   ./configure --with-yapp=pgplot --disable-shared
   make build check bench5
   source nemo_start.sh

Disabling SIP is not recommended, so we've been told.	On a Mac you will also need to have
Xcode installed, and gfortran (e.g. via brew). We need a special section on this

Pre-conditions
--------------

For a minimal install a number of packages need to be present on your system. Compilers, the make
utility, the csh shell, etc.   For some systems (e.g. Ubuntu) we keep a list of minimum
requirements of the packages that you will need for a minimal install.n

.. code-block:: bash

   cd $NEMO
   make install_apt
   cat src/scripts/requirements/apt.txt
   cat src/scripts/linux/ubuntu20.04

where the last ``ubuntu20.04`` file is a more complete list of packages.


Rebuilding
----------

If you have an existing installation, but many things have change, this is probably the preferred method:

.. code-block:: bash
   
   cd $NEMO
   git pull
   make rebuild

this will also preserve the possibly peculiar options for configure that you passed the first time it was installed.
Or more importantly, if you had edited the ``$NEMOLIB/makedefs`` file.

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
(anaconda, conda, venv etc.), we will not recommend a method, as this will
likely depend on your own situation. The installation examples below
should give you enough information how to adapt it for your python
installation.  

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


Package Managers
----------------

Most operating systems will have some package manager that controls how software
is installed. There is also a list in ``$NEMO/src/scripts/linux`` and ``$NEMO/src/scripts/brew``,
but here we list a few common ones:

.. tab:: Ubuntu

   The package manager is called ``apt``

   .. code:: bash

      sudo apt install ...

      build-essential	     
      gcc
      g++
      gfortran
	     
      pgplot5
      rman
      xorg-dev	    

.. tab:: Fedora

   The package manager is called ``dnf`` (formerly ``rpm``)

   .. code:: bash

      sudo dnf install ...
	     
      gcc
      gcc-gfortran
      gcc-g++
      tcsh
      make
      libtirpc-devel

      pgplot-devel
      cfitsio-devel
      netcdf-devel
      hdf-devel
      hdf5-devel


.. tab:: RedHat

   Not tested, probably same as Fedora

   .. code:: bash

      sudo dnf install ...

.. tab:: Homebrew

   The package manager is called ``brew``, but installation is done via https://brew.sh   

   Normally installed in the users own space. Prepend with the usual "sudo" if need be. Can be used
   on both Linux and Mac.  Recent versions have barred pgplot, because of licencing issues.

   .. code:: bash

      brew install ...


