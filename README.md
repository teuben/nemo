NEMO is a toolbox for stellar dynamics, particle simulations, stellar orbits,
image processing and tabular data manipulation. See also https://teuben.github.io/nemo 

This is the 4th major release of NEMO,  and although data are compatible
with earlier releases, source code may need to be tweaked a
bit to compile and link in the newer releases. Some compatibility with ZENO
is also advertised.

   * NEMO V1:	IAS release (Barnes, Hut & Teuben, 1987)
   * NEMO V2:	UMD release (Teuben, 1994)
   * NEMO V3:	UMD release (Teuben, 2001) in CVS, with autoconf support and
		hooks into manybody.org modules starlab and partiview
   * NEMO V4:   UMD/ESO release (2017) now maintained in github

A related package, ZENO, was spun off NEMO V1, and is still maintained by Josh Barnes.

	 NEMO:      ascl:1010.051
	 ZENO:      ascl:1102.027 (see also $NEMO/usr/zeno)

Packages we optionally use:

	 PGPLOT:    ascl:1103.002
	 CFITSIO:   ascl:1010.001
	 WCSLIB:    ascl:1108.003
	 glnemo2:   ascl:1110.008
	 gyrfalcON: ascl:1402.031 (included with NEMO)
	 HDF4
	 HDF5       https://www.hdfgroup.org
	 gsl
	 plplot
	 unsio
	 uns_project

See README.install for installation guidelines. In its simplest form
the following commands may work (use the .sh or .csh version depending
on your $SHELL):

	 wget https://teuben.github.io/nemo/install_nemo
	 chmod +x install_nemo
	 ./install_nemo
	 source nemo/nemo_start.sh

Some obvious and perhaps not so obvious tools you will need to have pre-installed:  C/C++/Fortran compiler,
(t)csh, git, cmake. Use your local package manager to install those before you attempt to run the install_nemo
script. For example, on a fresh Ubuntu distro, you're likely going to need to install two
packages before running the install script, and triggering a few more using the ubuntu=1 flag, viz.

	 sudo apt install tcsh git
 	 ./install_nemo ubuntu=1


And you don't have wget (e.g. MacOSX) use curl or any other program (even any browser works) to get that script.
E.g. on a Mac where you have installed Homebrew with pgplot, this should get you a working NEMO with plotting
enabled:

	 curl -O https://teuben.github.io/nemo/install_nemo
	 chmod +x install_nemo
	 ./install_nemo brew=1
	 source nemo/nemo_start.sh

If you have used MacPorts, I may owe you another hint.

At the simplest level, the installation uses an autoconf based scheme:

         git clone https://github.com/teuben/nemo
         cd nemo
         ./configure ...
         make build
         source nemo_start.sh

Once NEMO has been installed, here are
some examples of scripts and figures: https://teuben.github.io/nemo/examples/
or look at an example ipython notebook https://github.com/teuben/nemo/blob/master/nemo_start_example.ipynb



