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

A related package, ZENO, was spun off NEMO V1, and is maintained by Barnes.

	 NEMO:      ascl:1010.051
	 ZENO:      ascl:1102.027 (see also $NEMO/usr/zeno)

Packages we optionally use:

	 PGPLOT:    ascl:1103.002 (can be included with NEMO)
	 CFITSIO:   ascl:1010.001
	 WCSLIB:    ascl:1108.003
	 HDF:	    ascl:1502.009
	 glnemo2:   ascl:1110.008
	 gyrfalcON: ascl:1402.031 (included with NEMO)
	 gsl
	 plplot
	 unsio
	 uns_project

See README.install for installation guidelines. In it's simplest the following commands may work
(replace .csh with .sh if appropriate)

	 wget https://teuben.github.io/nemo/install_nemo
	 chmod +x install_nemo
	 ./install_nemo nemo=$HOME/opt/nemo
	 source $HOME/opt/nemo/nemo_start.csh

Some obvious and perhaps not so obvious tools you will need to have pre-installed:  C/C++/Fortran compiler,
csh, git, cmake. Use your local package manager to install those before you attempt to run the install_nemo
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



