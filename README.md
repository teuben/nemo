NEMO is a toolbox for stellar dynamics, particle simulations, stellar orbits,
image processing and tabular data manipulation. Documentation is maintained
in the github pages, https://teuben.github.io/nemo

This is the 4th major release of NEMO,  and although data are compatible
with earlier releases, old source code may need to be tweaked a
bit to compile and link in the newer releases. Some compatibility with ZENO
is also advertised. A brief history of NEMO:

   * NEMO V1:	IAS release (Barnes, Hut & Teuben, 1987)
   * NEMO V2:	UMD release (Teuben, 1994)
   * NEMO V3:	UMD release (Teuben, 2001) in CVS, with autoconf support and
		hooks into manybody.org modules starlab and partiview
   * NEMO V4:   UMD/ESO release (2017) now maintained in github

A related package, ZENO, was spun off NEMO V1, and is still maintained by Josh Barnes. Two
other packages that geneologically came after NEMO are StarLab and AMUSE
(see also https://ascl.net for code references):

	 NEMO:      ascl:1010.051
	 ZENO:      ascl:1102.027 (normally installed in $NEMO/usr/zeno)
	 STARLAB:   ascl:1010.076
	 AMUSE:     ascl:1107.007

Packages we optionally use (sometimes installed in $NEMO/opt via code in $NEMO/local):

	 PGPLOT:    ascl:1103.002
	 CFITSIO:   ascl:1010.001
	 WCSLIB:    ascl:1108.003
	 glnemo2:   ascl:1110.008
	 gyrfalcON: ascl:1402.031 (included with NEMO)
	 HDF4
	 HDF5       https://www.hdfgroup.org
	 netcdf4
	 gsl
	 plplot
	 unsio
	 uns_project

Tools you will need to have pre-installed: A C/C++/Fortran
compiler, (t)csh, and git.  For graphics it's probably
useful to have pgplot, but the default ps driver works
fine just to get started quickly.



There are a few ways to install NEMO.  Although there are
some install scripts with many options, and there is the README.install file
for background information, here is the basic method
for most Linux distros (assuming you have the preconditions):

         git clone https://github.com/teuben/nemo
         cd nemo
         ./configure --with-yapp=pgplot
         make build check bench 
         source nemo_start.sh

On the most recent apple controlled hardware, with SIP enabled, you're in for a rude
awakening. I use brew, and assuming you have gcc-10 (and related) and pgplot installed, this should
work (there are other ways to install tools on a mac,but don't get me started):

         git clone https://github.com/teuben/nemo
         cd nemo
         CC=gcc-10 CXX=g++-10 F77=gfortran-10 ./configure --disable-shared --with-yapp=pgplot
         make build check bench
         source nemo_start.sh



Once NEMO has been installed, here are some examples of scripts and
figures: https://teuben.github.io/nemo/examples/ or look at an example
ipython notebook
https://github.com/teuben/nemo/blob/master/nemo_start_example.ipynb
for something completely different.



