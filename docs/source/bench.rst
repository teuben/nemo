.. _bench:

Benchmarking
============

How fast does the N-body treecode run?
To what degree does optimization/vectorizing help? When do
programs become I/O dominated? Some of the numbers quoted below should
be taken with great care, since a lot of other factors can go into
the timing result. 

A number of programs in NEMO have a command line parameter such as
``nmodel=N``, ``nbench=N`` or ``iter=N`` (N normally set to 1)
but together with ``help=c`` or prefixing with ``/usr/bin/time`` will
give an accurate measurement how long
the code takes to execute ``N`` loops of a particular algorithm. For
some programas their respective man pages discuss a particular benchmark.

On the top level we have ``make bench5`` and ``make bench``, the latter
dynamically controlled with the scattered ``Benchfile``'s


N-body integration
------------------

The standard NEMO benchmark of the treecode integration is to
``hackcode1``
without any parameters.  It will generate a spherical
stellar system in virial equilibrium with 128 particles, and
integrate it for 64 timesteps (``tol=1 eps=0.05``).   
In Table~\ref{t:bench1} the  amount
of CPU (in seconds) needed for **one** timestep is listed
in column 2. When not otherwise mentioned,
the code used is the standard NEMO ``hackcode1`` with
default compilation on the machine quoted. Note that one can often
obtain significant performance increase (factor of 2 on some sparc
architectures) by studying the native compiler and
in particular its optimization options.


.. Treecode Benchmarks}

.. list-table::    Treecode Benchmarks
   :header-rows: 1

   * - Machine
     - cpu-sec/step
     - code      
     - comments

   * - i5-1135G6 @ 4.2 GHz
     - 0.000089
     - hackcode1
     - 2020 laptop

   * - i7-8550U @ 4 GHz
     - 0.000178
     - hackcode1  -
     - 2018 laptop

   * - Sun-3/60
     - 5.400
     -
     - -fswitch (orig development)
   * - 3b1 (10Mhz 68010)
     - 49.000
     -
     -
   * - 386SX (16Mhz)
     - 87.000
     -
     - (linux) software floating point
   * - core 2 duo @ 2.0 GHz
     - 0.0012
     - hackcode1
     - 2007 laptop
     
     
i7-3630QM @ 3.4 GHz          & 0.000177 & hackcode1 & 2014 laptop \\
i70-870 @ 2.93 GHz	     & 0.00030 & hackcode1 & 2010 desktop \\

Dec-alpha		     & 0.0042 & hackcode1 & -O4 -fast \\
Dec-alpha		     & 0.0048 & hackcode1 & default \\
CRAY X/YMP48                 & 0.0060 & TREECODE V3 & estimate (1989) \\
Onyx-2			     & 0.0088 & hackcode1 & default (1996) \\
ETA-10                       & 0.010 & TREECODE V2 & estimate (1987)  \\
Sun Ultra-140		     & 0.012 & hackcode1 & -xO4 -xcg92 -dalign -xlibmil \\
Sun 20/62                    & 0.013 & hackcode1 & default (1994) \\
Cyber 205                    & 0.018 & TREECODE V2 & estimate (1986) \\
Sun 20/61                    & 0.020 & hackcode1 & \\
HP/UX 700                    & 0.020 & hackcode1 &  \\
Sun Ultra-140		     & 0.024 & hackcode1 & default \\
Sun 20/??		     & 0.024 & hackcode1 & -xO4 -xcg92 -dalign -xlibmil \\
G3 PowerPC 250Mhz	     & 0.026 & hackcode1 & -O \\
Sun 10/51                    & 0.029 & hackcode1 & -O -fast -fsingle \\
Cray-2                       & 0.029 & TREECODE2   & REAL - Pitt, oct 91\\
% SGI ???                      & 0.030 & hackcode1   & John Wangs machine
DEC DS3000/400 alpha         & 0.036 & hackcode1   & default compilation \\
Pentium-100                  & 0.038 & hackcode1   & default \\
SGI Indigo		     & 0.045 & hackcode1   & default compilation \\
CRAY YMP                     & 0.059 & hackcode1   & default compilation \\
% bootes:
Sparc-10                     & 0.063 & hackcode1   & using {\tt acc -cg92} \\
486DX4-100 (linux)           & 0.068 & hackcode1   & default \\
486DX2-66 (linux)            & 0.093 & hackcode1   & -DSINGLEPREC \\
Sparc-2	                     & 0.099 & gravsim V1  & \\
IBM R/6000                   & 0.109 & hackcode1   & default cc compiler \\
Dec 5000/200		     & 0.116 & hackcode1   & \\
Sparc-2                      & 0.130 & hackcode1   & -DSINGLEPREC -fsingle \\
Sparc-2                      & 0.180 & hackcode1   & \\
Multiflow 14/300             & 0.190 & hackcode1   & \\
Convex C220                  & 0.290 & & \\
NeXT                         & 0.240 &             & [ganymede 68040, nov 91]\\
Sparcstation1+               & 0.340 & & \\
Sun-4/60 Sparcstation 1      & 0.420 & & \\
Alliant FX??                 & 0.430 & gravsim V1 & \\
Alliant FX4/w 3 proc's       & 0.590 & & \\
VAX workstation 3500         & 0.970 & & \\
Sun-4/60 Sparcstation 1      & 1.040 & treecode2   & cf. C-code @ 0.420 \\
Sun-3/110                    & 1.660 & hackcode1 & fpa.il \\
Sun-3/60                     & 2.280 & & \\




The ``gravsim``
code\footnote{C version code of the treecode
written by Mark Bellon - Urbana, IL} is better suited for a
multiprocessor machine.  Its user interface and
database format are however different from NEMO's and interface scripts
can be defined which make working with this code a little easier.  Both
these C versions of the treecode ({\tt hackcode1} and {\tt gravsim}) are
inherently slower because they are recursive and spend most of their CPU
time in treewalking (with a lot of integer arithmetic).  
The modified (vectorized) Hernquist
fortran code (referred to as {\tt TREECODE V3}) \index{TREECODE}
has an approximate speedup of about a factor 200-400 over
the original VAX/Sun-3 speed on a CRAY supercomputer.

It is also perhaps interesting to quote that replacing the
{\tt sqrt} function by a very fast machine dependant one
will increase the speed of the C version of the treecode by
about 20\%. Some recent HP\index{HP} computers have a special
hardware floating point operation to perform {\tt 1/sqrt()}.

Nbody0
~~~~~~

The program is Aarseth's simplest
nbody code (contained in Binney and Tremaine, 1987, no regularization or nearest neighbors).
The input is
a Hubble expanding cartesian lattice, w/ 925 pts, GMtot=1, expansion
factor = 6 (omega = 1.2).  Long version followed for 60 time units,
short version for 5. Results are summarized in table below. First
table compiled by D. Richstone.

It seems the input data have been lost.


.. list-table::    N-Body0 Benchmarks
   :header-rows: 1

   * - Machine
     - time1 (sec)
     - time2 (sec)
     - speed
   * - Sun-4/110(Pele - 8Mb)
     - 
     - 21,753
     - 0.41
   * - Vaxstation 3100(Miffy - M48, 24Meg)
     -
     - 1302
     - 0.65

Sparc 1			&		& 1023	&	0.83 \\

Sparc IPC(Courage - 16 Mb) &	9,015	&	  850	&	1.000 \\

Sparc 2 	&	4,483	&		&	2.01 \\

Sparc 2'	&		&	  417	&	2.04 \\

Dec 5000/200 	&		&	  318	&	2.67 \\

Stardent(ism)	&		&	  211 	&	4.03 \\

IBM Risc (Juno)	&	2,117	&	  198	&	4.27 \\
					
IBM Risc (wibm01)&	2,115	&		&	4.26 \\
Convex		&		&	  172 	&	4.94 \\
HP/UX 700     &                  &     26.2   &     \\

Cray YMP	&		&	  19.1	&	44.5


Orbit integration
-----------------

Benchmark is taking 100,000 leapfrog steps. For 2D optimized 
potentials the timing on
a Sparc-1 station is about 12" for ``log`` or ``plummer``, and 
23" for ``teusan85`` in the core region (orbit remaining within
the body of the bar).


