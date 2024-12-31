FALCON2
-------

Tmem he falcon V2 package is still under development for a public release. The main
integrator is now called ``griffin``, replacing the older ``gyrfalcON`` from falcon V1.
Although falcon V1 is distributed with NEMO, it is not maintained and may eventually
succumb to software rot. Users should switch to using ``griffin``.

A paper describing the ``griffin`` FMM code is in Dehnen (2014)

.. warning::
   Adding FALCON2 to NEMO will result in some programs that have duplicated names, e.g. **mkplummer**.
   The **sphinx** program is another example of a python program to build documentation.


Installation
~~~~~~~~~~~~

Hopefully soon we can use 

.. code-block::

   mknemo falcON2

but unless you are a close collaborator, this will likely not be working.

This type of install would place the package in ``$NEMO/local/falcon2``, the
file ``INSTALL.md`` describes the installation procedure.



Command Line Interface
~~~~~~~~~~~~~~~~~~~~~~

Although the CLI will look familiar to NEMO users, there are some salient differences.
Some highlights:

- the CLI is a set of ``key=val``, like in NEMO. If given in order, the ``key=`` portion can be
  omitted, leaving just a series of values. This not recommended in scripts.
- ``--help`` describes the keywords and some help, much like the ``help=h`` option in NEMO
- ``help=1`` shows hidden CLI options (NEMO doesn't have hidden options)
- ``debug=`` 




Examples
~~~~~~~~


For details on a specific program, type

.. code-block::

   program --help


- Create a Plummer sphere with 100 particlces

.. code-block::

   mkplummer p100.in 100

- View the contents of this HDF5 dataset

.. code-block::

   dump p100.in

- Convert falcon2 HDF5 files to NEMO snapshot, and review like the ``dump`` program

.. code-block::

   s2a p100.in  | tabcomment - - delete=t | tabtos - p100.bsf block1=m,pos,vel,skip nbody=100
   tsf p100.bsf

- Comparing performance of griffin in parallel (oneTBB) mode:

.. code-block::

   rm -f p10k.*
   mkplummer p10k.in 10000 time=0
   /usr/bin/time griffin p10k.in p10k.out step=1 tstop=10 eps=0.05 tau=2^-4
   threads=0    98.41user 19.02system 0:24.23elapsed 484%CPU 
   threads=1    21.19user  0.79system 0:22.03elapsed  99%CPU     
   threads=2    24.06user  2.90system 0:15.05elapsed 179%CPU
   threads=4    33.37user  4.73system 0:13.77elapsed 276%CPU  40.70user 5.90system 0:16.47elapsed 282%CPU 
   threads=8    62.97user 11.25system 0:21.23elapsed 349%CPU
           12   82.25user 14.55system 0:24.61elapsed 393%CPU
           16   93.44user 16.78system 0:25.63elapsed 430%CPU
           20   97.90user 18.67system 0:23.95elapsed 486%CPU
   
   
- Comparing falcon1 with falcon2

.. code-block::

   #  some parameters
   nbody=10
   eps=0.05
   kmax=8
   tstop=1
   #
   rm -f p1.*
   mkplummer p1.in $nbody time=0
   griffin p1.in p1.out step=1 tstop=$tstop eps=$eps tau=2^-$kmax
   # 8.7sec

   s2a p1.in | tabcomment - - delete=t | tabtos - p1.bsf block1=m,pos,vel,skip nbody=$nbody
   gyrfalcON p1.bsf p1.out2  step=1 tstop=$tstop eps=$eps kmax=$kmax 
   # 0.6 sec

   # comparing initials (notice p1.out does not contain times=0)
   echo "=== Initial conditions ==="
   s2a p1.in times=0 | tabcols - 2:7
   snapprint p1.out2 times=0 format=%11.9e


   # comparing final 
   echo "=== Final snapshot ==="
   s2a p1.out times=$tstop | tabcols - 2:7
   snapprint p1.out2 times=$tstop format=%11.9e


but griffin is more accurate and such a direct comparison may not be fair. In particular, is tau=2^-$kmax even fair?


Programs
~~~~~~~~

By default there are currently 21 programs installed:


* a2s -- ascii to snapshot converter
* calc -- a simple calculator  *[name conflicts with a common linux program]*
* corerad -- find core radius & density following Casertano & Hut (1985) and McMillan, Hut & Makino (1990)
* dump -- dump a hdf5 snapshot file to stdout  *[name conflicts with a common linux program]*
* gravity -- add gravity to snapshot(s)
* griffin -- N-body code
* join -- join falcON snapshot files   *[name conflicts with a common linux program]*
* manipulate -- use manipulators on falcON snapshots
* mkdisc -- make a simple circum-stellar gas disc
* mkgrid -- construct (possibly perturbed) fcc packing
* mkparker -- set up a parker wind
* mkplummer -- construct plummer sphere *[name conflicts with a NEMO program]*
* mkpolytrope -- construct a polytropic gas sphere
* mksphere -- initial conditions from an equilibrium distribution function
* mkstar -- make a single or binary star
* s2a -- snapshot to ascii converter
* s2s -- copy and manipulate falcON snapshots
* setH -- adapt SPH smoothing lengths
* snapprop -- evaluates bodies function over snapshot, reports to stdout
* sphinx -- SPH code  *[name conflicts with a common linux program]*
* symmetrize -- symmetrizes snapshots; can also be used to reduce N


Surprises
~~~~~~~~~

to be resolved

- dump --help

  does not show the help we usually see, but seems to think --help is a file


- "griffin --help"

  says: "please provide 'out', 'tau', 'eps'"
  why complain, we didn't attempt to run.

- mkplummer has a keyword 'q-ran', but itsn't it better to use q_ran, since that's
  a more common one used in all falcon programs. keep it consistent.

- mkplummer has a default time=1  (e.g. mkdisc does the expected time=0)

- ..
