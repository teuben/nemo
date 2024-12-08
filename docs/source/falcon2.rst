FALCON2
-------


The falcon V2 package is still under development for a public release. The main
integrator is now called ``griffin``, replacing the older ``gyrfalcON`` from falcon V1.
Although falcon V1 is distributed with NEMO, it is not maintained and may eventually
succumb to software rot. Users should switch to using ``griffin``.

A paper describing the ``griffin`` FMM code is in Dehnen (2014)

.. warning::
   Adding FALCON2 to NEMO will result in some programs that have duplicated names, e.g. **mkplummer**.


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

  s2a p100.in  | tabcomment - - delete=t | tabtos - p100.bsf block1=m,pos,vel,skip nbody=100
  tsf p100.bsf


- Comparing falcon1 with falcon2

.. code-block::

   #  some parameters
   nbody=10
   eps=0.05
   kmax=8
   tstop=1
   #
   rm p1.*
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

- Comparing X and Y


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

- mkplummer has a default time=1

- ..
