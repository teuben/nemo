FALCON2
-------


The falcon V2 package is still under development for a public release. The main
integrator is now called griffin, replacing the old gyrfalcON from falcon V1.

Installation
~~~~~~~~~~~~

Hopefully soon we can use 

.. code-block::

   mknemo falcON2

but unless you are a close collaborator, this will likely not be working.

This type of install would place the package in ``$NEMO/local/falcon2``.



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

- Convert falcon2 HDF5 files to NEMO snapshot, and recview junk like the dump program

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


but where timescales need to be checked

- Comparing X and Y


Surprises
~~~~~~~~~

- dump --help

  does not show the help we usually see, but seems to think --help is a file


- griffin does 
