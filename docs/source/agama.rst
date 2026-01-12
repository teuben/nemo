AGAMA
-----

AGAMA (Action-based GAlaxy Modelling Architecture) is a software
library intended for a broad range of tasks within the field of
stellar dynamics.  The library contains a powerful framework for
dealing with arbitrary density/potential profiles and distribution
functions (analytic, extracted from N-body models, or fitted to the
data), a vast collection of general-purpose mathematical routines, and
covers many aspects of galaxy dynamics up to the very high-level
interface for constructing self-consistent galaxy models.

AGAMS's NEMO plugin allows to use any Agama potential as an external
potential in gyrfalcON and other NEMO programs (in a similar context
as the AMUSE plugin). The potential may be specified either as a file
with coefficients (for potential expansions), or more generally, as an
INI file with parameters of possibly several components defined in
groups [Potential1], [Potential whatever], . . .


Installation
~~~~~~~~~~~~

AGAMA can usually be installed as follows (as a user):

.. code-block::

   pip install agama

As a developer, it may be easier to install the source code via github, and follow the
following recipe somewhere in your workflow (skipping details where your python/virtual environment
is)

.. code-block::

   git clone https://github.com/GalacticDynamics-Oxford/Agama
   pip install -e Agama

note there are some requirements (e.g. GSL, Eigen, ...) but by default they will be
installed within Agama's source tree. See also $NEMO/usr/agama for installing and testing
AGAMA in a NEMO environment.

Example
~~~~~~~

This example is reproduced from the AGAMA manual. It considers the
evolution of an N-body model of a Plummer sphere of mass m = 0.1 and
radius a = 0.1 (satellite), which moves on a nearly-circular orbit
inside a larger Plummer sphere of mass M = 1 and radius A = 1 (host
galaxy), represented by a smooth external potential. First create the
satellite, after which the model is shifted over into a near-circular
orbit:

.. code-block::
   
      mkplum out=sat.nemo nbody=10000 r_s=0.1 mass=0.1
      snapshift sat.nemo sat_shift.nemo rshift=-1,0,0 vshift=0,-0.6,0      

Then create an INI file (``pot.ini``) describing the external potential:

.. code-block::
   
   [Potential]
   type=Plummer
   mass=1
   scaleRadius=1

Finally run the simulation for a few orbital periods:


.. code-block::
   
    gyrfalcON sat_shift.nemo sat_out1.nemo eps=0.02 kmax=5 step=0.25 tstop=50 accname=agama accfile=pot.ini


this should take about 20 seconds on a typical year-2020 CPU.
The satellite nicely orbits around the origin, leaving a small tidal tail.

Although this example is using the falcon tools (mkplum, gyrfalcON, agama),
using original NEMO tools a simular simulation could be achieved using
mkplummer, snapscale and hackcode3. We leave this for another version of
this manual.
