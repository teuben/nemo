AGAMA
-----

AGAMA (All-purpose galaxy modelling architecture) 


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
   cd Agama

note there are some requirements (e.g. GSL, Eigen, ...)



---


Nemo [58] is a collection of programs for performing and analyzing
N-body simulations, which use common data exchange format and
UNIX-style pipeline approach to chain together several processing
steps.  The centerpiece of this framework is the N-body simulation
code gyrfalcON [23]. It computes the gravitational force between
particles using the fast multipole method, and can optionally include
an external potential.

The Nemo plugin allows to use any Agama potential as an external
potential in gyrfalcON and other Nemo programs (in a similar context
as the Amuse plugin). The potential may be specified either as a file
with coefficients (for potential expansions), or more generally, as an
INI file with parameters of possibly several components defined in
groups [Potential1], [Potential whatever], . . .


To build the plugin, one needs to have Nemo installed (obviously) and
the environment variable $NEMO defined; then make nemo will compile
the plugin and place it in $NEMO/obj/acc folder, where it can be found
by Nemo programs. For instance, this adds an extra potential in a
gyrfalcON simulation:

    $ gyrfalcON infile outfile accname=agama accfile=mypot.ini [accpars=1.0] ...

where the last optional argument specifies the pattern speed ¦
(frequency of rotation of the potential figure about z axis). All
units in the INI or coefs file here should follow the convention G
= 1.
