AMUSE and NEMO
==============

AMUSE (Astrophysical Multipurpose Software Environment) originates some ideas
from its predecessors: ACS, StarLab and NEMO, but uses the python language.
Another feature of AMUSE is that
python is also the *glue* shell between legacy codes that can orchestrate
simulations taking components from different codes, whereas in NEMO legacy codes
have a NEMO CLI interface, at best.

For seasoned
`AMUSE <https://amuse.readthedocs.io/en/latest>`_
users, here we highlight some differences between the two, and give some examples
how to achieve the same task in NEMO and AMUSE.


Differences
-----------

- **Shell**:
  NEMO uses a Unix shell, AMUSE uses python (ipython, jupyter, ...). 

- **Units**:
  NEMO uses dimensionless values, and units are implied. 
  Most programs actually use virial units where G=1, but there are a few programs
  (e.g. galaxy, nbodyX) that use other units. The
  `units(1NEMO)  <https://teuben.github.io/nemo/man_html/units.1.html>`_
  tries to help you converting.
  AMUSE (optionally?) attaches units to numbers , using a python trick, e.g.

.. code-block::

   from amuse.units import units

   mass   = 1.0 | units.MSun

**astropy** users might be a bit baffled, since this looks very different. But

.. code-block::

   m1 = mass.as_astropy_quantity() 

will look more familiar.   In pure **astropy** it might look as follows:

.. code-block::

   from astropy import units as u

   m = 1.0 * u.solMass
   m2 = m.to(u.kg).value

- **Community Code**:
  Both packages maintain a tight connection to legacy software and community codes. You can find
  them in 
  **$AMUSE/src/amuse/community** and
  **$NEMO/usr**
  resp.

Examples
--------

Creating a Plummer sphere
~~~~~~~~~~~~~~~~~~~~~~~~~

Here we create a Plummer sphere, in virial units, in NEMO, and display an X-VX projection on the sky
in a shell session:

.. code-block::

   source /opt/nemo/nemo_starts.sh

   mkplummer p100 100
   snapplot p100 xvar=x yvar=vx

or in the style of using pipes this can be a one liner

.. code-block::

   mkplummer - 100 | snapplot - xvar=x yvar=vx

And in AMUSE the following python session can do something similar:

.. todo:: figure out the right py-gnuplot

.. code-block::

   from amuse.units import units
   from amuse.units import nbody_system
   from amuse.ic.plummer import new_plummer_sphere

   convert_nbody = nbody_system.nbody_to_si(100.0 | units.MSun, 1 | units.parsec)
   plummer = new_plummer_sphere(1000, convert_nbody)


   plotter = Gnuplot.Gnuplot()
   plotter.splot(plummer.position.value_in(units.parsec))

Installation
~~~~~~~~~~~~

For the benefit of NEMO users, AMUSE can usually be installed *easily* as follows:

.. code-block::

   pip install amuse


but this can take a while as it finds the right dependencies and needs to compile
some. Or will also easily fail if you have not set up MPI or other tools needed.

A potentially faster way is to first install
the AMUSE frame work and then the selected module(s):

.. code-block::

   pip install amuse-framework
   pip install amuse-seba amuse-brutus


There are many more details in the
`AMUSE installation manual <https://amuse.readthedocs.io/en/latest/install/index.html>`_.


