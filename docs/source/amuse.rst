AMUSE
-----

AMUSE (Astrophysical MUltipurpose Software Environment) originates some ideas
from its predecessors: ACS, StarLab and NEMO, but uses the python language.
Another feature of AMUSE is that
python is also the *glue* shell between legacy codes that can orchestrate
simulations taking components from different codes, whereas in NEMO legacy codes
have a NEMO CLI interface, at best.

For seasoned
`AMUSE <https://amusecode.org>`_
users, here we highlight some differences between the two, and give some examples
how to achieve the same task in NEMO and AMUSE.


Differences
~~~~~~~~~~~

- **Shell**:
  NEMO uses a Unix shell, AMUSE uses python (ipython, jupyter, ...).  A neat way to start
  an interative amuse is via:  `ipython -profile amuse`

- **Community Code**:
  Both packages maintain a tight connection to legacy software and community codes. You can find
  them in 
  **$AMUSE/src/amuse/community** and
  **$NEMO/usr**
  resp., though the latter has some supporting script in **$NEMO/src/scripts/mknemo.d**

- **Units**:
  NEMO uses dimensionless values, and units are implied.
  Most programs actually use virial units (a.k.a N-body units, or Henon units) where G=1,
  but there are a few programs
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


Examples: Creating a Plummer sphere
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we create a Plummer sphere, in virial units, in NEMO, and display an X-VX projection on the sky
in a shell session:

.. code-block::

   source /opt/nemo/nemo_start.sh

   mkplummer p100 100
   snapplot p100 xvar=x yvar=vx

or in the style of using pipes this can be a one liner. Here is that example, and a few
followups with grey scale and contour plots:

.. code-block::

   mkplummer - 100 | snapplot - xvar=x yvar=vx
   mkplummer - 10000 | snapgrid - - xvar=x yvar=vx   | ccdplot  -
   mkplummer - 10000 | snapgrid - - xvar=x yvar=vx   | ccdplot  - 0.01,0.1,0.3,0.6,0.9
   mkplummer - 10000 | snapgrid - - xvar=x yvar=vx   | ccdsmooth -  - | ccdplot  - 0.01,0.1,0.3,0.6,0.9

And in AMUSE the following python session can do something similar:

.. todo:: figure out the right py-gnuplot

.. code-block::

   from amuse.units import units
   from amuse.units import nbody_system
   from amuse.ic.plummer import new_plummer_sphere

   convert_nbody = nbody_system.nbody_to_si(100.0 | units.MSun, 1 | units.parsec)
   plummer = new_plummer_sphere(1000, convert_nbody)

   # gnuplot
   plotter = Gnuplot.Gnuplot()
   plotter.splot(plummer.position.value_in(units.parsec))

   # matplotlib
   x=p[:,0].value_in(units.pc)
   y=p[:,1].value_in(units.pc)
   plt.plot(x,y,'ok')


The AMUSE manual has some
`NEMO I/O examples <https://amuse.readthedocs.io/en/latest/reference/fileformat.html#nemo>`_.

Installation
~~~~~~~~~~~~

AMUSE can usually be installed *easily* as follows (as a user):

.. code-block::

   pip install amuse

but this can take a while as it finds the right dependencies and needs to compile
massive amounts of code. Some of these can easily fail if you don't have the correct
`prerequisites <https://amuse.readthedocs.io/en/latest/install/howto-install-AMUSE.html>`_
(e.g. MPI).

A potentially faster way is to first install
the AMUSE frame work and then the selected module(s):

.. code-block::

   pip install amuse-framework
   pip install amuse-bhtree amuse-seba amuse-brutus

There are many more details in the
`AMUSE installation manual <https://amuse.readthedocs.io/en/latest/install/index.html>`_.


As a developer, it is easier to install the source code via github, and follow the
following recipe somewhere in your workflow (skipping details where your python/virtual environment
is)

.. code-block::

   git clone https://github.com/amusecode/amuse
   cd amuse
   pip install -e .
   make bhtree.code


this is the recommended way for NEMO users.
