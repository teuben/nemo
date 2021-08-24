Potentials and Accellerations (*)
=================================

Here we lists a number of potentials, taken from
from CTEX comments in the
**$NEMO/src/orbit/potential/data** source code
directory.
Most NEMO programs that deal with potentials have three 
program keywords associated with 
potentials:

- **potname=**
  describes the name

- **potpars=**
  optional parameters

- **potfile=**
  optional associated filenames (or other textual information)

Each section below details a potential
and  explains the usage of the **potpars=** and **potfile=**
keywords. 
The section title is the actual **potname=** to be used for
this potential.
Mostly **G=1**, unless otherwise mentioned.

.. todo::  describe the newer **accelerations** from *falcON*

.. todo::  potentials auto-build from source code ???
   
.. include:: potctex.rst

Accellerations
--------------

This is a falcON addition. They are defined in ``$FALCON/src/public/acc``,
where their list (for installation) is defined in ``$FALCON/makepub``.


For a new potential, say ``GasPotential.cc``, add to ``$FALCON/makepub``:

.. code-block::

   acc_pub         :=
                         ...
                         $(ACC)GasPotential.so

and a proper dependency as well:

.. code-block::

   $(ACC)GasPotential.so: $(SACC)GasPotential.cc $(ACCT) $(defacc_h) $(makefiles)
                          $(MAKE_ACC)
