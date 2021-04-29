Potentials
==========

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
