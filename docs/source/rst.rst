RST reminders
=============
	

-----
Lists
-----
- Bullet are made like this
- A link to the  `NEMO <https://github.com/teuben/nemo>`_ repo, or
  `github pages <https://teuben.github.io/nemo/>`_.
- Links to the man pages:
  `github <https://teuben.github.io/nemo/man_html/index.html>`_ or
  `local  <../../../man_html/index.html>`_
- Links to some old examples:
  `github <https://teuben.github.io/nemo/examples/index.html>`_ or
  `local  <../../../examples/index.html>`_

- Links to some Restructuredtext RST guides:
    * https://sphinx-tutorial.readthedocs.io/step-1/
    * https://www.writethedocs.org/guide/writing/reStructuredText/
    * https://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html
    
  
- Point levels must be consistent
    * Sub-bullets
        + Sub-sub-bullets
- Lists

Term
    Definition for term
Term2
    Definition for term 2

:List of Things:
    item1 - these are 'field lists' not bulleted lists
    item2
    item 3

:Something: single item
:Someitem: single item

------------
Code blocks
------------

There are three equivalents: ``code``, ``sourcecode``, and ``code-block``.

.. code:: python

   # some python code
   import os
   print(help(os))
   if True:
      print("yes")
   else:
      print("no")
   
.. sourcecode::

  # Equivalent

.. code-block::

  # Equivalent

.. code-block:: bash


   # some bash code
   if [ -e /tmp ]; then
      echo "You have /tmp"
   fi

------
Images
------

Here is an image ...

..  .. image::  ../examples/eagle_1.png

from the eagle example.



------
Tables
------

+--------+--------+--------+
| Time   | Number | Value  |
+========+========+========+
| 12:00  | 42     | 2      |
+--------+--------+--------+
| 23:00  | 23     | 4      |
+--------+--------+--------+

Math
----

Some example math, taken from latex

.. math::

   \ddot{ {\bf r}}_i \, = \, -G \sum_{j=1;\, j \not = \,i}^N {m_j \,({\bf r}_i - {\bf r}_j)  \over {(r_{ij}^2 + \epsilon^2)^{3/2} 


or inline  :math: `c^2 = a^2 + b^2`
