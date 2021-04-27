.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Mar 31 15:30:17 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to NEMO's documentation!
================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

----
Lists
----
- Bullet are made like this
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

   import os
   print(help(os))

.. sourcecode::

  # Equivalent

.. code-block::

  # Equivalent

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
----------------
.. math::

   \ddot{ {\bf r}}_i \, = \, -G \sum_{j=1;\, j \not = \,i}^N {m_j \,({\bf r}_i - {\bf r}_j)  \over {(r_{ij}^2 + \epsilon^2)^{3/2} } }



Indices and tables
==================

* :ref:`modindex`
* :ref:`search`