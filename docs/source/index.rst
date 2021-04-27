.. NEMO documentation master file, created by
   sphinx-quickstart on Wed Mar 31 15:30:17 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to NEMO's documentation!
================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. math::

   \ddot{ {\bf r}}_i \, = \, -G \sum_{j=1;\, j \not = \,i}^N {m_j \,({\bf r}_i - {\bf r}_j)  \over {(r_{ij}^2 + \epsilon^2)^{3/2} } }


User: Developing
----------------
.. toctree::
   :maxdepth: 2

   intro
   filestr
   examples
   glossary

Programmer: developing
----------------------
.. toctree::
   :maxdepth: 2

   glossary

	      
Background on RST
=================
	      

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

   \ddot{ {\bf r}}_i \, = \, -G \sum_{j=1;\, j \not = \,i}^N {m_j \,({\bf r}_i - {\bf r}_j)  \over {(r_{ij}^2 + \epsilon^2)^{3/2} } }


Organization
============

This is the start of a chapter, which will contain sections, below which subsections. maybe even deeper.


MySection
---------

And inside a section we will allow subsections.


MySubSection
~~~~~~~~~~~~

and do we need more in NEMO?



Indices and tables
==================

* :ref:`modindex`
* :ref:`search`
