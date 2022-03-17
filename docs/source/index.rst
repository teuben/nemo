Welcome to NEMO's documentation!
================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. math::

   \ddot{ {\bf r}}_i \, = \, -G \sum_{j=1;\, j \not = \,i}^N {m_j \,({\bf r}_i - {\bf r}_j)  \over {(r_{ij}^2 + \epsilon^2)^{3/2} } }


.. todo::  This draft is a work in progress, Most of the old
   *NEMO Users and Programmers Guide* has been converted here, although not all sections
   have been included yet.
	   
   Other online entry points for NEMO are:
   `github pages <https://teuben.github.io/nemo/>`_ and
   `github code <https://github.com/teuben/nemo>`_ .  We also keep an index to
   `(unix) man pages for all programs <https://teuben.github.io/nemo/man_html/index1.html>`_.

   Although this manual should be on
   https://astronemo.readthedocs.io
   we also keep a local copy 
   `readthedocs on astroumd <https://www.astro.umd.edu/~teuben/nemo/readthedocs/html>`_.

.. note::  **(*)** are sections that have not been fully cleaned up and can contain old latex markup

.. toctree::
   :numbered:
   :maxdepth: 2

   intro
   iface
   filestr
   graphics	      
   examples
   using
   .. setup
   install
   progr
   potname
   coordsys
   .. dirs
   .. gls
   trouble
   .. history
   .. bench
   refs
   codes
   glossary
   rse
   todo	      
   rst

Indices and tables
==================

* :ref:`modindex`
* :ref:`search`


.. include:: lastmod.rst
