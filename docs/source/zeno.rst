ZENO and NEMO
=============


`ZENO <http://www.ifa.hawaii.edu/faculty/barnes/zeno/index.html>`_
is Josh Barnes' version of an earlier version of NEMO that he
continues to develop.
We also keep a
`zeno manual page <https://teuben.github.io/nemo/man_html/zeno.1.html>`_
highlighting some differences.

.. warning::
   Adding ZENO to NEMO will result in some programs that have duplicated names.

Installation
~~~~~~~~~~~~

For the benefit of NEMO users, ZENO can usually be installed as follows:

.. code-block::

   cd $NEMO/usr/zeno
   make zeno

This will currently download two repos:   zeno_jeb and zeno_pjt. Pick one by
using a symlink to become the official one for the install:

.. code-block::

   ln -s zeno_pjt zeno
   source zeno_start.sh
   cd zeno
   make -f Zeno

Now various ZENO commands are available:

.. code-block::

   ls $ZENOPATH/bin

