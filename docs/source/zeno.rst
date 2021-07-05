ZENO and NEMO
=============


`ZENO <http://www.ifa.hawaii.edu/faculty/barnes/zeno/index.html>`_
is Josh Barnes' version of an earlier version of NEMO that he
continues to develop.
We also keep a
`zeno manual page <https://teuben.github.io/nemo/man_html/zeno.1.html>`_
highlighting some differences.

Installation
~~~~~~~~~~~~

For the benefit of NEMO users, ZENO can usually be installed as follows:

.. code-block::

   cd $NEMO/usr/zeno
   make zeno
   ln -s zeno_pjt zeno
   source zeno_start.sh
   cd zeno
   make -f Zeno

Notice that this will install two versions of the ZENO source tree, and with
a symbolic link one of them needs to be picked to be there official one.

