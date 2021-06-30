.. _using:

Using NEMO
==========

In order to use NEMO, you will need to modify your
shell environment, for example in the ``bash`` shell
this could be

.. code-block:: bash

	source /opt/nemo/nemo_start.sh

assuming NEMO was installed in ``/opt/nemo/``, and

.. code-block:: bash

	source /opt/nemo/nemo_start.csh

for users of a ``csh`` like shell. Normally this
line would be added to your ``~/.bashrc`` or ``/.cshrc`` file.

After this the following commands should work for you

.. code-block:: bash

	tsf --help
	man tsf

Updating NEMO
-------------

Although a full update is out of scope for the discussion here, a common case is
when a program is not available, or a collaborator had made a new or updated program
available via github.  The following procedure generally works, assuming you have
write permission in the ``$NEMO`` directory tree:

.. code-block:: bash

   mknemo -u tsf
   mknemo -h
   man mknemo		

and an updated version should now be available (check the value of the ``VERSION=``
value in the output of ``--help``).

Writing NEMO program programs is covered in :ref:`progr`, or see
also :ref:`install`.

	
