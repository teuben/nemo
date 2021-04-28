.. _using:

Using NEMO
==========

In order to use NEMO programs, NEMO will need to modify your
shell environment


.. code-block:: bash

	source /opt/nemo/nemo_start.sh

for users of the ``bash`` or ``zsh`` shell, and assuming NEMO
was installed in ``/opt/nemo/``

.. code-block:: bash

	source /opt/nemo/nemo_start.csh

for users of the ``csh`` or ``tcsh`` shell.

After this the following commands should work:

.. code-block:: bash

	tsf --help
	man tsf

Updating NEMO
-------------

Although a full update is out of scope for the discussion here, a common case is
when a program is not available, or a collaborator had made a new or updated program
available via github.  The following procedure generally works, assuming you have
write permission in $NEMO:

.. code-block:: bash

   mknemo -u tsf
   tsf --help

and an updated version should now be available (check the ``VERSION=`` default)

	
