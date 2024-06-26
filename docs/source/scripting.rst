Scripting
---------


.. todo::  more scripting examples needed here

We will focus on **bash** scripting here, basing comments
on the **mkmh97.sh** script.
A typical script will contain the following components

1. *hash-bang* the first line of the script, and optionally include nemo_functions.sh
   where various convenience functions are maintained

   
.. code-block:: bash

     #! /usr/bin/env bash
     #
     source nemo_functions.sh

and make the script executable		

.. code-block:: bash

     chmod +x script_name

   
2. parse the command line, ideally with a NEMO like ``keyword=value`` -
   optionally add *qtrun* directives to also allow to run it with a GUI frontend

.. code-block:: bash

      run=run1       # directory and run ID       #> ENTRY 
      nbody=1000     # number of bodies           #> RADIO 1000,10000,100000,1000000

and the parsing could be achieved as follows

.. code-block:: bash

      for arg in "$@"; do
         export "$arg"
      done

3. optionally, but highly recommended, some matching **--HELP** tags that can be retrieved using
a --help command line option


.. code-block:: bash

      #--HELP
      a=1      # the a value      #> SCALE 0:10:1
      b=2      # the b value      #> SCALE 0:20:2
      #--HELP

      if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
         set +x
         awk 'BEGIN{s=0} {if ($1=="#--HELP") s=1-s;  else if(s) print $0; }' $0
         exit 0
      fi
		


4. optionally maintain a set of *nemopars* or *nemovar*

.. code-block::

     _pars=nemopars.rc

     save_vars="run nbody"

     show_vars $save_vars >> $_pars
     source  $_pars

     echo "a1=$a1" >> $_pars



5. the body of the script, shell variables, nemopars variables etc. should be used
where possible. 
