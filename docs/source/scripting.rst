Scripting
---------

We will focus on **bash** scripting here, for an example see the
**mkmh97.sh** script.
A typical script will contain the following components

1. *hash-bang* the first line of the script, and optionally include ``nemo_functions.sh``
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
      eps=0.05       # softening                  #> SCALE 0:1:0.01

and then parsing could be achieved as follows

.. code-block:: bash

      for arg in "$@"; do
         export "$arg"
      done

3. optionally, but highly recommended, some matching **--HELP** tags that can be retrieved using
the --help command line option


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



Example:  Extracting results from run directories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A common workflow is to run a series of simulations and  walking over a multi-dimensional parameter space.
This usually results into running each simulations in its own run-directory, plus storing
parameters and possibly resulting values in a parameter file for later extraction. The directory name
could even contain the value of these parameters as well. Several common strategies have been seen in
the wild, we list three:

.. code-block::

     # 1. simply enumerated (e.g id=001,002,003,....) parameters stored inside
     run_${id}

     # 2. linear list of directories, parameters stored inside
     run_${a}_${b}_${c}

     # 3. hierarchical directories, parameters stored inside
     run/$a/$b/$c

Although easier to visually identify the values of the parameters in 2. and 3., they don't scale very
well if a new parameter is introduced.  In the first case a simple lookup table can be created using
``nemopars``,  thus making it easier to find which parameters are used in while run directory. Here's
an example:

.. code-block::

     nemopars id,a,b,c run_*/nemopars.rc  > run.pars



Summary
~~~~~~~

Summarizing, here are the recommended methods to maintain and extract NEMO variables.

1. **nemopars**: extract parameters from a bash-style rc file (python should also be able to use it)

2. **nemovar**:  get and set NEMOVAR variables

3. **show_vars**:   alias via nemo_functions.sh

   
