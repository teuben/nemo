STARLAB
-------



Installation
~~~~~~~~~~~~

For the benefit of NEMO users, STARLAB can usually be installed as follows:

.. code-block::

   mknemo starlab

Examples
~~~~~~~~

These examples were taken from http://www.sns.ias.edu/~starlab/examples/,
but an up-to-date list may be found in the file EXAMPLES in the Starlab
distribution. If installed within NEMO, this should be in
$NEMO/local/starlab/EXAMPLES.


For details on a specific program, type

.. code-block::

        program --help


- Create a linked list of 100 equal-mass nodes of unit total mass

.. code-block::


        makenode -n 100 -m 1

- Create a system of 100 nodes with a Salpeter mass spectrum with masses
in the range 0.5 to 10

.. code-block::

        makenode -n 100 | makemass -f 1 -x -2.35 -l 0.5 -u 10

- Create a system of 100 nodes with a mass spectrum and evolve the stars
without dynamics

.. code-block::

        makenode -n 100 | makemass -f 1 -x -2.35 -l 0.5 -u 10 \
                        | ...(to come)...

- Create a 500-particle Plummer model, with numbered stars, scaled to
standard dynamical units

.. code-block::

        makeplummer -n 500 -i

- Create a 500-particle W0 = 5 King model, with numbered stars, unscaled

.. code-block::
   
        makeking -n 500 -w 5 -i -u

- Create a 500-particle W0 = 5 King model with a Miller-Scalo mass
spectrum between 0.1 and 20 solar masses, then rescale to unit total
mass, total energy -0.25, and virial ratio 0.5 and display the results
graphically

.. code-block::
   
        makeking -n 500 -w 5 -i -u \
    	| makemass -F Miller_Scalo -l 0.1 -u 20 \
    	| scale -m 1 -e -0.25 -q 0.5 \
    	| xstarplot -l 5 -P .5

- Create a 500-particle W0 = 5 King model with a Miller-Scalo mass
spectrum between 0.1 and 20 solar masses, add in a 10 percent 1-10 kT
binary population, then rescale to unit total mass, total energy
(top-level nodes) -0.25, and virial ratio (top-level nodes) 0.5, and
finally verify the results by analyzing the final snapshot

.. code-block::

        makeking -n 500 -w 5 -i -u \
    	| makemass -f 2 -l 0.1 -u 20 \
    	| makesecondary -f 0.1 -l 0.25 \
    	| scale -m 1 -e -0.25 -q 0.5 \
    	| makebinary -l 1 -u 10 \
    	| sys_stats -n

- Evolve this model without stellar evolution for 100 dynamical times,
with log output every dynamical time and snapshot output every 10
dynamical times, with a self-consistent tidal field, removing escapers
when they are more than two Jacobi radii from the cluster center

.. code-block::
   
        makeking -n 500 -w 5 -i -u \
    	| makemass -f 2 -l 0.1 -u 20 \
    	| makesecondary -f 0.1 -l 0.25 \
    	| makebinary -l 1 -u 10 \
    	| scale -m 1 -e -0.25 -q 0.5 \
    	| kira -t 100 -d 1 -D 10 -Q -G 2

- Create a King model with a power-law mass spectrum and a binary
population, then evolve it with stellar and binary evolution

.. code-block::
   
        makeking -n 500 -w 5 -i -u \
    	| makemass -f 1 -x -2.0 -l 0.1 -u 20 \
    	| makesecondary -f 0.1 -l 0.1 \
    	| add_star -Q 0.5 -R 5 \
    	| scale -M 1 -E -0.25 -Q 0.5 \
    	| makebinary -f 1 -l 1 -u 1000 -o 2 \
    	| kira -t 100 -d 1 -D 10 -f 0.3 \
                          -n 10 -q 0.5 -Q -G 2 -B

- Perform a series of 100 3-body scattering experiments involving an
equal-mass circular binary and a double-mass incomer, with impact
parameter equal to the binary semimajor axis, relative velocity at
infinity half that needed for zero total energy, and all other
parameters chosen randomly, and display the results as a movie

.. code-block::
   
        scatter3 -m 0.5 -e 0 -M 1 -r 1 -v 0.5 \
                 -n 100 -C 5 -D 0.1 \
    	| xstarplot -l 4

- Compute cross-sections for interactions between a circular binary with
component masses 0.75 and 0.25 and an incoming star of mass 1 and
velocity at infinity 0.1, all stars having radius 0.05 binary
semimajor axes

.. code-block::
   
        sigma3 -d 100 -m 0.25 -e 0 -M 1 -v 0.1 \
               -x 0.05 -y 0.05 -z 0.05

- Create a scattering configuration involving a head-on collision
between a circular binary and a stable hierarchical triple, and verify
the result

.. code-block::

        makescat -M 1.5 -r 0 -v 1 -t -a 1 -e 0 \
                 -p -a 1 -e 0 -p1 -a 0.1 -e 0 \
    	| flatten | make_tree -D 1 | pretty_print_tree

- Create a scattering configuration involving a head-on collision
between a circular binary and a stable hierarchical triple, and
integrate it forward in time

.. code-block::
   
        scatter -i "-M 1.5 -r 0 -v 1 -t -a 1 -e 0 \
                -p -a 1 -e 0 -p1 -a 0.1 -e 0" \
    	        -t 100 -d 1 -v
