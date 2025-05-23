#! /usr/bin/env python
#

N=128
TIME=10
FILENAME='p.txt'

from amuse.units import nbody_system
from amuse.ic.plummer import new_plummer_model
from amuse.community.bhtree import Bhtree
from amuse.io import write_set_to_file,read_set_from_file 
stars = new_plummer_model(N)
gravity = Bhtree()
stars_in_gravity = gravity.particles.add_particles(stars)
gravity.evolve_model(TIME | nbody_system.time)
write_set_to_file(stars_in_gravity, FILENAME, format='txt')


# epsilon_squared
# theta_for_tree
# dt_dia
# timestep
