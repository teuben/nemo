#! /usr/bin/env python
#
# See also $AMUSE_DIR/doc/interactive_tutorial/06-Using_a_community_code.ipynb

N=128
TIME=10
FILENAME='p.txt'
FMT='txt'

from amuse.units import nbody_system
from amuse.ic.plummer import new_plummer_model
from amuse.community.bhtree import Bhtree
from amuse.io import write_set_to_file,read_set_from_file

stars = new_plummer_model(N)
gravity = Bhtree()

eps = 0.025 | nbody_system.length
tstop = 10  | nbody_system.time

gravity.parameters.opening_angle  = 0.75
#gravity.parameters.theta_for_tree  = 0.75     # this is confusing in interface.py
gravity.parameters.epsilon_squared = eps * eps
gravity.parameters.dt_dia          = 1.0      | nbody_system.time
gravity.parameters.timestep        = 0.015625 | nbody_system.time       # 1/64

stars_in_gravity = gravity.particles.add_particles(stars)
gravity.evolve_model(tstop)
write_set_to_file(stars_in_gravity, FILENAME, format=FMT)
print(f"Wrote {FILENAME} in {FMT} format")
