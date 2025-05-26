#! /usr/bin/env python
#
#   Example of integrating a Plummer sphere with AMUSE
#
# See also $AMUSE_DIR/doc/interactive_tutorial/06-Using_a_community_code.ipynb
#

NBODY = 128
SEED  = 123
TIME  = 10
BASE  = 'plummer'
FMT   = 'txt'

filename = f"{BASE}.{FMT}"

import numpy as np
from amuse.units import nbody_system
from amuse.ic.plummer import new_plummer_model
from amuse.community.bhtree import Bhtree
# from amuse.community.bhtree.interface import BHTree   # why this also possible?
from amuse.io import write_set_to_file

np.random.seed(SEED)
stars = new_plummer_model(NBODY)
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
write_set_to_file(stars_in_gravity, filename, format=FMT)
print(f"Wrote {filename} at {gravity.model_time} using {gravity.model_name} in {FMT} format, nbody={NBODY}, seed={SEED}")
gravity.stop()
