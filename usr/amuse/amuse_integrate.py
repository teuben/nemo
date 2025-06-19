#! /usr/bin/env python
#
#   Example of integrating a Plummer sphere with AMUSE
#   Also writes data with a few intermediate snapshots
#
# See also $AMUSE_DIR/doc/interactive_tutorial/06-Using_a_community_code.ipynb
#


import sys
import argparse
import numpy as np
from amuse.units import nbody_system
from amuse.io import write_set_to_file
from amuse.ic.plummer import new_plummer_model

my_version = "14-jun-2025"

my_help = [f"Version: {my_version}",
           "",
           "Integrate a snapshots using AMUSE, ",
           "The following codes (--code or -d) are allowed:",
           "   bhtree (default)",
           "   hermite",
           "   smalln",
           "   brutus",
           "Use ---debug (or -d) to see defaults",
           ]

def commandLine():
    parser = argparse.ArgumentParser(description="\n".join(my_help),
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--ifile',     '-i',  help="input snapshot name",    default=None)    # still ignored
    parser.add_argument('--ofile',     '-o',  help="output snapshot name",   default=None)
    parser.add_argument('--nbody',     '-n',  help="nbody if creating",      default=128)
    parser.add_argument('--seed',      '-s',  help='seed if creating',       default=123)
    parser.add_argument('--step',      '-t',  help='output dump time',       default=1)
    parser.add_argument('--stop',      '-T',  help='integration stop time',  default=10)
    parser.add_argument('--delta',     '-D',  help='integration step',       default=0.015625)
    parser.add_argument('--eps',       '-e',  help='softening',              default=0.025)
    parser.add_argument('--theta',     '-A',  help='treecode opening angle', default=0.75)
    parser.add_argument('--code',      '-c',  help='community code',         default='bhtree')
    parser.add_argument('--debug',     '-d',  help='add some debugging',     action="store_true")

    args = parser.parse_args()
    if args.debug:
        print(args)
    return args

#
args = commandLine()

# 
nbody = int(args.nbody)
seed  = int(args.seed)
code  = args.code
ofile = args.ofile
   
# code control
eps   = float(args.eps)   | nbody_system.length
tnow  = 0.0               | nbody_system.time
step  = float(args.step)  | nbody_system.time
tstop = float(args.stop)  | nbody_system.time
delta = float(args.delta) | nbody_system.time
theta = float(args.theta)

# create a Plummer sphere, and write the first snapshot
if seed != 0:
    np.random.seed(seed)
stars = new_plummer_model(nbody)
write_set_to_file(stars.savepoint(tnow), ofile, overwrite_file = True)
print("Writing snapshot",tnow)

# set up gravity solver

allowed_codes = ['bhtree', 'hermite', 'smalln', 'brutus']

if code == 'bhtree':
    from amuse.community.bhtree import Bhtree
    gravity = Bhtree()
    gravity.parameters.opening_angle   = theta
    gravity.parameters.timestep        = delta
    gravity.parameters.epsilon_squared = eps * eps
    gravity.parameters.dt_dia          = step
elif code == 'hermite':
    from amuse.community.hermite import Hermite
    gravity = Hermite()
    gravity.parameters.epsilon_squared = eps * eps
    gravity.parameters.dt_dia          = step
elif code == 'smalln':
    from amuse.community.smalln import Smalln
    gravity = Smalln()    
elif code == 'brutus':
    from amuse.community.brutus import Brutus
    gravity = Brutus()    
else:
    # @todo petar ph4
    print(f"cannot use {code}, allowed are {allowed_codes}")
    sys.exit(1)


print('GRAVITY PARAMETERS:',gravity.parameters)


stars_in_gravity = gravity.particles.add_particles(stars)

Qchan = True    # do we still need this ? -> yes, stars_in_gravity can't use savepoint()

# channel_from_gravity_to_framework = gravity.particles.new_channel_to(bodies)
# channel_from_gravity_to_framework.copy()

while gravity.model_time < tstop:
    gravity.evolve_model(gravity.model_time + (step))
    tnow = gravity.model_time
    if Qchan:
        stars_in_gravity.new_channel_to(stars).copy()
        stars.savepoint(timestamp=tnow)
        write_set_to_file(stars.savepoint(tnow), ofile, append_to_file = True)
    else:
        write_set_to_file(stars_in_gravity, ofile, append_to_file = True)        
    print("Writing snapshot",tnow)

print(f"Wrote {ofile} at {gravity.model_time} using {gravity.model_name} in amuse format, nbody={nbody}, seed={seed}")
gravity.stop()

if True:
    from amuse.support import literature
    literature.TrackLiteratureReferences.suppress_output()
