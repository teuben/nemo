The command "mknemo amuse" is now aimed to steer the install of AMUSE within NEMO

## spring 2025

Early 2025, although the new install (which is laborious) works, there is the old
less intrusive pip install.  Here are the steps in a fresh NEMO install:

```
cd $NEMO
make python
source anaconda3/python_start.sh
mknemo amuse
pip install amuse-framework
pip install amuse-bhtree
```

## older 2019 notes that need to be updated:

See https://amusecode.github.io/

On ubuntu some preconditions are needed (and a choice of MPI:  mpich vs.openmpi)

sudo apt-get install build-essential gfortran python-dev \
  libopenmpi-dev openmpi-bin \
  libgsl-dev cmake libfftw3-3 libfftw3-dev \
  libgmp3-dev libmpfr6 libmpfr-dev \
  libhdf5-serial-dev hdf5-tools \
  git

recommended are also

pip install numpy nose docutils mpi4py h5py

pip install scipy astropy jupyter pandas seaborn

# this can have issues
# pip install amuse

# then just the framework
pip install amuse-framework

# and whatever community codes you need

pip install amuse-brutus
pip install amuse-bhtree
pip install amuse-bse
pip install amuse-fractalcluster 
pip install amuse-gadget2 


Example comparing two codes:

from amuse.community.brutus.interface import Brutus
from amuse.community.bhtree.interface import BHTree
from amuse.datamodel import Particles
from amuse.units import nbody_system
from amuse.units import units

convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)


instance = BHTree(convert_nbody)
instance.parameters.epsilon_squared = 0.001 | units.AU**2

or:

instance = Brutus(convert_nbody)


stars = Particles(2)
sun = stars[0]
sun.mass = 1.0 | units.MSun
sun.position = [0.0,0.0,0.0] | units.m
sun.velocity = [0.0,0.0,0.0] | units.m / units.s
sun.radius = 1.0 | units.RSun
earth = stars[1]
earth.mass = 5.9736e24 | units.kg
earth.radius = 6371.0 | units.km 
earth.position = [1.0, 0.0, 0.0] | units.AU
earth.velocity = [0.0, 29783, 0.0] | units.m / units.s
instance.particles.add_particles(stars)

channel = instance.particles.new_channel_to(stars)

print(earth.position[0])
print(earth.position.as_quantity_in(units.AU)[0])
instance.evolve_model(1.0 | units.yr)
print(earth.position.as_quantity_in(units.AU)[0])
channel.copy()
print(earth.position.as_quantity_in(units.AU)[0])
instance.evolve_model(1.5 | units.yr)
channel.copy()
print(earth.position.as_quantity_in(units.AU)[0])

instance.stop()


