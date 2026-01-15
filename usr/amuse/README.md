# AMUSE

The command "mknemo amuse" is now aimed to steer the install of AMUSE within NEMO

## AMUSE in NEMO

Here are the components in NEMO that were added/modified to accomodate this
AMUSE "integration".

- scripts here:   amuse_convert.py and amuse_integrate.py
- bhtree (mknemo -a bhtree) to compare with the version(s) in AMUSE
- rebound via reb_integrate - can be compared to the version in AMUSE
- docs:  amuse.1 man page  and amuse.rst for the RtD manual
- amuse (helper) install via mknemo

## 2025 notes:

Early 2025, although the new developer install (which is laborious) works, there is the old
less intrusive user based *pip install* .  Here are the steps in a fresh NEMO install:

```
cd $NEMO
make python
source anaconda3/python_start.sh
mknemo amuse
pip install amuse-framework
pip install amuse-bhtree
cd $NEMO/usr/amuse
make test1
```

### benchmarking

These benchmarks (see Makefile) don't use the same input data yet. TBD. Neither confirmed
if they statistically are the same. But otherwise show some interesting numbers.

Timings are from `/usr/bin/time` on an Ultra 7 155H.


* bench1:   creating a large 1M Plummer sphere in memory (0), and also writing it (1)

                           d76                                      jansky
      a0:  3.95user 4.35system 0:08.31elapsed 99%CPU    3.08user 0.96system 0:04.05elapsed 99%CPU
      a1: 33.75user 2.82system 0:36.61elapsed 99%CPUU  28.94user 1.14system 0:30.12elapsed 99%CPU
      n0:  2.42user 0.47system 0:02.89elapsed 99%CPU    2.01user 0.11system 0:02.14elapsed 99%CPU
      n1:  2.38user 0.74system 0:03.13elapsed 99%CPU    2.02user 0.20system 0:02.23elapsed 99%CPU
      g0:  4.95user 0.10system 0:05.06elapsed 99%CPU
      g1:  4.81user 0.19system 0:05.01elapsed 99%CPU

This shows I/O in AMUSE is overly expensive, AMUSE also has a huge overall system overhead,
but it does not go away when sitting in an ipython session.
CPU wise AMUSE is almost 2x slower in creating a Plummer sphere,
but oddly enough the equivalent \fIfalcon\fP tool is more than 2x slower than NEMO.


* bench2:   integrating an N=2048 Plummer sphere for 640 timesteps,
  	    using **amuse-bhtree** (a), comparing to **hackcode1** (n),
            but also compare to **nemo-bhtree** itself (b) and the zeno version (z)


      a:  6.26user 8.58system 0:13.58elapsed 109%CPU   3.43user 0.04system 0:03.69elapsed 94%CPU  
      n:  3.61user 0.00system 0:03.61elapsed 99%CPU    2.98user 0.00system 0:02.99elapsed 99%CPU
      b:  1.45user 0.02system 0:01.47elapsed 100%CPU
      z:                                               2.27user 0.00system 0:02.27elapsed 99%CPU
      g: (0.91user 0.00system 0:00.91elapsed 100%CPU)

A more careful comparison is needed if the parameters all agree, notably the opening angle
and the use of quadrupole corrections. How do we compare. Do we
calibrate on the force errors, which we know as function of opening angle for a tree code.

The gyrfalcon(g) is added for show, since it's an O(N) code, not O(NlogN),
and should not be compared to the "true" treecodes.

The amuse compiled BHTC is slow (50sec) and looks like it's totally wrong (E not conserved).
The nemo compiled bhtree might seem better, but doesn't save snapshots at the correct time. Energy
is not conserved well either, and a Plummer seems not maintained properly.

### alternative python workflows

### 1. native python with venv

On ubuntu I needed to add packages:

```
sudo apt install python3-venv ipython3 python3-pip

cd $NEMO
rm -rf amuse_venv
python3 -m venv amuse_venv
source amuse_venv/bin/activate
pip install amuse-framework
pip install amuse-bhtree
cd $NEMO/usr/amuse
make test1

```

### 2. anaconda3 python with venv

```
cd $NEMO
src/scripts/install_anaconda3 
source anaconda3/python_start.sh

python3 -m venv amuse_venv2
source amuse_venv2/bin/activate

pip install amuse-framework
pip install amuse-bhtree
cd $NEMO/usr/amuse
make test1

```

### 3. anaconda3 python with new style setup

```

cd $NEMO
src/scripts/install_anaconda3 
source anaconda3/python_start.sh

python3 -m venv amuse_venv3
source amuse_venv3/bin/activate

mknemo amuse
cd $NEMO/local/amuse

python3 -m pip install pip wheel

./setup install amuse-framework

# fails building mpi4py


### 4. anaconda3 python with new style conda setup

```

cd $NEMO
src/scripts/install_anaconda3 
source anaconda3/python_start.sh

mk_conda `pwd`/anaconda3/bin/conda > conda.rc
source conda.rc

mknemo amuse
cd $NEMO/local/amuse


conda create --name Amuse-env -y
conda activate Amuse-env
conda install -c conda-forge pip wheel

<lots stuff>


./setup install amuse-framework

# fails building mpi4py

## gists

* conda : https://gist.github.com/teuben/1e42049f020c11149b6c68b4dc1a14a0
* pip:    https://gist.github.com/teuben/e60844254172bbe5d5e9e10417ad5277 

## 2019 notes:

See https://amusecode.github.io/

On ubuntu some preconditions are needed (and a choice of MPI:  mpich vs.openmpi)


```
sudo apt-get install build-essential gfortran python-dev \
  libopenmpi-dev openmpi-bin \
  libgsl-dev cmake libfftw3-3 libfftw3-dev \
  libgmp3-dev libmpfr6 libmpfr-dev \
  libhdf5-serial-dev hdf5-tools \
  git
```


