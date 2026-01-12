The command "mknemo amuse" is now aimed to steer the install of AMUSE within NEMO

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

These benchmarks (see Makefile) don't use the same input data yet. TBD.

Timings are from `/usr/bin/time` on an Ultra 7 155H.


* bench1:   creating a large 1M Plummer sphere in memory (0), and also writing it (1)

      a0:  6.37user 7.06system 0:11.59elapsed 115%CPU
      a1: 43.93user 7.55system 0:49.78elapsed 103%CPU  (writing seems expensive in AMUSE)
      n0:  2.42user 0.47system 0:02.89elapsed 99%CPU
      n1:  2.38user 0.74system 0:03.13elapsed 99%CPU 


* bench2:   integrating an N=2048 Plummer sphere for 640 timesteps using **amuse-bhtree** (a), comparing to **hackcode1** (n),
            but also compare to **nemo-bhtree** itself (b)


      a:  6.26user 8.58system 0:13.58elapsed 109%CPU 
      n:  3.61user 0.00system 0:03.61elapsed 99%CPU 
      b:  1.45user 0.02system 0:01.47elapsed 100%CPU 

Note that bench2b can't be compared yet, "as is" snapshots are saved at near-integer times
1.066532, .... 9.076927, and of course nothing over 10.


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


