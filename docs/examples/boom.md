# Orbits in N-body potentials ("BOOM")

BOOM = Bunch of Orbit Maps

The goal is to analyze a set of orbits. There are several ways how they can be made in NEMO:

1. In a (static) analytical potential (2D or 3D)
   1. Random initial conditions (mkorbit, orbint)
   2. Initial Conditions from a distribution function that represents the potential
   
2. In a static N-body potential. This would be useful to check the accuracy of an orbit
   compared to the analytical one
   
3. In a dynamical N-body simulation
   1. stable model (e.g. plummer)
   2. adiabatic changing model (e.g. larger softening)
   3. violent model (e.g. collapse simulation)

## Notes 

Using the Particle-Attribute-Time (P-A-T) notation

1. nemo snapshot:    [T][A][P]
2. tipsy, Body       [T][P][A]
3. orbits:           [P][A][T]

and an example:

1. gyrfalcon produces snap files, one huge file (the TAP, or should we call it PAT)

        stoo  p100.snap  p100.orb
        orboom p100.snap p100.boom
	  
        snappat ??

2. bonsai produces tipsy files, one file per snapshot


        ./tipsy_to_snap junk1.snap junk1_000*
         stoo junk1.snap junk1.orb
	   
or in one step:

        ./tipsy_to_snap - junk1_000* | stoo - junk1.orb nsteps=99


Note that the body stores a snapshot like tipsy, a "P"-"A" style,
but the I/O routines store the snapshot like a "A"-"P"

### nbody

Make a 1000 body Plummer sphere, and integrate for 2000 steps, saving each step

     mkplummer p3 1000 seed=123
     gyrfalcON p3 p3a.out eps=0.1 kmax=4 step=0.0625 tstop=125 give=mxvap > p3a.log
     

### Comparing simulations

One can also run a second simulation, e.g. with finer steps, or better force accuracy

     gyrfalcON p3 p3b.out eps=0.1 kmax=5 step=0.0625 tstop=125 give=mxvap > p3b.log

and compare the orbits visually, side by side in two panels:

     snapcopy p3a.out - i=10 | snapmerge - - | snapplot3 - yapp=1/xs
     snapcopy p3b.out - i=10 | snapmerge - - | snapplot3 - yapp=2/xs

![alt text](boom1.png "Comparing two orbits")

or another way to compare a selected observable and look at the whole ensemble of particles. Here the
three quarticles in "x" are plotted as function of time:

     snapcmp p3a.out p3b.out | tabplot - 1 3,4,5

Energy conservation can be viewed as follows:

     tabplot p3a.log 1 2
     tabstat p3a.log  2 qac=t

and energy conservation for a specific particle:

     snapplot p3a.out xvar=t yvar=etot visib=i==10

### convert to orbits

A snapshot can be "transposed" into a series of orbits, one for each particle. We have two formats for
this:

     stoo p3.out p3.orb
     orboom p3.out p3.boom odm=f

Plotting an orbit

     csf p3.boom - Orbit 11 | orbplot -
     csf p3.orb  - Orbit 11 | orbplot -
     snapcopy p3.out - i=10 | snapmerge - - | snapplot -

Plotting an orbit density map (ODM) :

     snapcopy p3.out - i=10 | snapmerge - - | snapgrid - - | ccdplot -

or if you just prefer a binary map with 0 and non-zero:

     snapcopy p3.out - i=10 | snapmerge - - | snapgrid - - mean=t | ccdplot -

