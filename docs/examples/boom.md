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

     mkplummer p3 1000
     gyrfalcON p3 p3.out eps=0.1 kmax=4 step=0.0625 tstop=125 give=mxvap > p3.log
     4.30user 0.09system 0:04.41elapsed 99%CPU

### convert to orbits

     stoo p3.out p3.orb
     orboom p3.out p3.boom odm=f

### plot-1

     csf p3.boom - Orbit 11 | orbplot -

### plot-2

     snapcopy p3.out - select=i==10 | snapmerge - - | snapplot -
     snapcopy p3.out - select=i==10 | snapmerge - - | snapplot3 -
     snapcopy p3.out - select=i==10 | snapmerge - - | snapgrid - - | ccdplot -

x
