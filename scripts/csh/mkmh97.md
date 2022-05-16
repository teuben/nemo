#    Interactions between Equal Mass Plummer Spheres

Based on the paper by Makino & Hut (1997) - https://ui.adsabs.harvard.edu/abs/1997ApJ...481...83M/abstract - this
example lets two equal mass Plummer (1911) interact with either other. Two options
are given: a (near) collision or a (near) circular orbit.

## Input parameters

All input parameters are of the form *keyword=value*, like most NEMO programs, with
the exception there is no syntax or validity of the keyword checking!

As in NEMO, units are virial units.

1. **run**:  an identifying name used in all the derived filenames, e.g. run0, which then
creates files like run0.1 and run0.4.log etc.  Warning: if you specify an already
used simulation, that run is overwritten!

2. **nbody**: number of bodies per plummer sphere. Notice that nbody=1 is also allowed,
allowing you to check on Kepler orbits. [1000]

3. **tstop**: stopping time. Should be several times r0/v0   [50]

4. **step**: step time when snapshots are stored. 1 is probably ok, for movies you probably need 0.1 [1.0]

5. **v0**: initial impact speed between the two spheres.
Use positive values for a (near) collision. Negative values are reserved for a (near) circular orbit, although
with the right choices of **v0** and **r0** a circular orbit can also be achieved.  [1.0]

6. **rp**: initial impact parameter. 0 means a head-on collision [0]

7. **r0**: initial impact distance between the two spheres [10]

8. **eps**: gravitational softening. There are some codes in NEMO that allow negative softening,
   in which case a Post-Newtonian (PN) approximation is used. This is outside the realm of this study. [0.05]

## Plotting Examples

Plot some initial conditions:

      snapplot run0.3 xrange=-16:16 yrange=-16:16
      snapplot run0.3 xrange=-16:16 yrange=-2:2 yvar=vx

Plot some of the evolution:

      # plot energy vs. time - how well is energy conserved?
      tabplot run0.4.log

      # histogram of energy, easier to see fractional conservation
      tabhist run0.4.log 2

      # evolution in XY projection
      snapplot run0.4 xrange=-4:4 yrange=-4:4

      # the orbit of star #100
      snapplot run0.4 xrange=-4:4 yrange=-4:4 trak=t visib==100

      # the evolution in Radius-Angular momentum space
      snapplot run0.4.log xrange=0:4 yrange=-2:2 xvar=r xvar=jz
      
## Circular orbits?

For **v0<0** we can set up the two sytems in a circular orbit by launching from
(r0,0,0) and (0,v0,0), i.e. in clock wise motion. For an exact circular orbit,
the following example should work

      ./mkmh97.sh v0=-1 r0=2 nbody=1 eps=0
      snapplot run0.4 trak=t 
	  
and you will see the two particles chase each other on the same circular orbit. Pick a 
different v0 or r0 and this will not be true.

## TODO

1. define v0 as being from infinity? (cf. Makino & Hut paper)

2. allow circular orbit using negative **v0**.   The initial offset **r0** is then used
as the diameter of the circular orbit. This could be a fun way to study dynamical
friction (cf. Bontekoe & v Albada, White, ...)

3. Exact Newtonian solutions vs. Order Of Magnitude Estimates (OOME)

4. various sanity tests
  e.g. energy conservation as function of integration step

5. for v0 at r0, what is vi (v_infty)

6. write down the eq. for a circular orbit: what is the relation ship between r0 and v0

7. find v0 for the circular orbit at r0=4

8. integrate plunging orbit for different values of softening.
   make sure it's not an escaping orbit, i.e. E < 0
   Plot for example X vs. VX, or time (t)
