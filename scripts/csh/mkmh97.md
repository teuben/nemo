#    Interactions between Equal Mass Plummer Spheres

Based on the paper by Makino & Hut (1997) - https://ui.adsabs.harvard.edu/abs/1997ApJ...481...83M/abstract - this
example lets two equal mass Plummer (1911) interact with either other. Two options
are given: a (near) collision or a (near) circular orbit.

## Input parameters

All input parameters are of the form *keyword=value*, like most NEMO programs, with
the exception there is no syntax or validity of the keyword checking!

As in NEMO, units are virial units.

1. **run=**:  an identifying name used in all the derived filenames, e.g. run0, which then
creates files like run0.1 and run0.4.log etc.  Warning: if you specify an already
used simulation, that run is overwritten!

2. **nbody=**: number of bodies per plummer sphere. Notice that nbody=1 is also allowed,
allowing you to check on Kepler orbits. [1000]

3. **tstop=**: stopping time. Should be several times r0/v0   [50]

4. **step=**: step time when snapshots are stored. 1 is probably ok, for movies you probably need 0.1 [1.0]

5. **v0=**: initial impact speed between the two spheres.
Use positive values for a (near) collision. Negative values are reserved for a (near) circular orbit, although
with the right choices of **v0** and **r0** a circular orbit can also be achieved.  [1.0]

6. **rp=**: initial impact parameter. 0 means a head-on collision [0]

7. **r0=**: initial impact distance between the two spheres [10]

8. **eps=**: gravitational softening. There are some codes in NEMO that allow negative softening,
   in which case a Post-Newtonian (PN) approximation is used. This is outside the realm of this study. [0.05]
   
9. **kmax=**: parameter to control the timestep = 1/(2^kmax) [8]

10. **hack=**: Set to 1 if hackcode1 (the original 1986 Barnes & Hut treecode) has to be used. It's a bit slower than
   the default gyrfalcON, but on some machines the latter may not compile. [1]

## Plotting Examples

In the examples below the name *run0" is used, since that's the default for **run=**, however, if you want to 
preserve the data, be sure to use another run name!

1. Plot some initial conditions. Here is an X-Y plot, and an X-VX plot:

      snapplot run0.3 xrange=-16:16 yrange=-16:16
      snapplot run0.3 xrange=-16:16 yrange=-2:2 yvar=vx

2. Plot some of the evolution:

      # plot energy vs. time - how well is energy conserved?
      tabplot run0.4.etot

      # histogram of energy, easier to see fractional conservation
      tabhist run0.4.etot 2

      # evolution in X-Y projection
      snapplot run0.4 xrange=-16:16 yrange=-16:16
	  
	  # evolution but coloring each galaxy different
      snapplot run0.4 xrange=-16:16 yrange=-16:16 color='i<1000?0.1:0.2'
	  
      # the orbit of star #100 (try a few stars, some are bound, some escape)
      snapplot run0.4 xrange=-4:4 yrange=-4:4 trak=t visib=i==100
	  
	  # the evolution in X-VX space (rather interesting)
      snapplot run0.4 xrange=-16:16 yrange=-2:2 xvar=x yvar=vx

      # the evolution in Radius-Angular momentum space
      snapplot run0.4 xrange=0:16 yrange=-8:8 xvar=r yvar=jz
	  
3. Make a CCD frame so we can compare them to telescope images

      # a CCD frame of the final snapshot, and convert to a FITS file, to 
	  # to view it in ds9
	  snapgrid run0.4 run0.4.ccd xrange=-16:16 yrange=-16:16 times=50
      ccdfits run0.4.ccd run0.4.fits ndim=2 radecvel=t
	  ds9 run0.4.fits
	  
The quality of the image in [ds9](https://sites.google.com/cfa.harvard.edu/saoimageds9/download)
will depend strongly on the number of particles
in the galaxy. The default image only has 64 x 64 pixels, and with only
2000 particles in this default simulation there will be lots of pixels
with 0 stars.

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

9. Estimate how many particles we need to see shells, and how many to trace the full
   dynamic range of the galaxy as for example can be seen in the CenA image.
