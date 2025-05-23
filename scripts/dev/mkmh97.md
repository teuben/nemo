#    Interactions between two Plummer Spheres

Loosely based on the paper by Makino & Hut (1997) -
https://ui.adsabs.harvard.edu/abs/1997ApJ...481...83M/abstract - the
**mkmh97.sh** script lets two Plummer (1911) interact with either
other. In the original Makino & Hut paper the masses were equal, but
the model, relative velocity and impact parameter were all parameters.

In this script we extend this with a few options:

1. a (near) collision or a (near) circular orbit. (controlled by the sign of **v0=**)

2. unequal masses, where the first galaxy is fixed at mass=1 (**m=**)

3. pick a different N-body integrator (**code=**)

4. different models, though for shell galaxies it will be useful to replace
   the Plummer with a colder stellar system [not yet implemented]

5. allow the first galaxy (with mass=1) to be a fixed analytical potential. This will
   however result in some other limitations to the analysis. (**fixed=**)

Some additional analysis is performed as well, see also **mh16.sh** for an example
for further time-dependant analysis.

## Input parameters

All input parameters are of the form *keyword=value*, like most NEMO programs, but with
the exception there is no check on syntax or validity of the keyword!  Use **--help**
or **-h**
as the first keyword to get up to date information on the keywords and their defaults.

As in many NEMO programs, units are virial (N-body) units.

1. **run=**:  run directory name, also the name used in all the derived filenames. For example with
   run=run0, you would see files like run0/run0.1 and run0/run0.4.log etc.  Warning: if you specify an already
   used simulation, the integration step is skipped, but the analysis is done with
   whatever analysis parameters are given again. [run0]

2. **nbody=**: number of bodies per plummer sphere. Notice that nbody=1 is also allowed,
   allowing you to check on Kepler orbits. [2048]

3. **m=**: mass of the 2nd galaxy where the first galaxy has fixed mass=1   [1]

4. **em=**: use equal mass particles? By default, each system has same number
   of particles, even if their massed are different. [0]

4. **zm=**: if set to 1, zero the masses of the companion and place a testparticle at it's origin [0]

4. **fixed=**: if set to 1, use a fixed potential for galaxy-1.  [0]

4. **potname=**: if **fixed=1**, this is the potname. [plummer]

4. **potpars=**: if **fixed=1**, this is the potpars. [0,1,3*pi/16]

4. **zerocm=**:  if set to true, apply center of mass to the input model *not implemented yet*

4. **step=**: step time when full snapshots are stored. 1 is probably ok,
   for movies you probably need 0.1.   For very large values of nbody a larger value for the step
   is probably adviced, unless you have a lot of disk space. Perhaps step=5.  [1.0]

5. **v0=**: initial impact speed between the two galaxies.  Use
   positive values for a (near) collision. Negative values are reserved
   for a (near) circular orbit, although with the right choices of **v0**
   and **r0** a circular orbit can also be achieved.
   [1.0]

6. **rp=**: initial impact parameter. 0 means a head-on collision. Note MH97 use
   two conventions. Pericenter distance and impact parameter. [0]

7. **r0=**: initial impact distance between the two spheres [16]

8. **eps=**: gravitational softening. There are some codes in NEMO
   that allow negative softening, in which case a Post-Newtonian (PN)
   approximation is used. This is outside the realm of this
   study. [1/32 = 0.03125]
   
9. **kmax=**: parameter to control the timestep = 1/(2^kmax) [7]

10. **code=**: Set to 0 if **hackcode1** (the original 1986 Barnes &
   Hut treecode) has to be used. It's a bit slower than the default
   **gyrfalcON** , but on some machines the latter may not
   compile. code=2 is reserved for GPU enabled machines where **bonsai2**
   has been compiled. This code is likely to be the fastests as long as it fits
   in the GPU memory.
   [1]

12. **seed=**: the usual NEMO seed value. 0=current time in seconds since 1970,
   -1=number of centiseconds since boot, -2=PID, -3=linux kernel entropy based.
   The command **date +%Y%m%d%H%M%S** gives a perhaps more memorable seed.
   [0]

13. **trim=**: set to 1 if you want to trim the simulation from all but the final
    (tstop) snapshot. [0]

13. **save=**: save the "final" (tstop=) plots in a subdirrectory called "movies"
    labeled with the tstop.   [1]

Parameters for (re)analysis (this happens if the **run=** already exists):

13. **tstop=**: stopping (or analysis on a re-run) time. Should be several times r0/v0
   [50]

13. **box=**: spatial plotting box size
   [32]

14. **r16=**: radius where m16 is measured.
   [16]

15. **vbox=**: velocity plotting size
   [2]

16. **npixel=**:  number of pixels in the CCD frames
   [128]

17. **power=**:  gamma factor for CCD plots
   [0.5]

18. **bsigma=**: asinh/log breakover point
   [0.0001]

19. **tplot=**: times to plot in the 3x3 evolution plot.
   [0,5,10,15,20,25,30,40,50]

19. **hackforce=**: which version of hackforce to use when computing force re-calculations

20. **yapp=**:   pick png, or vps (for yapp_pgplot) _ps for native ps

20. **progress=**:  spawn a GUI progress bar for the n-body evolution [0]

21. **debug=**:  not the usual NEMO debug=, but special for this script.
    0=nothing   1=set -x -e -u in bash
    [1]

## Running

The first time when you run the script and the **run** directory does not exist yet,
initial conditions are generated and the simulation is run to **tstop**.

Subsequent runs with an existing  **run** directory will re-analyse the simulation at
the value of **tstop**, which can be any of the dump-times that the simulation was
run with, so it does not need to be the initial **tstop** with which the simulation
was run.

The default run=run0 (with m=1 nbody=2048 tstop=50) will take around 30 seconds to run,
with the gyrfalcON code (code=1). With hackcode1 (code=0) it will be about 60 seconds
(on a 2020 style i7-1185G7, or where NEMOBENCH5 measured around 1000).

## Files

The following files should be present, the example is for run=run0:

     etot.hist.log           Logfile
     etot.hist.png           Histogram of total energy values
     etot.plot.png           Time evolution of total energy
     evolution-vr.plot.png   Time evolution of r-vr
     evolution-xy.plot.png   Time evolution of x-y
     final.3d.plot.png       UnFolded x-y-z view of the final ($tstop) 
     final.ccd               Final CCD of the whole system
     final.ccd.png
     final.fits
     final.plot.png
     final1.ccd              Final CCD of just Galaxy-1 (plus fits and point plot)
     final1.ccd.png
     final1.fits
     final1.plot.png
     final2.ccd              Final CCD of just Galaxy-2 (plus fits and point plot)
     final2.ccd.png
     final2.fits
     final2.plot.png
     final2.snap
     final2c.snap            Final snapshot of Galaxy-2 in its own center of mass frame
     final2c.tab             Radial profile (from radprof) ***
     final2cm.tab            Table of cumulative mass vs. radius
     final2u.ccd             Final bound particles of Galaxy-2
     final2u.snap
     fixed-path.tab          Integration of (single particle) radial orbit
     fixed1.tab
     fixed2.tab
     init.ccd.png            Initial conditions
     init.plot.png
     massg2g1.png            Cumulative mass of G2 around the center of G1
     nemopars.rc             Simulation parameters
     path-energy.png         Path of the binding energy
     path-g2.png             Path of mass of G2 in several ways
     path-pos.png            Path of the two galaxies in X;  1 in red, 2 in green
     path-vel.png            Path of the two galaxies in VX; 1 in red, 2 in green
     run0.1                  Galaxy-1, the one with mass=1 - snapshot
     run0.2                  Galaxy-2, the one with mass=$m - shapshot
     run0.3                  The initial conditions - snapshot 
     run0.4                  Simulation, as time intervals of $step - snapshot
     run0.4.etot             Table of time and total energy
     run0.4.g1.tab           Table of pos and vel of Galaxy-1           
     run0.4.g2.tab           Table of pos and vel of Galaxy-2
     run0.4.log              Logfile of the simulation
     run0.4.t.tab            Table of times
     run0.xv.tab             Table of time,x1,vx1,x2,vx2
     run0.xve.tab            Table of time,x1,vx1,x2,vx2,kin1,kin2,pot,etot
     run0.xvm.tab            Table of time,x1,vx1,x2,vx2,m16,m2    (only after m16.sh was run)

## Parameters

In addition to the input parameter, the following parameters are computed by the script,
and placed in the "rc" file:

1. **m1=**:  absolute mass of stars around G1 bound to G1. 

1. **m2=**:  absolute mass of stars around G2 bound to G1. 

2. **m16=**: relative mass of stars of G2 around the center of G1, within a radius of **r16**.

2. **etot=**:  orbital energy of bound particles of the last 10 snapshot

3. **detot=**:  variation of **etot**

2. **etot2=**:  orbital energy of bound particles corrected for mass-loss

2. **etot0=**:  orbital energy at time=0 (ideal)

4. **v=**:      correct input velocity (from v0), as computed from two systems at infinity

All parameters are stored in the **nemopars.rc** file, which can be queried with
the **nemopars** program.

## Pipeline Examples

Here is an example of running 16 instances of the same input parameters (m,v,n) but
with differently seeded models:

     m=0.05
     v=0.80
     n=8000
     r=$(nemoinp 1:16 format=%02d)
     for rr in $r; do
        mkmh97.sh run=run_${m}_${v}_${n}_${rr} m=$m v0=$v nbody=$n tstop=100 step=1 seed=-2
     done

after which querying all nemopars.rc files will result in a table, which could be plotted or further analyzed:

     nemopars m,v0,etot,detot,etot2 run_*_*_*_*/nemopars.rc
     # m v0 etot detot etot2
     0.05 0.80 0.00371626 9.70455e-05 0.00261469 
     0.05 0.80 0.00361105 0.000106836 0.0026038 
     0.05 0.80 0.00261421 0.000115724 0.00178161 
     0.05 0.80 0.00273166 7.53117e-05 0.00183988
     ...
  
## Plotting Examples

In the examples below the name **run0** is used, since that's the default for **run=**, however, if you want to 
preserve the data, be sure to use another run name! These are run inside of the run directory.

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
        snapplot run0.4 xrange=-16:16 yrange=-16:16 color='i<2048?0.1:0.2'
	  
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
in the galaxy. The default image only has 128 x 128 pixels, and with only
4096 particles in this default simulation there will be lots of pixels
with 0 stars. One can use **ccdfill** to patch this up a bit, perhaps in combination
with **ccdsmooth**.

## Circular orbits?

For **v0<0** we can set up the two sytems in a circular orbit by launching from
(r0,0,0) and (0,v0,0), i.e. in clock wise motion. For an exact circular orbit,
the following example should work

      rm -rf run0
      ./mkmh97.sh run=run0 v0=-1 r0=2 nbody=1 eps=0
      snapplot run0/run0.4 trak=t 
	  
and you will see the two particles chase each other on the same circular orbit. Pick a 
different v0 or r0 and this will not be true.

      ./mkmh97.sh run=run1 eps=0.05 nbody=1000
	  
Compare to Bontekoe & v Albada (1987MNRAS.224..349B) work.

## EXERCISES

1. define v0 as being from infinity? (cf. Makino & Hut paper)
   Noting that for position values of v0 the value of r0 will determine
   a minimum value.

2. allow circular orbit using negative **v0**.   The initial offset **r0** is then used
   as the diameter of the circular orbit. This could be a fun way to study dynamical
   friction (cf. Bontekoe & v Albada 1987MNRAS.224..349B   White, ...)

3. Exact Newtonian solutions vs. Order Of Magnitude Estimates.   Can you use this
   script to get an idea of the accuracy of the forces of tree vs. exact newton?

4. various sanity tests
   e.g. energy conservation as function of integration step

5. for v0 at r0, what is the asymptotic value 'v' at infinity?

6. write down the eq. for a circular orbit: what is the relation ship between r0 and v0

7. find v0 for the circular orbit at r0=4

8. integrate plunging orbit for different values of softening.
   make sure it's not an escaping orbit, i.e. E < 0
   Plot for example X vs. VX, or time (t)

9. Estimate how many particles we need to see shells, and how many to trace the full
   dynamic range of the galaxy as for example can be seen in the CenA image.

10. Run the script with two systems with initial conditions r0=0 and v0=0. What does this do? Hint: you
    may need a special seed= to avoid an interesting numerical problem (some codes will even crash)

## Example analysis

###  Equal masses?

Here are two commands comparing the analysis of two runs with different masses for the particles.

     ./mkmh97.sh run=run100 tstop=50 em=0 m=0.1 code=1 kmax=6 nbody=100000
     ./mkmh97.sh run=run101 tstop=50 em=1 m=0.1 code=1 kmax=6 nbody=100000

and plotting

     snapplot run100/run100.4 times=50 xrange=-32:32 yrange=-32:32 yapp=2/xs visib='i>100000'
     snapplot run101/run101.4 times=50 xrange=-32:32 yrange=-32:32 yapp=1/xs visib='i>100000'

a rough view shows that selecting x<0 might discimincate:

     snapcopy run100/run100.4 - times=50 select='i>100000 && x<0' | tsf -
     snapcopy run101/run101.4 - times=50 select='i>100000 && x<0' | tsf -

This gives  5578 vs.  587   or about 5.6 vs. 5.9% of particles of G1 are lost to G1.

and perhaps a better way

     snaptrim run100.4 snap1 times=50
     snapcenter snap1 snap2 'i<100000?-phi*phi*phi:0'
     snapcopy snap2 snap3 select='i>=100000'
     radprof snap3 mode=mass tab=t > tab1
     tabplot tab1 1 4 0 16 0 0.008

## CCD images

A few words on CCD images created with **snapgrid**:
For a given total mass of the input snapshot, the resulting
CCD image will conserve the total mass, or its "brightness".
The image stores a surface brightness, i.e. the mass per square pixel. This will
allow one to compare the intensities in images with different pixel sizes.
The program **ccdhead** will show the pixel size, and **ccdstat** will
report the total mass (Sum.Dx.Dy), which should be 1.0 for a default Plummer sphere:

     mkplummer - 1000 | snapgrid - - p1ka.ccd
     ccdhead p1ka.ccd
          Size:      64 64 1
          Cell:      0.0625 0.0625 2e+20  
          X-range:   -2 2
          Y-range:   -2 2
     ccdstat p1ka.ccd
          Sum and Sum*Dx*Dy*     : 238.592000 0.932000

the reason the total mass is not quite 1.0 is because a default Plummer sphere has particles
out to infinity, although by default **mkplummer** has rfrac=22.8. By setting this to the
default box size of 2 (or less) in snapgrid, you should now see all the mass in the CCD image,
although now one should call this a *Truncated Plummer Sphere*:

     mkplummer - 1000 rfrac=2 | snapgrid - - p1kb.ccd
     ccdstat p1kb.ccd
          Sum and Sum*Dx*Dy*     : 256.000000 1.000000

     ./mkmh97.sh run=run6 tstop=50 code=1 nbody=100000 m=0.01 kmax=6 npixel=128    0
     ./mkmh97.sh run=run7 tstop=50 code=1 nbody=100000 m=0.02 kmax=6 npixel=128    0.0002512	
     ./mkmh97.sh run=run5 tstop=50 code=1 nbody=100000 m=0.04 kmax=6 npixel=128    0.003552
     ./mkmh97.sh run=run9 tstop=50 code=1 nbody=100000 m=0.05 kmax=6 npixel=128    0.00585   
     ./mkmh97.sh run=run8 tstop=50 code=1 nbody=100000 m=0.08 kmax=6 npixel=128    0.01608

     ./mkmh97.sh run=run51 tstop=50 nbody=1000000 m=0.01 step=5
     ./mkmh97.sh run=run52 tstop=50 nbody=1000000 m=0.02 step=5
     ./mkmh97.sh run=run53 tstop=50 nbody=1000000 m=0.04 step=5
     ./mkmh97.sh run=run54 tstop=50 nbody=1000000 m=0.06 step=5
     ./mkmh97.sh run=run55 tstop=50 nbody=1000000 m=0.08 step=5
     ./mkmh97.sh run=run56 tstop=50 nbody=1000000 m=0.10 step=5

     ./mkmh97.sh run=run53 tstop=50 nbody=1000000 m=0.04 step=5 r0=5 v0=0

## Using images to measure mass using ds9 regions

A snapshot that was converted to an image using **snapgrid** stores it's data in mass per square pixel. If the image
was converted to FITS, the default units are degrees.  In **ds9** you can then select a region around the source
of interest (*Edit -> Region*) after which (*Region -> Get Information -> Analysis -> Statistics*) will give the sum
of the values in that region.   For a given pixel size **p** in arcsec (which ds9 reports) the total mass in the region is then
** sum * (p/3600)^2**

## orbits vs. orbits

      mkorbit p2.orb -10 0 0 0 0 0 potname=plummer potpars=0,1,3*pi/16
      orbint p2.orb p2a.orb 50000 0.01 ndiag=10
      orblist p2a.orb  | tabplot - 2 3
      orblist p2a.orb  | tabplot - 2 3 0 600 -12 12

      ./mkmh97.sh run=run202 tstop=500 nbody=10000 m=0.04
      ./mkmh97.sh run=run203 tstop=500 nbody=10000 m=0.04  v0=0
      snapcenter run203.4 run204.4c "weight=i>=10000?-phi*phi*phi:0" one=t
      snapprint run204.4c t,x | tabplot - 1 2 0 500  line=1,1
      snapprint run204.4c t,x > run204.tab

      ./mkmh97.sh run=run205 tstop=500 nbody=100000 m=0.04  v0=0.5
      28823.446 arcsec
      4.37375e-05     left , smaller (2nd splash)
      g2a = 4.603125e-05    = 0.07365
      0.00045568125	      = 0.7290    center (original G2)
      g2b = 0.0001215125    = 0.19442

      14402.926 ar
                fraction
      0.00021045     0.084
      0.00050480001  0.201
      0.001759025    0.704
      0.0025

## Deflection Angle

Instead of using **rp=** one can also use the deflection angle to parametrize the interaction.
For a pointmass approximation this is fine, but for softened models like a Plummer sphere
this is two-valued as for small **rp** the deflection angle goes back to 0.
Here's a simple way to measure the deflection angle, based on the fact that
for unbound encounters the orbit becomes a straight line again, so we measure
that angles as follows:

     rm -rf run2001
     ./mkmh97.sh run=run2001 v0=2 tstop=20 nbody=1000 step=0.5 rp=0.5 seed=123
     tabplot run2001/run2001.4.g2.tab xcoord=0 ycoord=0
     tail -10  run2001/run2001.4.g2.tab  | tabnllsqfit - | txtpar - 'deg(atan(%1)),deg(%2/(1+%1**2))' p0=b=,1,2 p1=b=,1,3
     15.8207 0.113264
     tail -10  run2001/run2001.4.g1.tab  | tabnllsqfit - | txtpar - 'deg(atan(%1)),deg(%2/(1+%1**2))' p0=b=,1,2 p1=b=,1,3     
     15.6923 0.0946659

which reports the deflection angle (and the error) in degrees, based on a linear fit of the last 10 points.   For 1000 bodies
the density and center of mass will not coincide, so interactions do not really have the same **rp**.   Repeating this experiment
100 times gave a mean and dispersion of 14.3 and 1.2 resp.



## Theory

Ignoring mass loss etc. the critical velocity as function of mass (0,1) of the
intrudor should be:

     nemoinp 0.01:1:0.01 | tabmath - - "sqrt((1+1/%1)*(1+%1**1.5)/2)" > criticalv.tab

but this shows a very unexpected trend.

## Versions

In an earlier version (Sep-2022 and before) we didn't use the MH defaults yet. The following snippet
of code produces the same results with the two versions of the script

     mkmh97_old.sh run=run0 seed=123 m=0.02 v0=0.6
     mkmh97_new.sh run=run1 seed=123 m=0.02 v0=0.6 r0=10 kmax=6 eps=0.05 nbody=1000

where the new defaults are now r0=16, kmax=7, eps=0.03125 nbody=2048.

For this seed the relevant output parameters are:

     m16=0.17137
     etot=0.000908545

## Issues

- mh97 claim the energy jumps during an interaction are neglection of tidal interaction terms....
  but their code is an N^2 code on GRAPE-3A.
  Checking with this script with gravidy (code=3) also showed the jumps,
  More accurate forces will lower the energy jump, so it lowers when hackcode1_qp is used
  (run100q), or when gyrfalcON is used (run101).
  For gravidy the jump also appears, but depends on the softening. Run104  vs. run106 clearly shows this.
  Example:
  
  ./mkmh97.sh tstop=20 nbody=1000 seed=123 code=0  kmax=6 step=0.05 run=run100
  ./mkmh97.sh tstop=20 nbody=1000 seed=123 code=0q kmax=6 step=0.05 run=run100q  
  ./mkmh97.sh tstop=20 nbody=1000 seed=123 code=1  kmax=6 step=0.05 run=run101 
  ./mkmh97.sh tstop=20 nbody=1000 seed=123 code=3  kmax=6 step=0.05 run=run103
  ./mkmh97.sh tstop=20 nbody=1000 seed=123 code=3  kmax=6 step=0.05 run=run104 eta=0.001
  ./mkmh97.sh tstop=20 nbody=1000 seed=123 code=3  kmax=6 step=0.05 run=run105 eta=0.001 eps=0.01
  ./mkmh97.sh tstop=20 nbody=1000 seed=123 code=3  kmax=6 step=0.05 run=run106 eps=0.01  

  ./mkmh97.sh tstop=20 nbody=1000 seed=123 code=4  kmax=6 step=0.05 run=run107
  ./mkmh97.sh tstop=20 nbody=1000 seed=123 code=4  kmax=7 step=0.05 run=run108
  ./mkmh97.sh tstop=20 nbody=1000 seed=123 code=4  kmax=7 step=0.05 run=run109 eps=0.01


## Operations

The **save_vars** variable in the script lists the variables that are saved in the **nemopars.rc** file, so they original
values as set on the command line (or the defaults in the script) can be used on subsequent runs.
It is not uncommon to add new variables, which are then not in the nemopars.rc file. Instead the default value would be
seen. This is a problem, and the only solution would be to add those variables manually to the nemopars.rc file, e.g.

      echo zm=1 >> run0e/nemopars.rc

We can currently run 4 types of interactions between G1 (fixed mass=1) and G2 (mass=m):

1. TT72 style:

2. MH97: fully self-consistent

3. MH97 fixed=1 :   G1 is an analytical potential, but G2 is self-consistent, 

4. HQ88 style: MH97 fixed=1 zm=0 :   G1 is an analytical potential, G2 is a point mass with test particles

we cannot run the JSPAM style interactions, this requires knowning the orbit of two smooth potentials under
their own gravity. Also JSPAM has an option to add tidal effect and 


## Simplified setup

Extracting an example of an interaction:

mkplummer run.1 10000
mkplummer - 10000 | snapscale - run.2 mscale=0.05 rscale=0.05**0.5 vscale=0.05**0.25
snapstack run.1 run.2 run.3 deltar=16,0,0 deltav=-0.5,0,0
gyrfalcON run.3 run.4 eps=0.05 kmax=7 step=1 tstop=100 give=mxvap  > run.4.log
