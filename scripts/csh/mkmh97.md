#    Interactions between equal mass Plummer spheres

Based on the paper by Makino & Hut (1997) - https://arxiv.org/pdf/2205.04623.pdf

## Input parameters

All input parameters are of the form *keyword=value*, like most NEMO programs, except
there is no syntax or validity of the keyword checking!

As in NEMO, units are virial units.

1. **run**:  an identifying name used in the filenames, e.g. run1. If you re-specify the same
same, that run is removed!

2. **nbody**: number of bodies per plummer sphere [1000]

3. **tstop**: stopping time. Should be several times r0/v0   [50]

4. **step**: step time when snapshots are stored. 1 is probably ok, for movies you probably need 0.1 [1.0]

5. **v0**: initial impact speed (todo:  should define this at infinity) .
Use positive values for a collision. Negative values are reserved for a circular orbit, although
with the right choices of **v0** and **r0** a circular orbit can also be achieved.  [1.0]

6. **rp**: initial impact parameter. 0 means a head-on collision [0]

7. **r0**: initial impact distance [10]

8. **eps**: gravitational softening. There are some codes in NEMO that allow negative softening,
   in which case a Post-Newtonian approximation is used. This is outside the real of this study. [0.05]


## Some plotting

initial conditions:

      snapplot run1.3 xrange=-16:16 yrange=-16:16

      snapplot run1.3 xrange=-16:16 yrange=-2:2 yvar=vx

      snapplot run1.4 xrange=-16:16 yrange=-2:2 yvar=vx trak=t

some evolution plots

      tabplot run1.4.log
      tabhist run1.4.log 2
      tabplot run1.4.log xrange=-4:4 yrange=-4:4
      tabplot run1.4.log xrange=0:4 yrange=-2:2 xvar=r xvar=jz
      


## Circular orbits?

For **v0<0** we can set up the two sytems in a circular orbit by launching from
(r0,0,0) and (0,v0,0), i.e. in clock wise motion. You will need to figure out the
right v0 for a circular orbit. E.g.

      ./mkmh97.sh v0=-1 r0=1 nbody=1
      snapplot run0.4 trak=t

is clearly not the right velocity.

## TODO

* define v0 as being from infinity? (cf. Makino & Hut paper)

* allow circular orbit using negative **v0**.   The initial offset **r0** is then used
as the diameter of the circular orbit. This could be a fun way to study dynamical
friction (e.g. Bontekoe & v Albada, White, ...)

* allow softened potentials (also useful for nbody=1 experiments)


