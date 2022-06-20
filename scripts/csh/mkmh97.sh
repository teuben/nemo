#!/bin/bash
#
#    Example of setting up a Makino & Hut 1997 simulation
#    See also: https://ui.adsabs.harvard.edu/abs/1997ApJ...481...83M/abstract
#
#    @todo
#    - allow circular orbit
#    - recode such that v0 is speed at infinity
#
# version: 12-may-2022   initial version with just (near) head-on collision
#          15-may-2022   added option to make (near) circular orbit
#          20-jun-2022   add potentials & acc

set -x

#            parameters for the script that can be overriden via the commandline
run=run0        # basename of the files belonging to this simulation
nbody=1000      # number of bodies in one model
tstop=50        # stop time of the integration
step=1          # step when to dump snapshots
v0=1.0          # initial impact/circular speed
rp=0.0          # impact offset radius (not used for v
r0=10.0         # initial offset position for v0 > 0
eps=0.05        # softening
kmax=8          # integration timestep is 1/2**kmax
hack=1          # hackcode? else use the faster gyrfalcON code

#             simple keyword=value command line parser for bash
for arg in $*; do
  export $arg
done

#             delete the old run
rm -f $run.*

# make two random plummer spheres in virial units and stack them
mkplummer $run.1 $nbody
mkplummer $run.2 $nbody
if [ $(nemoinp "ifgt($v0,0,1,0)") = 1 ]; then
    # (near) head-on collision
    snapstack $run.1 $run.2 $run.3 deltar=$r0,$rp,0 deltav=-$v0,0,0  zerocm=t
else
    # (near) circular orbit
    snapstack $run.1 $run.2 $run.3 deltar=$r0,0,0   deltav=0,$v0,0   zerocm=t    
fi

# integrate (hackcode1 is slower for large Nbody systems)
if [ $hack = 1 ]; then
    hackcode1 $run.3 $run.4 eps=$eps freq=2**$kmax freqout=1/$step fcells=2 tstop=$tstop options=mass,phase,phi,acc > $run.4.log
    snapdiagplot $run.4 tab=$run.4.etot
else
    gyrfalcON $run.3 $run.4 eps=$eps kmax=$kmax step=$step tstop=$tstop give=mxvap > $run.4.log
    tabcols $run.4.log 1,2 > $run.4.etot
fi


