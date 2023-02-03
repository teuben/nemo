#!/bin/bash
#
#    Set up a collision between two spherical (Plummer) galaxies, similar to Makino & Hut (1997)
#
# history: 30-jan-2023   recoded old snippets of code (following mkmh97.sh) in bash

_version=3-feb-2023

#--HELP
#            parameters for the integration
run=run0        # directory and basename of the files belonging to this simulation 
nbody=2048      # number of bodies in one model
step=1          # step in time when to dump snapshots
v0=1.0          # initial impact/circular speed
rp=0.0          # impact offset radius
r0=16.0         # initial offset position for v0 > 0   [note v0 > 2/sqrt(r0)]
eps=0.03125     # gravitational softening
kmax=7          # integration timestep is 1/2**kmax
seed=0          # random seed (use seed=123 for the benchmark)
tstop=50        # stop time of the integration (or analysis time when doing a re-run)
box=32          # spatial box size for plotting and CCD frames
vbox=2          # velocity box size
npixel=128      # number of pixels in xy CCD frame
tplot=0,5,10,15,20,25,30,40,50     # times to plot the evolution
yapp=/xs                           # /xs  [yapp_pgplot preferred]
debug=1                            # 1=set -x,-e,-u   0=nothing
#
#--HELP

#  display the '#--HELP' portion(s) if requested
if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    set +x
    awk 'BEGIN{s=0} {if ($1=="#--HELP") s=1-s;  else if(s) print $0; }' $0
    exit 0
fi
#  simple keyword=value command line parser for bash
for arg in "$@"; do
    eval "$arg"
done
#   add various helpful things for bash
if [ $debug == 1 ]; then
    set -x
    set -e
    set -u
fi
#   make sure the directory does not exist yet
if [ -d $run ]; then
    echo Run directory run=$run already exists, cannot overwrite
    exit 0
fi

# do all the work in a run directory
mkdir $run
cd $run

# make two random plummer spheres in virial (N-body) units 
mkplummer $run.1 $nbody seed=$seed
mkplummer $run.2 $nbody seed=$seed

# (near) head-on collision
snapstack $run.1 $run.2 $run.3 deltar=$r0,$rp,0 deltav=-$v0,0,0  zerocm=t

# integrate, and keep a table of time and energy
hackcode1 $run.3 $run.4 eps=$eps freq=2**$kmax freqout=1/$step fcells=2 tstop=$tstop options=mass,phase,phi,acc > $run.4.log    
snapdiagplot $run.4 tab=$run.4.etot

# analysis: energy conservation
tabplot $run.4.etot 1 2  yapp=1$yapp
tabhist $run.4.etot 2 yapp=/null

# plot total and each galaxy by itself
r=-$box:$box
snapplot $run.4 xrange=$r yrange=$r times=$tstop                   yapp=2$yapp
snapplot $run.4 xrange=$r yrange=$r times=$tstop visib="i<$nbody"  yapp=3$yapp
snapplot $run.4 xrange=$r yrange=$r times=$tstop visib="i>=$nbody" yapp=4$yapp

# plot color by binding energy. In yapp_pgplot 0.1 is red for unbound, 0.2 is green for bound
snapplot  $run.4 color='phi+0.5*(vx*vx+vy*vy+vz*vz) > 0 ? 0.1 : 0.2' psize=0.05 xrange=$r yrange=$r           times=$tstop yapp=5$yapp
snapplot3 $run.4 color='phi+0.5*(vx*vx+vy*vy+vz*vz) > 0 ? 0.1 : 0.2' psize=0.05 xrange=$r yrange=$r zrange=$r times=$tstop yapp=6$yapp
snapplot  $run.4 color='phi+0.5*(vx*vx+vy*vy+vz*vz) > 0 ? 0.1 : 0.2' psize=0.05 xvar=r yvar=vr xrange=0:$box yrange=-$vbox:$vbox times=$tplot nxy=3,3 yapp=7$yapp

# make a new snapshot putting the bound particles first, and the unbound particles with half mass
# to make it easier for glnemo2 to visualize them.
snaptrim run0.4 - times=$tstop | unbind - run0.b
snaptrim run0.4 - times=$tstop | unbind - - bind=f | snapscale - run0.u mscale=0.5
snapstack run0.b run0.u run0.bu
select=$(snapmstat run0.bu | tail -1)
echo "Do something like: glnemo2 run0.bu $select"
