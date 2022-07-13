#!/bin/bash
#
#    Example of setting up a Makino & Hut 1997 simulation
#    See also: https://ui.adsabs.harvard.edu/abs/1997ApJ...481...83M/abstract
#
#    @todo
#    - recode such that v0 is speed at infinity
#
# bench:
#          ./mkmh97.sh code=0 run=run100 nbody=10000 tstop=20 seed=123 
#          ./mkmh97.sh code=1 run=run101 nbody=10000 tstop=20 seed=123
#                     i7-1185G7  Xeon E5-2687W      i9-12900K
#          code=0     5:12 min       10:48            3:31
#          code=1     2:10 min        3:20            1:25
#          code=2     ----            0:16            ----

#
# version: 12-may-2022   initial version with just (near) head-on collision
#          15-may-2022   added option to make (near) circular orbit
#          15-jun-2022   added m= mass ratio parameter, changed hack= -> code= and implemented bonsai
#          17-jun-2022   added seed= and defined a benchmark
#          20-jun-2022   add potentials & acc
#          27-jun-2022   $run is now a directory
#           7-jul-2022   $em=1 as option to have equal mass particles
#          11-jul-2022   able to rerun if the run directory exists

set -x
set -e

#            parameters for the integration
run=run0        # directory and basename of the files belonging to this simulation
nbody=1000      # number of bodies in one model
em=0            # equal mass particles?
step=1          # step when to dump snapshots
v0=1.0          # initial impact/circular speed
rp=0.0          # impact offset radius (not used for v
r0=10.0         # initial offset position for v0 > 0
eps=0.05        # softening
kmax=8          # integration timestep is 1/2**kmax
code=0          # 0=hackcode1 2=gyrfalcON  3=bonsai2
m=1             # mass of second galaxy (mass of first will always be 1)
seed=0          # random seed
#             parameters for the analysis
box=32          # box size for plotting and CCD frames
npixel=256      # number of pixels in CCD frame
tstop=50        # stop time of the integration (or analysis time)

#             simple keyword=value command line parser for bash
for arg in $*; do
  export $arg
done

version=12-jul-2022

#             delete the old run, work within the run directory
if [ -d $run ]; then
    restart=0
else
    restart=1
fi
mkdir -p $run
cd $run
echo "# version=$version" >> mkmk97.rc
echo "$*"                 >> mkmk97.rc

source  mkmk97.rc


if [ $restart = 1 ]; then
    # make two random plummer spheres in virial units and stack them
    # optionally: use fewer particles in the impactor so all particles have the same mass (em=1)
    mkplummer $run.1 $nbody      seed=$seed
    if [ $em = 0 ]; then
	mkplummer -      $nbody      seed=$seed | snapscale - $run.2 mscale="$m" rscale="$m**0.5" vscale="$m**0.25"
    else
	mkplummer -      "$nbody*$m" seed=$seed | snapscale - $run.2 mscale="$m" rscale="$m**0.5" vscale="$m**0.25"
    fi
    if [ $(nemoinp "ifgt($v0,0,1,0)") = 1 ]; then
	# (near) head-on collision
	snapstack $run.1 $run.2 $run.3 deltar=$r0,$rp,0 deltav=-$v0,0,0  zerocm=t
    else
	# (near) circular orbit
	snapstack $run.1 $run.2 $run.3 deltar=$r0,0,0   deltav=0,$v0,0   zerocm=t    
    fi

    # integrator:  0:  hackcode1 is O(NlogN) code
    #              1:  gyrfalcON is O(N)
    #              2:  bonsai2 is O(N) but scales faster for "small" N
    echo "Use:   tail -f $run/$run.4.log     to monitor progress of the integrator"
    if [ $code = 0 ]; then
	hackcode1 $run.3 $run.4 eps=$eps freq=2**$kmax freqout=1/$step fcells=2 tstop=$tstop options=mass,phase,phi,acc > $run.4.log    
	snapdiagplot $run.4 tab=$run.4.etot
    elif [ $code = 1 ]; then
	gyrfalcON $run.3 $run.4 eps=$eps kmax=$kmax step=$step tstop=$tstop give=mxvap > $run.4.log
	tabcols $run.4.log 1,2 > $run.4.etot
    elif [ $code = 2 ]; then
	# See $NEMO/usr/bonsai
	snaptipsy  $run.3 $run.3t
	bonsai2 -i $run.3t --snapname $run. --logfile $run.gpu --dev 0 --snapiter $step -T $tstop -t $(nemoinp 1/2**$kmax)  1>$run.4.log 2>$run.4.err
	grep Etot= run0.4.log  |tabcols - 4,6 > $run.4.etot
	set +x
	for f in $run._*; do
	    echo Processing $f
	    (tipsysnap $f - | csf - - item=SnapShot 1>>$run.4 ) 2>>$run.4.t2s
	done
	echo Wrote final combined snapshot in $run.4
    else
	set +x
	echo Unknown code=$code, valid are:
	echo 1 = hackcode1
	echo 2 = gyrfalcON
	echo 3 = bonsai2
    fi
else
    echo "============================================================================="
    echo "Skipping integration because run=$run already existed"
fi

# now some analysis will follow

bsigma=0.0001
tplot=0,5,10,15,20,25,30,40,50

tabplot $run.4.etot 1 2  yapp=etot.plot.png/png > /dev/null 2>&1
tabhist $run.4.etot 2    yapp=etot.hist.png/png > etot.hist.log 2>&1
snapplot $run.3 xrange=-$box:$box yrange=-$box:$box              yapp=init.plot.png/png
snapplot $run.4 xrange=-$box:$box yrange=-$box:$box times=$tstop yapp=final.plot.png/png
snapgrid $run.3 - xrange=-$box:$box yrange=-$box:$box              nx=$npixel ny=$npixel |\
    ccdmath - - "log(1+%1/$bsigma)" |\
    ccdplot - yapp=init.ccd.png/png headline="Initial Conditions"
snapgrid $run.4 - xrange=-$box:$box yrange=-$box:$box times=$tstop nx=$npixel ny=$npixel |\
    ccdmath - - "log(1+%1/$bsigma)" |\
    ccdplot - yapp=final.ccd.png/png headline="Conditions tstop=$tstop"
snapgrid $run.4 - xrange=-$box:$box yrange=-$box:$box times=$tstop nx=$npixel ny=$npixel evar="i<$nbody?m:0" |\
    ccdmath - - "log(1+%1/$bsigma)" |\
    ccdplot - yapp=final1.ccd.png/png headline="Galaxy-1 at tstop=$tstop"
snapgrid $run.4 - xrange=-$box:$box yrange=-$box:$box times=$tstop nx=$npixel ny=$npixel evar="i>=$nbody?m:0" |\
    ccdmath - - "log(1+%1/$bsigma)" |\
    ccdplot - yapp=final2.ccd.png/png headline="Galaxy-2 at tstop=$tstop"

snapplot  $run.4 xrange=-$box:$box yrange=-$box:$box                   times=$tplot nxy=3,3 yapp=evolution.plot.png/png
snapplot3 $run.4 xrange=-$box:$box yrange=-$box:$box zrange=-$box:$box times=$tstop         yapp=final.3d.plot.png/png
