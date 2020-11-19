#!/bin/bash
#
#    Example of setting up a Toomtre & Toomre type interaction, and some
#    further analysis
#
#    Many todo's:
#    - add example to get snapfit to work (self-consistency check)
#    - recode the parameters to match the 16 input parameters that indentikit advocates
#
#

set -x

#            parameters for the script that can be overriden via the commandline
run=run1
box=2
nbody=100
nx=64

rscale=0.25
vscale=0.5
theta=60
axis=y

inc=0
pa=0

#             simple keyword=value command line parser for bash
for arg in $*; do\
  export $arg
done

#             delete this old run
rm -f $run.*


#             derived parameters, do not change anything below
ny=$nx
nz=$nx
range=-$box:$box

echo "1:4 mass ratio sequence"
#   the radii must depend on $nbody such that particles have near neighbors in angle and radius about equal
mktt72 - nbody=$nbody radius=0.025:.25:"0.025*2*pi/$nbody"    central=t eps=0.2 | snapmerge - - | snapmass - $run.1 mass='exp(-r/0.5)' aux=t
mktt72 - nbody=$nbody radius=0.1:1:"0.1*2*pi/$nbody" mass=4.0 central=t eps=0.2 | snapmerge - - | snapmass - $run.2 mass='exp(-r/0.5)' aux=t
snaprotate $run.1 $run.1r theta=$inc,$pa order=yz
snapstack $run.1r $run.2 $run.3 5,2,0 -1,0,0  zerocm=t


/usr/bin/time hackcode1 $run.3 $run.4 eps=0.2 freqout=10 freq=100 tstop=10 > $run.4.log



# pick a time
# scale it by something so we can check if snapfit finds it
# rotate i a bit so we can check if snapfit finds it
# grid and smooth into a cube
# make a moment0 (total gas) and moment1 (mean velocity) map for reference
 
snaptrim $run.4 - times=5 | snapscale - - rscale=$rscale vscale=$vscale |\
    snaprotate - - $theta $axis | snapgrid - $run.5.cube xrange=$range yrange=$range zrange=$range nx=$nx ny=$ny nz=$nz evar=1 svar=0.05 szvar=0.05
ccdmom $run.5.cube $run.5.mom0 mom=0
ccdmom $run.5.cube $run.5.mom1 mom=1
ccdplot $run.5.mom0 power=0.4 yapp=1/xs
nds9 $run.5.cube

exit 0


# very time consuming program, still has bugs/problems
snapfit $run.4 $run.5.cube theta1=-50:50:3 theta2=-170:170:3   contour=$run.cnt  iter=2    > $run.snapfit.log

# look if over time there is a log10(sum) that looks like a minimum
grep ^Min $run.snapfit.log | tabplot - 9 7 line=1,1  yapp=2/xs
grep ^Min $run.snapfit.log | tabplot - 9 11,12,13,14 line=1,1 color=2,3,4,5 yapp=3/xs

