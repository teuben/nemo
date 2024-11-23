#! /usr/bin/env bash
#
#  edge_gbt.sh   simulate an EDGE GBT spectrum, derived from edge_aca.sh
#
#  bench:    /usr/bin/time ./edge_gbt.sh logn=7
#  7.65user 5.29system 0:12.88elapsed 100%CPU    orig
#  5.15user 3.21system 0:08.37elapsed 99%CPU     skip rotate, no moments
#  3.74user 2.21system 0:05.98elapsed 99%CPU     skip snapsort (mkdisk was sorted)
#  3.11user 2.30system 0:05.40elapsed 100%CPU    faster mkdisk

_script=edge_gbt.sh
_version=23-nov-2024
_pars=nemopars.rc
_date=$(date +%Y-%m-%dT%H:%M:%S)

#--HELP
#
#  edge_gbt.sh   simulate an EDGE GBT spectrum
#
run=model1             # identification, and basename for all files       #> ENTRY    
logn=6                 # log of number of bodies per model                #> SCALE 1:8:0.2
r0=10                  # characteristic radius (arcsec)                   #> SCALE 1:100:1
v0=200                 # characteristic velocity (km/s)                   #> SCALE 10:500:5
m0=0.5                 # brandt power  (0=flat)                           #> SCALE 0:10:0.1
r1=0.1                 # central unresolved bulge, bar or black hole
v1=0                   # representative rotation speed at r1
re=20                  # exponential scalelength of disk  (arcsec)        #> SCALE 1:100:1
n=-1                   # n<1 exp disk; n>=0 PLEC disk                     #> SCALE -1:10:0.1
rm=-1                  # if > 0: use it has R_mol^2 - AMT model           #> SCALE -1:50:0.1
rmax=60                # edge of disk  (arcsec)                           #> SCALE 1:80:1
inc=60                 # INC of disk                                      #> SCALE 0:90:1
z0=0                   # scaleheight [@todo buggy]                        #> SCALE 0:10:0.5

vbeam=5                # FWHM of spectral smoothing beam   (arcsec)       #> SCALE 1:40:2
vrange=400             # velocity gridding -vrange:vrange (km/s)
nvel=200               # number of spectral pixels

seed=0                 # random seed
sigma=0                # random motion in plane  (km/s)                   #> SCALE 0:20:0.1

noise=0                # add optional noise to cube                       #> ENTRY
vlsr=0                 # optional VLSR if non-zero                        #> ENTRY

debug=-1               # add debugging                                    #> SCALE -1:9:1

plot=profile,rotcur    # show which plots can be made                     #> CHECK rotcur,profile,density

#--HELP

#  basic help if --help or -h is given
if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    set +x
    awk 'BEGIN{s=0} {if ($1=="#--HELP") s=1-s;  else if(s) print $0; }' $0
    echo "# Script version: $_script $_version"
    exit 0
fi

#             simple keyword=value command line parser for bash
for arg in "$@"; do
  export "$arg"
done

#  handle debugging
debug=$(nemoinp $debug format=%d)
if [ $debug -gt 0 ]; then
    set -x
fi
export DEBUG=$debug
echo DEBUG=$DEBUG


#  fixed now compared to edge_aca.sh
range=$rmax          # gridding size
nsize=1              # number of spatial pixels (px=py=2*range/nx)


#  derive some parameters that appear common or logically belong together

grid_pars="xrange=-${range}:${range} yrange=-${range}:${range} nx=$nsize ny=$nsize"
cell=`nemoinp "2*$range/$nsize*60"`
cen=`nemoinp $nsize/2-0.5`
restfreq=230.53800     # CO(2-1) in GHz
nbody=`nemoinp "10**$logn" format=%d`
sininc=`nemoinp "sind($inc)"`

function nemo_stamp {
    if [ $debug -ge 0 ]; then
	echo "$(date +%H:%M:%S.%N) $*"
    fi
}


# ================================================================================ START

#  Announce:
echo "$0 version $_version"
nemo_stamp start

#  Clear old model
rm -f $run.* >& /dev/null

#  keep a log
echo "`date` :: $*" > $run.history


mmode=$(nemoinp "ifge($n,0,1,0)")
if [ $mmode = 0 ]; then
    mass="exp(-r/$re)"
else
    mass="pow(r/$re*exp(1-r/$re),$n)"
fi
mmode=$(nemoinp "ifge($rm,0,2,$mmode)")
if [ $mmode = 2 ]; then
    mass="exp(-r/$re)/(1+$rm*exp(-1.6*r/$re))"
fi
echo MMODE=$mmode MASS=$mass
echo NBODY=$nbody
nemo_stamp begin
if [ $m0 != 0 ]; then
  echo "Creating brandt disk with $nbody particles."
  mkdisk out=- nbody=$nbody seed=$seed z0=$z0,$z0 \
       potname=brandt potpars=0,$v0,$r0,$m0 mass=1 sign=-1 frac=$sigma abs=t rmax=$rmax |\
      snapmass - $run.20 "mass=$mass" norm=1
  nemo_stamp mkdisk
  rotcurves  name1=brandt pars1=0,$v0,$r0,$m0  tab=t radii=0:${rmax}:0.01 plot=f |\
      tabmath - - "%1,1,%2,$sigma,exp(-%1/$re)" all  > $run.tab
  nemo_stamp rotcurves
  vmax=$(tabstat $run.tab 3 | txtpar - p0=max:,1,2)
  echo "Max rotcur: vmax=$vmax m0=$m0"
  nemo_stamp tabstat
elif [ $v1 = 0 ]; then
  echo "Creating homogeneous disk with $nbody particles.  (need optimizing test)"
  mkdisk out=- nbody=$nbody seed=$seed z0=$z0,$z0 \
       potname=rotcur0 potpars=0,$v0,$r0 mass=1 sign=-1 frac=$sigma abs=t rmax=$rmax |\
    snapmass - - "mass=$mass" |\
    snapscale - $run.20 mscale=$nbody
  rotcurves  name1=rotcur0 pars1=0,$v0,$r0  tab=t radii=0:${rmax}:0.01 plot=f |\
      tabmath - - %1,1,%2,$sigma all > $run.tab
else
  m1=`nemoinp "$r1*($v1/0.62)**2"`
  echo "Creating homogeneous disk with $nbody particles and a nuclear component m1=$m1 (needs optimizing test)"
  rotcurves  name1=rotcur0 pars1=0,$v0,$r0  name2=plummer pars2=0,$m1,$r1 tab=t radii=0:${rmax}:0.01 |\
      tabmath - - %1,1,%2,$sigma all > $run.tab

  mktabdisk $run.tab - nbody=$nbody seed=$seed rmax=$rmax sign=-1 |\
    snapmass - - "mass=$mass" |\
    snapscale - $run.20 mscale=$nbody
fi

#snapsort $run.20 - r help=c | snapshell - radii=0:${rmax}:0.01 pvar=m help=c > $run.shell
snapshell $run.20 radii=0:${rmax}:0.01 pvar=m help=c > $run.shell
nemo_stamp snapshell

echo "Creating a velocity field - method 2"
snapgrid $run.20 $run.30 $grid_pars \
	 zrange=-${vrange}:${vrange} nz=$nvel mean=f evar=m zvar=vy*$sininc
nemo_stamp snapgrid
ccdstat $run.30
nemo_stamp ccdstat

if [ $vbeam = 0 ]; then
  ccdmath $run.30 $run.31 "%1+rang(0,$noise)"
else
  ccdmath $run.30 - "%1+rang(0,$noise)" | ccdsmooth - $run.31 $vbeam dir=z
fi
nemo_stamp ccdmath

# single dish profile
if [ 1 = 1 ]; then
ccdmom $run.31 - axis=1 mom=0 |\
    ccdmom - - axis=2 mom=0 |\
    ccdprint - x= y= z= label=z newline=t |\
    tabcomment - - punct=f delete=t > $run.spec
else
    # faster, but normalization is off now
    ccdprint $run.31 x= y= z= label=z newline=t |\
    tabcomment - - punct=f delete=t > $run.spec
fi
nemo_stamp ccdprint

echo -n "Total integral over spectrum: "
tabint $run.spec
nemo_stamp tabint

# export for other programs, in decent units
# this way the input spatial scale is in arcsec and km/s
# SKIP for production
# ccdfits $run.31 $run.fits radecvel=t scale=1/3600.0,1/3600.0,1.0 crval=$crval restfreq=$restfreq


echo "PLOT: $plot"
if [[ "$plot" == *"rotcur"* ]]; then
    tabplot $run.tab 1 3 headline="Rotation Curve" yapp=1/xs
fi
if [[ "$plot" == *"profile"* ]]; then    
    tabplot $run.spec line=1,1  headline="Spectrum around VLSR=$vlsr" yapp=2/xs
fi
if [[ "$plot" == *"density"* ]]; then    
    tabplot $run.shell 1 4 line=1,1  ymin=0 headline="Surface Density" yapp=3/xs
fi

nemo_stamp end
