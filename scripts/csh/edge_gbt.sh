#! /usr/bin/env bash
#
#  edge_gbt.sh   simulate an EDGE GBT spectrum, derived from edge_aca.sh
#
#
_script=edge_gbt.sh
_version=19-jun-2024
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


#  fixed now compared to edge_aca.sh
range=$rmax          # gridding size
nsize=1              # number of spatial pixels (px=py=2*range/nx)
pa=0                 # PA of receding side disk (E through N)           


#  derive some parameters that appear common or logically belong together

grid_pars="xrange=-${range}:${range} yrange=-${range}:${range} nx=$nsize ny=$nsize"
cell=`nemoinp "2*$range/$nsize*60"`
cen=`nemoinp $nsize/2-0.5`
restfreq=230.53800     # CO(2-1) in GHz
nbody=`nemoinp 10**$logn format=%d`


# ================================================================================ START

#  Announce:
echo "$0 version $_version"

#  Clear old model
rm -f $run.* >& /dev/null

#  keep a log
echo "`date` :: $*" > $run.history

if [ $m0 != 0 ]; then
  echo "Creating brandt disk with $nbody particles times"    
  mkdisk out=- nbody=$nbody seed=$seed z0=$z0,$z0 \
       potname=brandt potpars=0,$v0,$r0,$m0 mass=1 sign=-1 frac=$sigma abs=t rmax=$rmax |\
    snapmass - - "mass=exp(-r/$re)" |\
    snapscale - - mscale=$nbody |\
    snaprotate - $run.20 "$inc,$pa" yz
  rotcurves  name1=brandt pars1=0,$v0,$r0,$m0  tab=t radii=0:${rmax}:0.01 plot=f |\
      tabmath - - "%1,1,%2,$sigma,exp(-%1/$re)" all  > $run.tab
  vmax=$(tabstat $run.tab 3 | txtpar - p0=max:,1,2)
  echo "Max rotcur: vmax=$vmax m0=$m0"
elif [ $v1 = 0 ]; then
  echo "Creating homogeneous disk with $nbody particles times"
  mkdisk out=- nbody=$nbody seed=$seed z0=$z0,$z0 \
       potname=rotcur0 potpars=0,$v0,$r0 mass=1 sign=-1 frac=$sigma abs=t rmax=$rmax |\
    snapmass - - "mass=exp(-r/$re)" |\
    snapscale - - mscale=$nbody |\
    snaprotate - $run.20 "$inc,$pa" yz
  rotcurves  name1=rotcur0 pars1=0,$v0,$r0  tab=t radii=0:${rmax}:0.01 plot=f |\
      tabmath - - %1,1,%2,$sigma all > $run.tab
else
  m1=`nemoinp "$r1*($v1/0.62)**2"`
  echo "Creating homogeneous disk with $nbody particles times and a nuclear component m1=$m1"
  rotcurves  name1=rotcur0 pars1=0,$v0,$r0  name2=plummer pars2=0,$m1,$r1 tab=t radii=0:${rmax}:0.01 |\
      tabmath - - %1,1,%2,$sigma all > $run.tab

  mktabdisk $run.tab - nbody=$nbody seed=$seed rmax=$rmax sign=-1 |\
    snapmass - - "mass=exp(-r/$re)" |\
    snapscale - - mscale=$nbody |\
    snaprotate - $run.20 "$inc,$pa" yz
fi

echo "Creating velocity moments"
snapgrid $run.20 $run.21 $grid_pars moment=0 evar=m
snapgrid $run.20 $run.22 $grid_pars moment=1 evar=m
snapgrid $run.20 $run.23 $grid_pars moment=2 evar=m
# skip smoothing
ccdmath $run.21,$run.22,$run.23 $run.20d %1
ccdmath $run.21,$run.22,$run.23 $run.20v "%2/%1"
ccdmath $run.21,$run.22,$run.23 $run.20s "sqrt(%3/%1-%2*%2/(%1*%1))"
#  clipping
ccdmath $run.20d,$run.21 $run.21d "ifgt(%2,0,%1,0)"
ccdmath $run.20v,$run.21 $run.21v "ifgt(%2,0,%1,0)"
ccdmath $run.20s,$run.21 $run.21s "ifgt(%2,0,%1,0)"


echo "Creating a velocity field - method 2"
snapgrid $run.20 - $grid_pars \
	 zrange=-${vrange}:${vrange} nz=$nvel mean=f evar=m |\
    ccdflip - $run.30 flip=$flip wcs=t
ccdstat $run.30
if [ $vbeam = 0 ]; then
  ccdmath $run.30 $run.31 "%1+rang(0,$noise)"
else
  ccdmath $run.30 - "%1+rang(0,$noise)" | ccdsmooth - $run.31 $vbeam dir=z
fi

# single dish profile
ccdmom $run.31 - axis=1 mom=0 |\
    ccdmom - - axis=2 mom=0 |\
    ccdprint - x= y= z= label=z newline=t |\
    tabcomment - - punct=f delete=t > $run.spec

# export for other programs, in decent units
# this way the input spatial scale is in arcsec and km/s
ccdfits $run.31 $run.fits radecvel=t scale=1/3600.0,1/3600.0,1.0 crval=$crval restfreq=$restfreq


echo "PLOT: $plot"
if [[ "$plot" == *"rotcur"* ]]; then
    tabplot $run.tab 1 3 headline="Rotation Curve" yapp=1/xs
fi
if [[ "$plot" == *"profile"* ]]; then    
    tabplot $run.spec line=1,1  headline="Spectrum around VLSR=$vlsr" yapp=2/xs
fi
if [[ "$plot" == *"density"* ]]; then    
    tabplot $run.tab 1 5 line=1,1  ymin=0 headline="Surface Density" yapp=3/xs
fi

