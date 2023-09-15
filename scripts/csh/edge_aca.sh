#! /usr/bin/env bash
#
#  edge_aca.sh   simulate an EDGE ACA cube at CO(2-1)
#                Derived from mkvelfie3 and mkgalcube
#
#
_script=edge_aca.sh
_version=8-sep-2023
_pars=nemopars.rc
_date=$(date +%Y-%m-%dT%H:%M:%S)

#--HELP
#
#  edge_aca.sh   simulate an EDGE ACA cube at 230.538 GHz CO(2-1)
#
run=model1             # identification, and basename for all files       #> ENTRY    
nbody=1000000          # number of bodies per model                       #> RADIO 1000,10000,100000,1000000
r0=10                  # turnover radius (arcsec)                         #> SCALE 1:100:1
v0=200                 # peak velocity (km/s)                             #> SCALE 50:400:10
r1=0.1                 # central unresolved bulge, bar or black hole
v1=300                 # representative rotation speed at r1
re=20                  # exponential scalelength of disk  (arcsec)        #> SCALE 1:100:1
rmax=60                # edge of disk  (arcsec)                           #> SCALE 1:80:1

pa=90                  # PA of receding side disk (E through N)           #> SCALE 0:360:1               
inc=60                 # INC of disk                                      #> SCALE 0:90:1

beam=5                 # FWHM of spatial smoothing beam (arcsec)          #> SCALE 1:20:1
vbeam=5                # FWHM of spectral smoothing beam   (arcsec)       #> SCALE 1:40:2
range=80               # gridding from -range:range  (arcsec)                       
vrange=320             # velocity gridding -vrange:vrange (km/s)
nsize=160              # number of spatial pixels (px=py=2*range/nx)
nvel=120               # number of spectral pixels

z0=0                   # scaleheight [buggy]                              #> SCALE 0:10:0.5

seed=0                 # random seed
sigma=0                # random motion in plane  (km/s)                   #> SCALE 0:20:0.1

noise=0.1              # add optional noise to cube                       #> ENTRY
clip=0.1               # clipping level for cube                          #> ENTRY

refmap=0               # reference map for WCS and masking                #> IFILE 0

show=1                 # display some results (ds9, plots)                #> RADIO 0,1
#--HELP


if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    set +x
    awk 'BEGIN{s=0} {if ($1=="#--HELP") s=1-s;  else if(s) print $0; }' $0
    exit 0
fi

#             simple keyword=value command line parser for bash
for arg in "$@"; do
  export "$arg"
done

#  derive some parameters that appear common or logically belong together

grid_pars="xrange=-${range}:${range} yrange=-${range}:${range} nx=$nsize ny=$nsize"
cell=`nemoinp "2*$range/$nsize*60"`
cen=`nemoinp $nsize/2-0.5`
restfreq=230.53800     # CO(2-1) in GHz

#  keep a log, incase we call this routine multiple times
echo `date` :: $* >> $run.history

# ================================================================================ START

rm -f $run.* >& /dev/null


if [ $v1 = 0 ]; then
  echo "Creating homogeneous disk with $nbody particles times"
  mkdisk out=- nbody=$nbody seed=$seed z0=$z0,$z0 \
       potname=rotcur0 potpars=0,$v0,$r0 mass=1 sign=-1 frac=$sigma abs=t rmax=$rmax |\
    snapmass - - "mass=exp(-r/$re)" |\
    snapscale - - mscale=$nbody |\
    snaprotate - $run.20 "$inc,$pa" yz
else
  m1=`nemoinp "$r1*($v1/0.62)**2"`
  echo "Creating homogeneous disk with $nbody particles times and a nuclear component m1=$m1"
  rotcurves  name1=rotcur0 pars1=0,$v0,$r0  name2=plummer pars2=0,$m1,$r1 tab=t radii=0:${rmax}:0.01 |\
      tabmath - - %1,1,%2,$sigma all > $run.tab
  tabplot $run.tab 1 3 yapp=1/xs

  mktabdisk $run.tab - nbody=$nbody seed=$seed rmax=$rmax sign=-1 |\
    snapmass - - "mass=exp(-r/$re)" |\
    snapscale - - mscale=$nbody |\
    snaprotate - $run.20 "$inc,$pa" yz
fi

echo "Creating velocity moments"
snapgrid $run.20 $run.21 $grid_pars moment=0 evar=m
snapgrid $run.20 $run.22 $grid_pars moment=1 evar=m
snapgrid $run.20 $run.23 $grid_pars moment=2 evar=m
ccdsmooth $run.21 $run.21c $beam
ccdsmooth $run.22 $run.22c $beam
ccdsmooth $run.23 $run.23c $beam
ccdmath $run.21c,$run.22c,$run.23c $run.20d %1
ccdmath $run.21c,$run.22c,$run.23c $run.20v "%2/%1"
ccdmath $run.21c,$run.22c,$run.23c $run.20s "sqrt(%3/%1-%2*%2/(%1*%1))"
#  clipping
ccdmath $run.20d,$run.21 $run.21d "ifgt(%2,0,%1,0)"
ccdmath $run.20v,$run.21 $run.21v "ifgt(%2,0,%1,0)"
ccdmath $run.20s,$run.21 $run.21s "ifgt(%2,0,%1,0)"


echo "Creating the beam"
mkplummer - 1 | snapgrid - $run.p1 $grid_pars
ccdsmooth $run.p1 $run.beam $beam dir=xy

echo "Creating a velocity field - method 2"
snapgrid $run.20 $run.30 $grid_pars \
    zrange=-${vrange}:${vrange} nz=$nvel mean=f evar=m
ccdstat $run.30
if [ $vbeam = 0 ]; then
  ccdmath $run.30 $run.31 "%1+rang(0,$noise)"
else
  ccdmath $run.30 - "%1+rang(0,$noise)" | ccdsmooth - $run.31 $vbeam dir=z
fi
ccdsmooth $run.31 $run.32 $beam dir=xy
ccdmom $run.32 $run.33d axis=3 mom=0 clip=$clip
ccdmom $run.32 $run.33v axis=3 mom=1 clip=$clip rngmsk=true
ccdmom $run.32 $run.33s axis=3 mom=2 clip=$clip

# see if we need WCS from a reference map

if [ -e $refmap ]; then
    ra=$(fitshead $refmap | grep CRVAL1 | awk '{print $3}')
   dec=$(fitshead $refmap | grep CRVAL2 | awk '{print $3}')
  vref=$(fitshead $refmap | grep CRVAL3 | awk '{print $3}')
  echo "REFMAP: $ra $dec $vref"
  # @todo figure out of m/s,km/s or freq
  vscale=1000
  vlsr=$(nemoinp $vref/$vscale)
  crval=$ra,$dec,$vlsr 
else
  vlsr=0    
  crval=180,0,0
fi
     
# single dish profile
ccdmom $run.32 - axis=1 mom=0 |\
    ccdmom - - axis=2 mom=0 |\
    ccdprint - x= y= z= label=z newline=t |\
    tabcomment - - punct=f delete=t > $run.spec

# export for barolo or so, in decent units (could also use ccdsky)
# this way the input spatial scale is in arcsec and km/s
ccdfits $run.32 $run.fits radecvel=t scale=1/3600.0,1/3600.0,1.0 crval=$crval restfreq=$restfreq

if [ $show = 1 ]; then
    xpaset -p ds9 frame frameno 1
    nds9 $run.33d
    xpaset -p ds9 frame frameno 2
    nds9 $run.33v
    xpaset -p ds9 frame frameno 3
    nds9 $run.33s
    xpaset -p ds9 frame frameno 4
    nds9 $run.fits

    tabplot $run.spec line=1,1  headline="Spectrum around VLSR=$vlsr" yapp=2/xs
fi


echo "CCDSTAT"
ccdstat $run.32
echo "QAC_STATS              mean        rms         min        max"
ccdstat $run.32 qac=t          label='cube/full  '
ccdstat $run.32 qac=t robust=t label='cube/robust'
echo -n "Approx. Signal/Noise:  "
ccdstat $run.32 qac=t robust=t | txtpar - %1/%2 p0=1,6 p1=1,4


