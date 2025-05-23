#! /usr/bin/env bash
#
#   fits ellipse to an inclined galaxy image
#     31-mar-2022   Created   - Peter Teuben
#
#
#   Some examples of making a disk galaxy
#
# mkexpdisk disk1 10000000  z0=0.025   (these all have r_e=0.25)
# mkexpdisk disk2 10000000  z0=0.05
# mkexpdisk disk3 10000000  z0=0.1
# mkexpdisk disk4 10000000  z0=0.1 mode=2 zmode=3
# mkexphot expmodel.inp - 10000000 | snapscale - disk5 rscale=0.05
#      edit the expmodel.inp with scaleheight to 1.0 (from 0.2)
# mkplummer - 10000000 | snapscale - disk6 rscale=1,1,0.2
#
# GalIC ->
# DICE  ->

# -----------------------------------------------------------------
# Example output :    70 0.18 0.450205 63.2432 65.1965 0.311546
# -----------------------------------------------------------------
#   <file_name> 
#   INC    input known INC
#   Q0     input Q0  (e.g. measured via INC=90),or set to 0
#   CONT   contour level used 
#   Q      measured Q at given contour
#   INC_Q  predicted INC based on thin disk
#   INC_Q0 predicted INC based on Hubble's correction with Q0
#   Q1     predicted Q0 based on INC_Q0 = INC
# -----------------------------------------------------------------
#
# Need: snaprotate snapgrid ccdsmooth ccdplot tablsqfit txtpar


in=disk1   #>  IFILE   in=disk1
inc=70     #>  SCALE   inc=70        0:90:1
mag=2      #>  SCALE   mag=2         0:10:0.05
size=1     #>  SCALE   size=1        0.05:10:0.05
n=128      #>  SCALE   n=128         64:1024:64
sp=2       #>  SCALE   sp=2          0.5:5:0.5
debug=-1   #>  SCALE   debug=-2      -2:9:1
yapp=/xs   #>  ENTRY   yapp=/xs      
q0=0       #>  ENTRY   q0=0

#> HELP in       input snapshot
#> HELP inc      inclination
#> HELP contour  contour
#> HELP size     -size:size
#> HELP n        pixels
#> HELP sp       smoothing in pixels
#> HELP debug    NEMO's $DEBUG
#> HELP q0       optional


for arg in $*; do
    export $arg
done

export DEBUG=$debug
export YAPP=$yapp

s=$(nemoinp 2*$sp*$size/$n)
cnt=$(nemoinp "10**(-$mag/2.5)")

if [ ! -e $in ]; then
    echo File $in does not exist. Example of making one:
    echo mkexpdisk disk1 100000 z0=0.1
    exit 0
fi
    

snaprotate $in - $inc x |\
    snapgrid - - xrange=-$size:$size yrange=-$size:$size nx=$n ny=$n |\
    ccdsmooth - - $s |\
    ccdplot -  $cnt gray=t power=0.4 tab=- headline=$in abs=f |\
    tablsqfit - fit=ellipse > tmp.fit
if [ $debug -ge -1 ]; then
    echo "mag: $mag $cnt"
    echo "status rot/grid/smooth/plot/fit: ${PIPESTATUS[@]}"
    cat tmp.fit
fi
echo -n "$in $inc $q0 $contour "
txtpar tmp.fit p0=eccentr,1,5 \
       expr="%1,acos(%1)*180/pi,acos(sqrt((%1**2-$q0**2)/(1-$q0**2)))*180/pi,sqrt((%1**2-cosd($inc)**2)/(1-cosd($inc)**2))" 

# why should the vertical stellar distribution be isothermal, when dark matter is in play as well
