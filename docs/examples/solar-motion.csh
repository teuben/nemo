#! /bin/csh -f
#
#  emulate the solar motion in the rotating frame of reference in any of the GalPot models
#  See also 1998MNRAS.294..429D for a description of the models.
#
#    r0    solar radius
#    n     extension of any of the GalPot models (e.g. n=1 or n=2e
#
#  History:
#    21-mar-2005   Created                       Peter Teuben

set r0=8
set n=1
set gif=0
set clean=1

foreach arg ($*)
  set $arg
end


set potfile=$NEMODAT/GalPot/pot.$n
set rc=(`rotcurves GalPot 0 $potfile radii=$r0 tab=t debug=-1 yapp=/null`)
set omega=`nemoinp $rc[2]/$rc[1]`
set potential=(potname=GalPot potpars=-$omega potfile=$potfile)

if ($gif) then
  set yapp1=plot1.gif/gif
  set yapp2=plot2.gif/gif
else
  set yapp1=1/xs
  set yapp2=2/xs
endif  
  


echo r0=$r0 n=$n omega=$omega

set tmp=tmporb.$$

mkorbit - x=-$r0 y=0 z=0 vx=10.0 vy=5.2 vz=7.2 $potential |\
 orbint - $tmp.1  1000 0.002 ndiag=100 nsave=1
orbplot $tmp.1 xrange=-16:16 yrange=-16:16 yapp=$yapp1
otos $tmp.1 $tmp.2
snapplot $tmp.2 xrange=-32:32 yrange=-32:32 xvar=vr yvar=vt xlabel=-U ylabel=V yapp=$yapp2

if ($clean) rm -f $tmp.*
