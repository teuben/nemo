#! /bin/csh -f
#
#   Construct a Plummer Sphere of Little Plummer Spheres
#
#   14-mar-2019   Created - Peter Teuben
#
#   Essential input parameters:                            PoP-math.tex
#   ---------------------------------------------------------------
#   n1   number of little plummers in the big plummer        N
#   n0   number of stars in a little plummer                 n
#   r0   radius of a little plummer                          \epsilon
#        (make sure n1*r0 not too much bigger than 1)


if ($?NEMO == 0) exit 1

# Big Plummer 

set n1 = 50
set m1 = 1        # don't change
set r1 = 1        # don't change
set s1 = 0        # could fix?

# Little Plummer (make sure that n1*r0 < 1 to ensure isolated plummers)

set n0 = 50
set m0 = `nemoinp $m1/$n1`
set r0 = 0.1
set s0 = -2       # to ensure they're all different?  (-1 and 0 too risky if CPU fast)

# 
set out = run0
set top = 0

#
set debug=-1

# parse command line, it has to be a series of "key=val"
foreach arg ($*)
   set $arg
end   


# don't change anything below
echo Writing to out=$out
echo n1=$n1 n0=$n0 r0=$r0
sleep 2
rm -rf $out.*

mkplummer $out.1 $n1 seed=$s1
snapprint $out.1 format=%.10f debug=$debug > $out.1pv

set vscale = `nemoinp "1/sqrt($n1*$r0)"`

foreach i (`seq $n1`)
   set pv = (`getline $out.1pv $i`)
   if ($top == 0) then
      set rshift = $pv[1],$pv[2],$pv[3]
      set vshift = $pv[4],$pv[5],$pv[6]
   else
      set rshift = 0,0,0
      set vshift = 0,0,0
   endif
   if ($debug > 0) then
     echo $i $rshift $vshift
   endif
   set out0 = `printf $out.0.%05d $i`
   mkplummer - $n0 seed=$s0 | snapscale - - rscale=$r0 mscale=$m0 vscale=$vscale | snapshift - $out0 $rshift $vshift
   echo $out0 >> $out.list
end

snapadd @$out.list $out.2 debug=$debug

snapstat $out.1 all=t exact=t eps=0.0 > $out.1.stat
snapstat $out.2 all=t exact=t eps=0.0 > $out.2.stat
grep r_v $out.1.stat
grep r_v $out.2.stat
set rvf=`nemoinp "1/(1+1/($n1*$r0))"`
echo Predicted final virial radius r_v : $rvf
grep 2T/W $out.1.stat | awk '{print "2T/W r_v : ",$9}'
grep 2T/W $out.2.stat | awk '{print "2T/W r_v : ",$9}'

echo vscale=$vscale

if (0) then
  echo "Now try something like:"
  echo "gyrfalcON run1.2 run1.2.out tstop=100 step=0.2 eps=0.005 kmax=7 give=mxvap"
  echo "snapcenter run1.2.out run1.2.center weight='-phi*phi*phi'"
  echo "snapmradii run1.2.center 0.01,0.1:0.9:0.1,0.99 log=t > run1.2.mr"
  echo "tabplot run1.2.mr 1 2:12 line=1,1 color=2,3::9,2"
endif
