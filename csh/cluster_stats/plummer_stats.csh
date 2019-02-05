#! /bin/csh -f
#
#
#  some fun with Plummer Spheres -   22-aug-2018   Peter Teuben
#  NEMO's mkplummer
#
#
#  plummer_stats.csh tmp=exp1 nexp=100 mlo=1 mhi=1
#  plummer_stats.csh tmp=exp2 nexp=100 mlo=1 mhi=1 tstop=10
#  plummer_stats.csh tmp=exp3 nexp=100 
#  plummer_stats.csh tmp=exp4 nexp=100 tstop=10
#
# @todo  rethink scaling procedure:  should it be done first on snapshot level, or last on plotting level
#        current is hybrid. awk.  
#  

set nbody  = 1000000       # nbody for large cluster baseline plot
set nsmall = 100           # nbody for small cluster
set nexp   = 2             # number of experiments with small cluster
set mlo    = 0.3           # lower mass cutoff
set mhi    = 5             # lower mass cutoff
set pow    = -2.0          # power low of IMF (as per massexpr=pow(m,p))
set rcut   = 5.0           # cutoff in virial radii of small cluster
set tstop  = 0             # use 10 or so if evolution is needed

set rscale = 0.5           # scale l to pc    [Cluster M=100M_solar r_v=0.5 pc]
set vscale = 1.07          # scale v to km/s
set rscale = 1             # 
set vscale = 1             # 
set rbin   = 0:2.5:0.25    # bins, in the $rscale (pc) units

set quick  = 0             # skip r_v and r_c

set tmp    = tmplummer     # baseline for files

# poor man's command line parser
foreach arg ($*)
   set $arg
end

# administrativia
rm -rf $tmp.*
setenv DEBUG -1



#  for a large number of nbody, make a "baseline" plummer and contour plot of (r,vt)
echo "=== Make a baseline large cluster/plot"
set scaling = (xscale=$rscale yscale=$vscale xlab='R2(pc)' ylab='V2(km/s)')

if ($nbody > 0) then
    mkplummer $tmp.0.dat $nbody 
    snapgrid $tmp.0.dat $tmp.rvt xvar=r2 yvar=v2 yrange=0:2 xrange=0:8
    ccdplot $tmp.rvt 0.001,0.003,0.01,0.03,0.1,0.3,1 $scaling yapp=$tmp.rvt.ps/vps
    snapgrid $tmp.0.dat $tmp.rvz xvar=r2 yvar=vz yrange=-2:2 xrange=0:8
    ccdplot $tmp.rvz 0.001,0.003,0.01,0.03,0.1,0.3,1 $scaling yapp=$tmp.rvz.ps/vps

    snapsort $tmp.0.dat - r  | snapshell - $rbin v   r  > $tmp.0.vtab
    snapsort $tmp.0.dat - r  | snapshell - $rbin vt  r  > $tmp.0.vttab
    snapsort $tmp.0.dat - r2 | snapshell - $rbin vr2 r2 > $tmp.0.vr2tab
    snapsort $tmp.0.dat - r2 | snapshell - $rbin vt2 r2 > $tmp.0.vt2tab
endif    


#
echo "=== Loop over $nexp experiments of $nsmall bodies each"
set scaling = (xscale=$rscale yscale=$vscale)
foreach i (`seq 1 $nexp`)
   echo $i
   if ($mlo == $mhi) then
     mkplummer $tmp.$i.dat $nsmall rfrac=$rcut seed=$i >& /dev/null   
   else
     mkplummer $tmp.$i.dat $nsmall rfrac=$rcut seed=$i massname='n(m)' masspars=p,$pow massrange=$mlo,$mhi >& /dev/null
   endif
   if ($tstop > 0) then
     hackcode1 $tmp.$i.dat $tmp.$i.edat tstop=$tstop > $tmp.$i.edat.log
     rm $tmp.$i.dat 
     snaptrim $tmp.$i.edat $tmp.$i.dat  times=$tstop
   endif

   if ($quick == 0) then
     hackforce_qp $tmp.$i.dat - | snapstat - all=t | grep r_v >> $tmp.rvtab
   endif
     
   snapsort $tmp.$i.dat  $tmp.$i.sdat r2
   snapprint $tmp.$i.sdat r2,m,vx,vy,vz,v2 >  $tmp.$i.tab
   tabplot $tmp.$i.tab 1 2 xbin=$rbin tab=t $scaling yapp=/null > $tmp.$i.2htab     
   tabplot $tmp.$i.tab 1 3 xbin=$rbin tab=t $scaling yapp=/null > $tmp.$i.3htab  
   tabplot $tmp.$i.tab 1 4 xbin=$rbin tab=t $scaling yapp=/null > $tmp.$i.4htab
   tabplot $tmp.$i.tab 1 5 xbin=$rbin tab=t $scaling yapp=/null > $tmp.$i.5htab
   tabplot $tmp.$i.tab 1 6 xbin=$rbin tab=t $scaling yapp=/null > $tmp.$i.6htab
   # tabplot $tmp.$i.htab 1 2 dycol=4 point=2,0.1 xmin=0 xmax=8 ymin=0 ymax=2 yapp=$tmp.$i.htab.ps/vps
   # tabhist $tmp.$i.tab 1 0 4 10 24 yapp=99/xs gauss=f residual=f
end

if ($quick == 0) then
  echo "=== Finding r_v statistics"
  awk '{print $5}' $tmp.rvtab | tabhist - debug=0 |& grep Mean
endif  

echo "=== Accumulating bin averages for all experiments"
#   the right way to average is via the bin avarages
#   this matrix transposing is a very slow step
set nbin=`getline $tmp.1.3htab`
foreach i (`seq 1 $nexp`)
    foreach j (`seq 1 $nbin`)
	getline $tmp.$i.2htab $j >> $tmp.$j.2ctab    
	getline $tmp.$i.3htab $j >> $tmp.$j.3ctab
	getline $tmp.$i.4htab $j >> $tmp.$j.4ctab
	getline $tmp.$i.5htab $j >> $tmp.$j.5ctab
	getline $tmp.$i.6htab $j >> $tmp.$j.6ctab	
    end
end

echo "=== Creating final table for ensemble average"
foreach j (`seq 1 $nbin`)
    set m = (`tabstat $tmp.$j.2ctab 1,2,5 | grep ^mean `)
    set d = (`tabstat $tmp.$j.2ctab 1,2,5 | grep ^disp `)
    echo $m[2] $m[3] $m[4]  $d[2] $d[3] $d[4] >> $tmp.2stab
    set m = (`tabstat $tmp.$j.3ctab 1,4,5 | grep ^mean `)
    set d = (`tabstat $tmp.$j.3ctab 1,4,5 | grep ^disp `)
    echo $m[2] $m[3] $m[4]  $d[2] $d[3] $d[4] >> $tmp.3stab
    set m = (`tabstat $tmp.$j.4ctab 1,4,5 | grep ^mean `)
    set d = (`tabstat $tmp.$j.4ctab 1,4,5 | grep ^disp `)
    echo $m[2] $m[3] $m[4]  $d[2] $d[3] $d[4] >> $tmp.4stab
    set m = (`tabstat $tmp.$j.5ctab 1,4,5 | grep ^mean `)
    set d = (`tabstat $tmp.$j.5ctab 1,4,5 | grep ^disp `)
    echo $m[2] $m[3] $m[4]  $d[2] $d[3] $d[4] >> $tmp.5stab
    set m = (`tabstat $tmp.$j.6ctab 1,2,5 | grep ^mean `)
    set d = (`tabstat $tmp.$j.6ctab 1,2,5 | grep ^disp `)
    echo $m[2] $m[3] $m[4]  $d[2] $d[3] $d[4] >> $tmp.6stab
    set m = (`tabstat $tmp.$j.6ctab 1,4,5 | grep ^mean `)
    set d = (`tabstat $tmp.$j.6ctab 1,4,5 | grep ^disp `)
    echo $m[2] $m[3] $m[4]  $d[2] $d[3] $d[4] >> $tmp.7stab
end

echo "=== Plot"
set scaling = (xscale=1.0     yscale=1.0     xlab='R2(pc)' ylab='V2(km/s)')  # if scaling done earlier
tabplot $tmp.2stab 1 3 dycol=6 point=2,0.1 xmin=0 xmax=4 ymin=0         $scaling yapp=$tmp.2stab.ps/vps
tabplot $tmp.3stab 1 2 dycol=5 point=2,0.1 xmin=0 xmax=4 ymin=0  ymax=2 $scaling yapp=$tmp.3stab.ps/vps
tabplot $tmp.4stab 1 2 dycol=5 point=2,0.1 xmin=0 xmax=4 ymin=0  ymax=2 $scaling yapp=$tmp.4stab.ps/vps
tabplot $tmp.5stab 1 2 dycol=5 point=2,0.1 xmin=0 xmax=4 ymin=0  ymax=2 $scaling yapp=$tmp.5stab.ps/vps
tabplot $tmp.6stab 1 2 dycol=5 point=2,0.1 xmin=0 xmax=4 ymin=0  ymax=2 $scaling yapp=$tmp.6stab.ps/vps
tabplot $tmp.7stab 1 2 dycol=5 point=2,0.1 xmin=0 xmax=4 ymin=0  ymax=2 $scaling yapp=$tmp.7stab.ps/vps

echo "=== Done."




###
#  A word on units. Simulation data are often in virial units
#                                    (Heggie & Mathieu (1986)
#  where G=1  M=1 E=-1/4
#  U_m = M                   
#  U_l = -G M^2 / (4E)       
#  U_t = G M^2.5 / (-4E)^1.5 
#
# Common starcluster units:
#  U_m = M_solar
#  U_l = pc
#  U_v = km/s
#  ->   G = 1 / 2.32385e2
#  ->   U_t ~ 1e6.yr
# Also:
# 1 pc = 3.086e18 cm
# 1 yr = 3.156e7 sec
