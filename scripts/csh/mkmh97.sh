#!/bin/bash
#
#    Example of setting up a Makino & Hut 1997 simulation, but with unequal masses
#    See also: https://ui.adsabs.harvard.edu/abs/1997ApJ...481...83M/abstract
#
# bench:
#          ./mkmh97.sh code=0 run=run100 nbody=10000 tstop=20 seed=123 
#          ./mkmh97.sh code=1 run=run101 nbody=10000 tstop=20 seed=123
#                     i7-1185G7  Xeon E5-2687W      i9-12900K
#          code=0     5:12 min       10:48            3:31
#          code=1     2:10 min        3:20            1:25
#          code=2     ----            0:16            ----
#
# Note G1 (the one with mass=1) starts at x>0 and launched with vx<0, moving to the left
#      G2 (the one with mass=m) starts at x<0 and launched with vx>0, moving to the right
# For fixed=1 G1 will be fixed at (0,0), where G2 is launched at x<0 

#
# version: 12-may-2022   initial version with just (near) head-on collision
#          15-may-2022   added option to make (near) circular orbit
#          15-jun-2022   added m= mass ratio parameter, changed hack= -> code= and implemented bonsai
#          17-jun-2022   added seed= and defined a benchmark
#          20-jun-2022   add potentials & acc
#          27-jun-2022   $run is now a directory
#           7-jul-2022   $em=1 as option to have equal mass particles
#          11-jul-2022   able to rerun if the run directory exists
#          18-jul-2022   add r-vr evolution plot; new defaults for some parameters
#           8-aug-2022   store a table with time,x1,vx1,x2,vx2
#          11-aug-2022   using more generic nemopars.rc
#          20-aug-2022   add --help option in a neat self-documenting way
#          13-sep-2022   set m16=0 when no stars of G2 near G1
#          11-nov-2022   align integration parameters w/ MH97
#           5-dec-2022   label the table columns where needed
#          14-dec-2022   try to unbind G2 in it's own C.O.M., store etot
#          13-feb-2023   add trim option, using nemo_functions.sh now, add $_date
#          25-feb-2023   add detot
#          27-feb-2023   compute m1, etot2, plus add xve2.tab and xve0.tab, and
#          29-mar-2023   add option to use a fixed potential for galaxy1
#          29-apr-2023   fix orbital energies for fixed=1, also save etot0,v
#          22-jun-2023   added code=3 for princeton hackathon comparisons

_script=mkmh97
_version=29-nov-2023
_pars=nemopars.rc
_date=$(date +%Y-%m-%dT%H:%M:%S)     # now

#            text between #--HELP and #--HELP is displayed when --help is used
#--HELP
#            parameters for the integration   [defaults are following Makino & Hut 1997]
run=run0        # directory and basename of the files belonging to this simulation  [ empty means self] @todo
nbody=2048      # number of bodies in one model
m=1             # mass of second galaxy (mass of first will always be 1)
em=0            # equal mass particles? (em=0 means nbody same for both galaxies, thus individual masses not equal)
fixed=0         # fixed potential for G1 ?
potname=plummer       # fixed potential name for G1
potpars=0,1,3*pi/16   # fixed potential parameters for G1
step=1          # step in time when to dump snapshots
v0=1.0          # initial impact/circular speed
rp=0.0          # impact offset radius
r0=16.0         # initial offset position for v0 > 0   [note v0 > 2/sqrt(r0)]
eps=0.03125     # gravitational softening
kmax=7          # integration timestep is 1/2**kmax
eta=0.01        # extra timestep accuracy (only used by gravidy)
code=1          # 0=hackcode1  0q=hackcode1_qp  1=gyrfalcON   2=bonsai2 (GPU)  3=rungravidy  4=directcodeb
seed=0          # random seed (use seed=123 for the benchmark)
trim=0          # trim the data after analysis so only the last snapshot is kept
save=1          # save the "final" plots in a subdirectory "movies" labeled with tstop
#             parameters for the analysis (which can be re-run)
tstop=50        # stop time of the integration (or analysis time when doing a re-run)
box=32          # spatial box size for plotting and CCD frames
r16=16          # radius within which to measure G2 on top of G1
vbox=2          # velocity box size
npixel=128      # number of pixels in xy CCD frame
power=0.5       # gamma factor for CCD plots
bsigma=0.0001                      # asinh/log breakover point
tplot=0,5,10,15,20,25,30,40,50     # times to plot in evolution
hackforce=hackforce                # which hackforce to use in the force re-calculation
yapp=png                           # pick png, or vps (for yapp_pgplot) _ps for native ps
progress=0                         # spawn a progress bar for the nbody run
debug=1                            # 1=set -x,-e,-u   0=nothing
#
#--HELP

save_vars="run nbody m em fixed step v0 rp r0 eps kmax eta code seed trim tstop box r16 vbox npixel power bsigma tplot yapp"

if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    set +x
    awk 'BEGIN{s=0} {if ($1=="#--HELP") s=1-s;  else if(s) print $0; }' $0
    exit 0
fi

#             simple keyword=value command line parser for bash
for arg in "$@"; do
  export "$arg"
done


# portable yapp
yapp() {
    if test $yapp = "xs"; then
        echo $1/$yapp
    elif test $yapp = "_ps"; then
        echo $1.ps
    else
        echo $1.$yapp/$yapp
    fi
}

if [ $debug == 1 ]; then
    set -x
    set -e
    set -u
fi
    

#  @todo  if an explicit  (not advertised) restart was requested

#             delete the old run, work within the run directory
if [ -d $run ]; then
    restart=0
else
    restart=1
fi
if [ $run != . ]; then
    mkdir -p $run
    cd $run
else
    if [ ! -e $_pars ]; then
	echo "No $_pars in $(pwd)"
	exit 0
    fi
    
fi
# backwards compatible!
if [ -e mkmh97.rc ]; then
    echo This script is not compatible with this run
    exit 0
fi
# keep track of history  @todo need to bootstrap current $run
[[ -e $_pars ]] && source $_pars
if [ ! -e $_pars ]; then
    source nemo_functions.sh
    show_vars $save_vars >> $_pars
fi
echo "# $0 version=$_version"  >> $_pars
echo "# date=$_date"           >> $_pars
echo "$*" run=$run             >> $_pars
source  $_pars

# derived/fixed parameters
r=-$box:$box

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
	if [ $fixed = 0 ]; then
	    snapstack $run.1 $run.2 $run.3 deltar=$r0,$rp,0 deltav=-$v0,0,0  zerocm=t
	else
	    # 1 system only
	    snapshift $run.2 $run.3 rshift=-$r0,$rp,0 vshift=$v0,0,0
	    # 2 systems, but mass of G1 is 0
	    # snapscale $run.1 - mscale=0  | snapstack - $run.2 $run.3 deltar=$r0,$rp,0 deltav=-$v0,0,0  zerocm=t
	fi
    else
	# (near) circular orbit
	snapstack $run.1 $run.2 $run.3 deltar=$r0,0,0   deltav=0,$v0,0   zerocm=t
    fi

    # compute the orbit of a massless G2 in a fixed potential
    dt=0.01
    echo "Computing fixed-path orbit with dt=$dt"
    mkorbit - -$r0 $rp 0 $v0 0 0 potname=$potname potpars=$potpars |\
	orbint - - dt=$dt nsteps=$tstop/$dt nsave=$step/$dt |\
	orblist -  > fixed-path.tab 2>&1
    

    # integrator:  0:  hackcode1 is O(NlogN) code 
    #              1:  gyrfalcON is O(N)
    #              2:  bonsai2 is O(N) but scales faster for "small" N
    #              3:  rungravidy O(N^2) code
    echo "Use:   tail -f $run/$run.4.log     to monitor progress of the integrator"
    if [ $fixed = 1 ]; then
	echo "Fixed potential case; only code=0,1 are supported"
	if [ $code = 0 ]; then
	    hackcode3 $run.3 $run.4 eps=$eps freq=2**$kmax freqout=1/$step fcells=2 tstop=$tstop options=mass,phase,phi,acc \
		      potname=$potname potpars=$potpars   > $run.4.log
	    snapdiagplot $run.4 tab=$run.4.etot
	elif [ $code = 1 ]; then
	    gyrfalcON $run.3 $run.4 eps=$eps kmax=$kmax step=$step tstop=$tstop give=mxvap \
		      accname=$potname accpars=$potpars    > $run.4.log
	    tabcols $run.4.log 1,2 > $run.4.etot
	fi
	# exit 0
    else
	if [ $code = 0 ]; then
	    hackcode1 $run.3 $run.4 eps=$eps freq=2**$kmax freqout=1/$step fcells=2 tstop=$tstop options=mass,phase,phi,acc > $run.4.log    
	    snapdiagplot $run.4 tab=$run.4.etot
	elif [ $code = 0q ]; then
	    hackcode1_qp $run.3 $run.4 eps=$eps freq=2**$kmax freqout=1/$step fcells=2 tstop=$tstop options=mass,phase,phi,acc > $run.4.log    
	    snapdiagplot $run.4 tab=$run.4.etot
	elif [ $code = 1 ]; then
	    progress_bar.sh $progress $run.4.log 1 $tstop &
	    gyrfalcON $run.3 $run.4 eps=$eps kmax=$kmax step=$step tstop=$tstop give=mxvap > $run.4.log
	    tabcols $run.4.log 1,2 > $run.4.etot
	elif [ $code = 2 ]; then
	    # See $NEMO/usr/bonsai
	    snaptipsy  $run.3 $run.3t
	    ulimit -s unlimited
	    bonsai2 -i $run.3t --snapname $run. --logfile $run.gpu --dev 0 --snapiter $step -T $tstop -t $(nemoinp 1/2**$kmax)  1>$run.4.log 2>$run.4.err
	    grep Etot= $run.4.log  |tabcols - 4,6 > $run.4.etot
	    # since bonsai writes a new snapshot (in tipsy format) for each time, need to convert them
	    # to single NEMO snapshot for the same analysis as the other codes.
	    set +x
	    for f in $run._*; do
		echo Processing $f
		(tipsysnap $f - | csf - - item=SnapShot 1>>$run.4 ) 2>>$run.4.t2s
	    done
	    rm -f $run._*
	    if [ $debug == 1 ]; then
		set -x
	    fi	    
	    echo Wrote final combined snapshot in $run.4
	elif [ $code = 3 ]; then
	    progress_bar.sh $progress $run.4.d/$run.4.d.out.log 2 $tstop &	    
	    /usr/bin/time rungravidy $run.3 $run.4.d tcrit=$tstop deltat=$step eps=$eps eta=$eta > $run.4.log
	    ln -s $run.4.d/OUT3.snap $run.4
	    tabcols $run.4.d/$run.4.d.out.log 2,5 > $run.4.etot
	    # @todo ### Fatal error [snapcenter]: total weight is zero
	elif [ $code = 4 ]; then
	    directcode $run.3 $run.4 eps=$eps freq=2**$kmax freqout=1/$step tstop=$tstop options=mass,phase,phi,acc > $run.4.log    
	    snapdiagplot $run.4 tab=$run.4.etot
	else
	    set +x
	    echo "Unknown code=$code, valid are:"
	    echo "0 = hackcode1"
	    echo "1 = gyrfalcON"
	    echo "2 = bonsai2"
	    echo "3 = rungravidy"
	    echo "4 = directcode"
	fi
    fi
else
    echo "============================================================================="
    echo "Skipping integration because run=$run already existed, or: rm -r $run"
fi

#  special case for tstop=0:   just make initial conditins

if [ $tstop = 0 ]; then
    exit
fi


#  compute the path of G1 and G2.   G2 needs a special treatment if mass G2 << G1

if [ ! -e $run.xv.tab ]; then
    #
    echo "# t"                              >$run.4.t.tab
    snapcopy $run.4 - i=0 | snapprint - t  >>$run.4.t.tab

    #
    echo "# x1 y1 z1 vx1 vy1 vz1"                                                                                   >$run.4.g1.tab
    if [ $code == 2 ] || [ $code == 3 ]; then
	snapcopy $run.4 - "select=i<$nbody?1:0" | $hackforce - - debug=-1 | snapcenter - . "-phi*phi*phi" report=t >>$run.4.g1.tab
    else
	snapcenter $run.4 . "weight=i<$nbody?-phi*phi*phi:0" report=t                                              >>$run.4.g1.tab
    fi
    
    #
    if [ $fixed = 0 ]; then
	echo "# x2 y2 z2 vx2 vy2 vz2"                                                                                  >$run.4.g2.tab    
	snapcopy $run.4 - "select=i>=$nbody?1:0" | $hackforce - - debug=-1 | snapcenter - . "-phi*phi*phi" report=t   >>$run.4.g2.tab
    else
	# fake it
	cp $run.4.g1.tab $run.4.g2.tab
    fi
    #
    echo "# t x1 v1 x2 v2"                                                                        >$run.xv.tab
    paste  $run.4.t.tab $run.4.g1.tab $run.4.g2.tab | grep -v ^# | awk '{print $1,$2,$5,$8,$11}' >>$run.xv.tab

    
    #  plot the path in Pos and Vel separately
    tabplot $run.xv.tab 1 2,4 line=1,1 color=2,3 ycoord=0 yapp=$(yapp path-pos)
    tabplot $run.xv.tab 1 3,5 line=1,1 color=2,3 ycoord=0 yapp=$(yapp path-vel)

    # compare to the fixed-path
    tabcols $run.xv.tab 1,4,5    > fixed1.tab
    tabcols fixed-path.tab 2,3,6 > fixed2.tab
    paste fixed1.tab fixed2.tab | tabplot - 1 2,5 line=1,1 color=3,2 ycoord=0 yapp=$(yapp path-pos-fixed)
    paste fixed1.tab fixed2.tab | tabplot - 1 3,6 line=1,1 color=3,2 ycoord=0 yapp=$(yapp path-vel-fixed)

else
    echo "Skipping path computation since it's already done, or: rm $run/$run.xv.tab"
fi


#  mean and sigma of trended energy
estats=$(tabtrend $run.4.etot 2 | tabstat - qac=t  | txtpar -  p0=QAC,1,3 p1=QAC,1,4)

tabplot $run.4.etot 1 2  headline="mean/sigma $estats" yapp=$(yapp etot.plot) > /dev/null 2>&1
tabhist $run.4.etot 2                                  yapp=$(yapp etot.hist) > etot.hist.log 2>&1
snapplot $run.3 xrange=$r yrange=$r                                yapp=$(yapp init.plot)
snapplot $run.4 xrange=$r yrange=$r times=$tstop                   yapp=$(yapp final.plot)
snapplot $run.4 xrange=$r yrange=$r times=$tstop visib="i<$nbody"  yapp=$(yapp final1.plot)
snapplot $run.4 xrange=$r yrange=$r times=$tstop visib="i>=$nbody" yapp=$(yapp final2.plot)
    
snapgrid $run.3 - xrange=-$box:$box yrange=-$box:$box              nx=$npixel ny=$npixel |\
    ccdmath - - "log(1+%1/$bsigma)" |\
    ccdplot - power=$power yapp=$(yapp init.ccd) headline="Initial Conditions"
snapgrid $run.4 - xrange=-$box:$box yrange=-$box:$box times=$tstop nx=$npixel ny=$npixel evar=m |\
    tee final.ccd |\
    ccdmath - - "log(1+%1/$bsigma)" |\
    ccdplot - power=$power yapp=$(yapp final.ccd) headline="Conditions at tstop=$tstop"
snapgrid $run.4 - xrange=-$box:$box yrange=-$box:$box times=$tstop nx=$npixel ny=$npixel evar="i<$nbody?m:0" |\
    tee final1.ccd |\
    ccdmath - - "log(1+%1/$bsigma)" |\
    ccdplot - power=$power yapp=$(yapp final1.ccd) headline="Galaxy-1 at tstop=$tstop"
snapgrid $run.4 - xrange=-$box:$box yrange=-$box:$box times=$tstop nx=$npixel ny=$npixel evar="i>=$nbody?m:0" |\
    tee final2.ccd |\
    ccdmath - - "log(1+%1/$bsigma)" |\
    ccdplot - power=$power yapp=$(yapp final2.ccd) headline="Galaxy-2 at tstop=$tstop"

#  convert the tee'd ccd files to fits
rm -f final.fits final1.fits final2.fits
ccdfits final.ccd  final.fits  radecvel=t
ccdfits final1.ccd final1.fits radecvel=t
ccdfits final2.ccd final2.fits radecvel=t

#  final plotting
snapplot  $run.4 xrange=-$box:$box yrange=-$box:$box                   times=$tplot nxy=3,3 yapp=$(yapp evolution-xy.plot)
snapplot  $run.4 xrange=0:$box yrange=-$vbox:$vbox xvar=r yvar=vr      times=$tplot nxy=3,3 yapp=$(yapp evolution-vr.plot)
snapplot3 $run.4 xrange=-$box:$box yrange=-$box:$box zrange=-$box:$box times=$tstop         yapp=$(yapp final.3d.plot)

#  final G2 snapshot for analysis:
#  center G2 on G1, and plot and show the cumulative mass fraction as function of radius
if [ $fixed == 0 ]; then
    snaptrim $run.4 - times=$tstop | snapcopy - - "select=i>=$nbody" > final2.snap
else
    snaptrim $run.4 - times=$tstop                                   > final2.snap
fi
x1=$(grep -w ^$tstop $run.xv.tab | txtpar - p0=1,2)
v1=$(grep -w ^$tstop $run.xv.tab | txtpar - p0=1,3)
snapshift final2.snap - $x1,0,0 $v1,0,0 mode=sub > final2c.snap
radprof final2c.snap tab=t > final2c.tab
echo '# radius mass'                           >final2cm.tab
tabmath final2c.tab - %1,%4/$m all format=%f  >>final2cm.tab
# If the first radius (star) is not within 16, there's no G2 stars near G1
r1=$(txtpar final2cm.tab 'iflt(%1,16,1,0)' p0=1,1)
if [ $r1 = 1 ]; then
    set +e
    tabspline final2cm.tab    x=1:15:2
    m16=$(tabspline final2cm.tab    x=$r16 | txtpar - p0=1,2)
else
    m16=0
fi
tabplot final2cm.tab  1 2 0 16 xlab=Radius ylab=Mass  headline="x1=$x1 v1=$v1 m16=$m16"  yapp=$(yapp massg2g1)
echo "m16=$m16" >> $_pars

# center on G1, unbind stars
snaptrim $run.4 - times=$tstop | snapcopy - - "select=i<$nbody" > final1.snap
x1=$(grep -w ^$tstop $run.xv.tab | txtpar - p0=1,2)
v1=$(grep -w ^$tstop $run.xv.tab | txtpar - p0=1,3)
snapshift final1.snap - $x1,0,0 $v1,0,0 mode=sub |\
   $hackforce - - > final1c.snap
unbind final1c.snap - > final1u.snap
m1=$(snapmstat final1u.snap | txtpar - p0=TotMas,1,8) 
snapplot final1u.snap xrange=-$box:$box yrange=-$box:$box yapp=$(yapp final1u.plot)
snapgrid final1u.snap - xrange=-$box:$box yrange=-$box:$box nx=$npixel ny=$npixel evar=m |\
    tee final1u.ccd |\
    ccdmath - - "log(1+%1/$bsigma)" |\
    ccdplot - power=$power yapp=$(yapp final1u.ccd) headline="Galaxy-1 bound at tstop=$tstop"
echo m1=$m1 >> $_pars
# binding plots
snapplot  final1c.snap color='phi+0.5*(vx*vx+vy*vy+vz*vz) > 0 ? 0.1 : 0.2' xrange=$r yrange=$r                                       yapp=$(yapp final1c.plot)
snapplot3 final1c.snap color='phi+0.5*(vx*vx+vy*vy+vz*vz) > 0 ? 0.1 : 0.2' xrange=$r yrange=$r zrange=$r                             yapp=$(yapp final1c.3d.plot)
snapplot  final1c.snap color='phi+0.5*(vx*vx+vy*vy+vz*vz) > 0 ? 0.1 : 0.2' xvar=r yvar=vr xrange=0:$box yrange=-$vbox:$vbox nxy=3,3  yapp=$(yapp final1c.evolution-binding-vr.plot)
snapplot  final1c.snap color='phi+0.5*(vx*vx+vy*vy+vz*vz) > 0 ? 0.1 : 0.2' xvar=x  yvar=y xrange=$r yrange=$r nxy=3,3                yapp=$(yapp final1c.evolution-binding-xy.plot)



# center on G2, unbind stars
x2=$(grep -w ^$tstop $run.xv.tab | txtpar - p0=1,4)
v2=$(grep -w ^$tstop $run.xv.tab | txtpar - p0=1,5)
echo "PJT: $x2 $v2"
snapshift final2.snap - $x2,0,0 $v2,0,0 mode=sub |\
    $hackforce - - > final2c.snap
unbind final2c.snap - > final2u.snap
m2=$(snapmstat final2u.snap | txtpar - p0=TotMas,1,8) 
snapplot final2u.snap xrange=-$box:$box yrange=-$box:$box yapp=$(yapp final2u.plot)
snapgrid final2u.snap - xrange=-$box:$box yrange=-$box:$box nx=$npixel ny=$npixel evar=m |\
    tee final2u.ccd |\
    ccdmath - - "log(1+%1/$bsigma)" |\
    ccdplot - power=$power yapp=$(yapp final2u.ccd) headline="Galaxy-2 bound at tstop=$tstop"
echo m2=$m2 >> $_pars
# binding plots
snapplot  final2c.snap color='phi+0.5*(vx*vx+vy*vy+vz*vz) > 0 ? 0.1 : 0.2' xrange=$r yrange=$r                                       yapp=$(yapp final2c.plot)
snapplot3 final2c.snap color='phi+0.5*(vx*vx+vy*vy+vz*vz) > 0 ? 0.1 : 0.2' xrange=$r yrange=$r zrange=$r                             yapp=$(yapp final2c.3d.plot)
snapplot  final2c.snap color='phi+0.5*(vx*vx+vy*vy+vz*vz) > 0 ? 0.1 : 0.2' xvar=r yvar=vr xrange=0:$box yrange=-$vbox:$vbox nxy=3,3  yapp=$(yapp final2c.evolution-binding-vr.plot)
snapplot  final2c.snap color='phi+0.5*(vx*vx+vy*vy+vz*vz) > 0 ? 0.1 : 0.2' xvar=x  yvar=y xrange=$r yrange=$r nxy=3,3                yapp=$(yapp final2c.evolution-binding-xy.plot)


# compare m16 and m2 as function of time
if [ -e $run.xvm.tab ]; then
    sed s/nan/000/ $run.xvm.tab |\
	tabmath - - %1,%6*$m,%7 all |\
	tabplot - 1 2,3 line=1,1 color=2,3 ycoord=0,$m yapp=$(yapp path-g2)
fi

# binding energy argument, plot total binding energy as function of time
echo "# time x1 v1 x2 v2 kin1 kin2 pot12 etot"                              > $run.xve.tab
if [ $fixed = 0 ]; then
    tabmath $run.xv.tab - "0.5*%3**2,0.5*${m}*%5**2,-${m}/abs(%2-%4),%6+%7+%8" >> $run.xve.tab
else
    tabmath $run.xv.tab - "0,0.5*${m}*%5**2,-${m}/abs(%2),%6+%7+%8" >> $run.xve.tab    
fi
tabplot $run.xve.tab 1 9 line=1,1 ycoord=0 yapp=$(yapp path-energy) xlab=Time ylab=Energy



#plot by binding energy
r=-$box:$box
snapplot $run.4 color='phi+0.5*(vx*vx+vy*vy+vz*vz) > 0 ? 0.1 : 0.2' xrange=$r yrange=$r times=$tstop yapp=$(yapp final-binding.plot)
snapplot3 $run.4 color='phi+0.5*(vx*vx+vy*vy+vz*vz) > 0 ? 0.1 : 0.2' xrange=$r yrange=$r zrange=$r times=$tstop yapp=$(yapp final-binding.3d.plot)
snapplot $run.4 color='phi+0.5*(vx*vx+vy*vy+vz*vz) > 0 ? 0.1 : 0.2' xvar=r yvar=vr xrange=0:$box yrange=-$vbox:$vbox nxy=3,3 times=$tplot yapp=$(yapp evolution-binding-vr.plot)
snapplot $run.4 color='phi+0.5*(vx*vx+vy*vy+vz*vz) > 0 ? 0.1 : 0.2' xvar=x  yvar=y xrange=$r yrange=$r nxy=3,3 times=$tplot yapp=$(yapp evolution-binding-xy.plot)



# alternative orbital energy based on time=0
x10=$(nemoinp "$r0*$m/(1+$m)")
v10=$(nemoinp "-$v0*$m/(1+$m)")
x20=$(nemoinp "-$r0/(1+$m)")
v20=$(nemoinp "$v0/(1+$m)")
echo "0 $x10 $v10 $x20 $v20" | tabmath - - "0.5*%3**2,0.5*${m}*%5**2,-${m}/abs(%2-%4),%6+%7+%8"       > $run.xve0.tab
etot0=$(tabcols $run.xve0.tab 9)

if [ $fixed = 0 ]; then
    tail -1 $run.xv.tab  | tabmath - - "0.5*${m1}*%3**2,0.5*${m2}*%5**2,-${m1}*${m2}/abs(%2-%4),%6+%7+%8" > $run.xve2.tab
else
    tail -1 $run.xv.tab  | tabmath - - "0,0.5*${m1}*%5**2,-${m1}/abs(%2),%6+%7+%8" > $run.xve2.tab    
fi
etot2=$(tabcols $run.xve2.tab 9)

echo "Final binding energy behavior:"
etot=$(tail -10 $run.xve.tab | tabstat - 9 qac=t label=Etot | txtpar -  p0=QAC,1,3)
detot=$(tail -10 $run.xve.tab | tabstat - 9 qac=t label=Etot | txtpar -  p0=QAC,1,4)
echo "etot=$etot"   >> $_pars
echo "detot=$detot" >> $_pars
echo "etot2=$etot2" >> $_pars
echo "etot0=$etot0" >> $_pars


echo "m16=$m16"
echo "etot=$etot +\- $detot"

#   compute the v from infinity, give it a negative sign if "imaginary"
v2=$(nemoinp "$v0**2-2*(1+$m)/$r0")
v=$(nemoinp "sqrt(abs($v2))*sign($v2))")
echo "v=$v" >> $_pars

if [ $save == 1 ]; then
    label=$(printf %04.0f $(nemoinp 10*$tstop))
    echo "$run/movies label: $label"
    mkdir -p movies
    for p in final*png; do
	n=$(echo $p | sed s/.png/_${label}.png/)
	cp $p movies/$n
    done
fi

if [ $trim == 1 ]; then
    mv $run.4 $run.40
    snaptrim $run.40 $run.4 times=$tstop
    rm $run.40
fi

#--HELP

#    BENCHMARK:
#    rm -rf run0
#    run=run0 seed=123      should give:    m16=0.821802 for phi*phi*phi*phi
#    run=run0 seed=123      should give:    m16=0.861046 for -phi*phi*phi  (now the default)
#    nemopars m16 run0/nemopars.rc



#    NemoPlots with relevant scales for this script

# 1. plummer scaling to retain virial equilibrium
#PLOT nemoinp 0:1:0.01 | tabmath - - 'sqrt(%1),sqrt(%2)' | tabplot -  1 2,3 0 1 0 1 yapp=30/xs line=1,1 color=2,3 nxticks=9 nyticks=9 xlab=m ylab="Scale Factor" headline="Scale factor for position (R) and velocity (G)"

# 2. v0 vs. v (at infinity) for r0=16
#PLOT nemoinp 0:1.6:0.01 | tabmath - - "0.01,0.1,1,sqrt(%1**2+2*(1+%2)/16),sqrt(%1**2+2*(1+%3)/16),sqrt(%1**2+2*(1+%4)/16)"|tabplot - 1 5,6,7,1 0 1.6 0 1.6 color=2,3,4,1 line=1,1 xlab=V ylab=V0 headline="V0 for m=0.01(r), 0.1(g) and 1.0(b)"

#--HELP
