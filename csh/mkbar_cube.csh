#! /bin/csh -f
#
#  This scripts takes an HDF output snapshot file from cmhog
#  (bar hydro, polar coordinates), projects it to a requested
#  sky view, and using WIP summarizes the results
#
#  In the model the bar is assumed oriented along Y axis (PA=0), and 
#  flows CCW as seen from the positive Z axis.
#  For CW rotating galaxies the code adds 180 to 'inc' and/or 'pa'
#
#  Observations are assumed to have their reference pixel to
#  coincide with the center (ra0,dec0,vsys), although we expect
#  to relax this condition in a future version of this script.
#
#  The new convention for rot=1 means CCW, rot=-1 means CW, denoting
#  the sign of the galaxy angular momentum vector, where the positive
#  Z axis points to the observer (hence doppler recession is -vz).
#
#  23-sep-95	Created				Peter Teuben
#  15-nov-95    beam=0.2 frang=45
#  24-nov-95    defaults for more central region, fixed dependancies
#   9-jul-96    hacking for N5383 
#  24-jul-00    BIMA proposal N4303 et al.
#   5-sep-00    modified to write cubes instead of moment maps
#  12-mar-01    radecvel=t to make karma swallow these fits files
#  13-mar-01    use phi,inc,pa (no more +/- 180) and documented geometry
#               (notice that earlier versions had sign of radial vel wrong)
#  23-mar-01    added a refmap and fixed refscale; this assumes that the refmap
#		(often a cube) has the reference pixel defined to the be 
#		center of the galaxy (bimasong data often don't do the VELO axis correct)
#  17-apr-01    added velocity referencing using 'vsys' to be at v=0 in the model cube
#               (bugs when model and data have different delta-V)
#  10-apr-02    more generic version without a reference cube
#   6-jul-02    added back reference cube code with lots of changes, added inden= keyword

if ($#argv == 0) then
  echo Usage: $0 in=HDF_DATASET out=BASENAME ...
  echo Gridding and projecting 2D CMHOG hydro models to specified bar viewing angles
  exit 0
endif

# 			Required Keywords
unset in
unset out
# 			Geometry (defaults are for some reasonable galaxy)
set pa=30
set inc=60
set phi=30
set rot=1
#                       Spatial gridding (n cells from -rmax : rmax)
set rmax=6
set n=200
set beam=0.25
set color=1
set clean=1
set cube=0
set denlog=0

#			Velocity gridding for cube  (nvel cells from -vmax : vmax)
set nvel=50
set vmax=250

#                       WCS definition (and mapping) if derived from an observation (cube)
set wcs=0
set refmap=""
set pscale=1
set vscale=1
set vsys=0


#
set par=""
set inden=""
#			Parse commandline to (re)set keywords
foreach a ($*)
  set $a
  if (-e "$par") then
    source $par
  else if (X != X$par) then
    echo Warning, par=$par does not exist
  endif
  set par=""
end

#               Fixed constants (p: arcsec -> degrees    v: km/s -> m/s)
#               probably should not be changed
set puscale=2.77777777777e-4
set vuscale=1e3

#
#  fix inc/pa for ccw(rot=1) or cw(rot=-1) cases for NEMO's euler angles
#
if ($rot == -1) then
   set inc=$inc+180
else if ($rot == 1) then
   set pa=$pa+180
else
   echo "Bad rotation, must be 1 (ccw) or -1 (cw)"
   exit 1
endif
#		Report
echo     Files: in=$in out=$out 
echo -n "Grid: rmax=$rmax n=$n beam=$beam "
if ($wcs) then
   echo \(`nemoinp "$beam*$pscale"` arcsec\)
else
   echo ""
endif
echo Projection: phi=$phi inc=$inc pa=$pa
#               Derived quantities
set cell=`nemoinp "2*$rmax/$n"`
set range="-${rmax}:${rmax}"
echo -n "      Derived: cell=$cell"
if ($wcs) then
   echo \(`nemoinp "2*$rmax/$n*$pscale"` arcsec\)
else
   echo ""
endif

if ($wcs) then
  if (-e "$refmap") then
    #   referencing for the 3rd axis
    set nz=(`fitshead $refmap | grep ^NAXIS3 | awk '{print $3}'`)
    set pz=(`fitshead $refmap | grep ^CRPIX3 | awk '{print $3}'`)
    set vz=(`fitshead $refmap | grep ^CRVAL3 | awk '{print $3}'`)
    set dz=(`fitshead $refmap | grep ^CDELT3 | awk '{print $3}'`)
    
    # step in model cube, if vscale=1
    set dz1=`nemoinp "2*$vmax/$nz"`
    # 
    set vref=`nemoinp "($vz-$vsys*$vuscale)/($dz1)+$nvel/2+0.5"`
    #set vscale=`nemoinp "$vscale*(2*$vmax/$nvel)/($dz/1000)"`
    set refcen=`nemoinp $n/2+0.5`
    set vref=`nemoinp $nvel/2+0.5`
    set crpix=$refcen,$refcen,$vref

    #  cdelt in new units
    set dv=`nemoinp "2*$vmax*$vuscale/$nvel*$vscale"`
    set dp=`nemoinp "2*$rmax*$puscale/$n*$pscale"`
    set cdelt=-$dp,$dp,$dv

    # RA/DEC: just report how much the reference pixel in the reference cube is offset from out (ra0,dec0)
    set ra00=`echo $ra0 | awk -F: '{printf("%.10f\n",((($3/60+$2)/60)+$1)*15)}'`
    set dec00=`echo $dec0 | awk -F: '{printf("%.10f\n",(($3/60+$2)/60)+$1)}'`
    echo RA0,DEC0=$ra0,$dec0
    echo RA00,DEC00=$ra00,$dec00
    set vx=(`fitshead $refmap | grep ^CRVAL1 | awk '{print $3}'`)
    set vy=(`fitshead $refmap | grep ^CRVAL2 | awk '{print $3}'`)
    echo Offset in RA and DEC: `nemoinp "(($vx)-($ra00))*3600"`  `nemoinp "(($vy)-($dec00))*3600"` arcsec
    set crval=$ra00,$dec00,`nemoinp "$vsys*$vuscale"`
    
    echo $nz $pz $vz $dz 
    echo Vsys at OBS pixel: `nemoinp "($vsys*$vuscale-$vz)/$dz+$pz"` 
    echo CRPIX:   $crpix
    echo CRVAL:   $crval
    echo CDELT:   $cdelt
    set wcspars=(refmap=$refmap crpix=$crpix crval=$crval cdelt=$cdelt radecvel=t)
    #  BUG:   the mod-cube is assumed to be centered at their own WCS of (0,0,0)
    #         the obs-cube has a positional reference pixel which is assumed
    #         to coincide with (0,0), however
  else
    echo BUG: need to rewrite this section for when no refmap given..... since you did not
    exit 1
  endif    
else
  set wcspars=()
endif


#> nemo.need tabtos tabmath snaptrans snaprotate snapadd snapgrid ccdsmooth ccdmath ccdfits fitshead

set tmp=tmp$$
if (! -e $out.den.fits) then

    # convert the half-plane HDF file to a full plane snapshot file

    tsd in=$in out=$tmp.tab coord=t
    if (-e "$inden") then
      echo "Taking densities from $inden instead in $in"
      tsd in=$inden out=$tmp.tab.den coord=t
      mv $tmp.tab  $tmp.tab.vel
      tabmath $tmp.tab.vel,$tmp.tab.den $tmp.tab %1,%2,%3,%4,%10 all
    endif
    if ($status) goto cleanup
    tabtos in=$tmp.tab out=$tmp.s0 block1=x,y,vx,vy,mass
    snaptrans in=$tmp.s0 out=$tmp.s1 ctypei=cyl ctypeo=cart
    snaprotate in=$tmp.s1 out=$tmp.s2 theta=180 order=z
    snapadd $tmp.s1,$tmp.s2 $tmp.s3

    # project for skyview, and create a intensity and velocity field

    snaprotate $tmp.s3 $tmp.snap \
        "atand(tand($phi)/cosd($inc)),$inc,$pa" zyz

    echo -n "Projected model velocities:"
    snapprint $tmp.snap -vz | tabhist - tab=t |& grep min

    foreach mom (0 1 2)
         snapgrid in=$tmp.snap out=$tmp.$mom \
                xrange=$range yrange=$range nx=$n ny=$n moment=$mom mean=t
         ccdsmooth in=$tmp.$mom out=$tmp.$mom.s gauss=$beam
    end
    ccdmath $tmp.1.s,$tmp.0.s $tmp.vel %1/%2
    ## BUG: ifgt() doesn't work
    ##    ccdmath $tmp.0.s - "ifgt(%1,0,log(%1),-10)" | ccdfits - $out.den.fits
    ccdmath $tmp.0.s - "log(%1)" | ccdmath - - "ifeq(%1,0,-10,%1)" |\
        ccdfits - $out.den.fits  \
        object=$in comment="$0 $in $out $pa,$inc,$phi,$range,$n,$beam" $wcspars 
    ccdmath $tmp.2.s,$tmp.0.s,$tmp.vel - "sqrt(%1/%2-%3*%3)" |\
        ccdfits - $out.sig.fits \
        object=$in comment="$0 $in $out $pa,$inc,$phi,$range,$n,$beam" $wcspars 
    ccdfits $tmp.vel $out.vel.fits \
        object=$in comment="$0 $in $out $pa,$inc,$phi,$range,$n,$beam" $wcspars 

    # create the (smoothed) cube;    can be quite memory intensive

    if ($cube) then
      if ($denlog) then
        snapgrid in=$tmp.snap out=- \
          xrange=$range yrange=$range zrange=-${vmax}:${vmax} \
	  xvar=x yvar=y zvar=-vz \
	  nx=$n ny=$n nz=$nvel moment=0 mean=t |\
        ccdsmooth - - gauss=$beam |\
        ccdmath - - "log(%1)" |\
        ccdmath - - "ifeq(%1,0,-10,%1)" |\
        ccdfits - $out.cube.fits \
          object=$in comment="$0 $in $out $pa,$inc,$phi,$range,$n,$beam" $wcspars 
      else
        snapgrid in=$tmp.snap out=- \
          xrange=$range yrange=$range zrange=-${vmax}:${vmax} \
	  xvar=x yvar=y zvar=-vz \
	  nx=$n ny=$n nz=$nvel moment=0 mean=t |\
        ccdsmooth - - gauss=$beam |\
        ccdfits - $out.cube.fits \
          object=$in comment="$0 $in $out $pa,$inc,$phi,$range,$n,$beam" $wcspars 
      endif
    endif

    rm -fr $tmp.*
else
    echo Warning: skipping gridding and projecting
endif

exit 0

cleanup:
    echo Some error occured
    if ($clean) rm $tmp.*



