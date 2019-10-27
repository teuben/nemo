#! /bin/csh -f
#
#    example of an exponential disk with varying velocity dispersion
#    as function of radius to show the effects in a x-y-vz cube
#
#    24-may-2011      Created                  Peter Teuben
#
# DEFAULT cpu: 71.540u 10.640s 1:42.98 79.8%  on nemo (core2duo T7300 @ 2.0 GHz)
#	       28.331u 3.203s 0:30.09 104.7%  on chara (i7 870 @ 2.93 GHz)


# prefix for all derived dataset names
set run=run1

# model parameters
set nbody=100000
set nmodel=100
set qt=0.01
set z0=0.01
set mode=2
set zmode=3

# viewing and gridding parameters
set inc=60
set gmax=2
set vmax=2
set nx=256
set nv=256
set sx=0.05
set sv=0.05

# parse command line to override previously set parameter defaults

foreach arg ($*)  
  set $arg
end

# --------------------------------------------------------------------------

set xrange=-${gmax}:${gmax}
set zrange=-${vmax}:${vmax}

rm -rf $run.*

# create lots of models, to beat the noise
foreach i (`seq $nmodel`)
  mkexpdisk - nbody=$nbody Qtoomre=$qt z0=$z0 mode=$mode zmode=$zmode \
              seed=-1 >> $run.1
end

# rotate and grid them into a cube
csf $run.1 - item=SnapShot |\
  snaprotate -    -      theta=$inc order=x |\
  snapgrid -      $run.3 xrange=$xrange yrange=$xrange zrange=$zrange \
	                 nx=$nx         ny=$nx         nz=$nv \
			 evar=m stack=t

# smooth the cube
ccdsmooth  $run.3 -      gauss=$sx dir=xy cut=0.001 |\
 ccdsmooth -      $run.5 gauss=$sv dir=z  cut=0.001

# take a PV slice along major axis
ccdsub     $run.5 $run.6 y=$nx/2 dummy=f


