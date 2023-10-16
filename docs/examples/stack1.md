# Stacking

Stacking is a technique  to re-align spectra (or places where you might expect spectra) in order to
make the spectra show up with higher signal to noise. A good, but very straightforward example is
single dish observations, where we don't need to shift the spectra according to a model. This is 
actually the application we want to show here. However, first we show an example where we do not have
to shift the spectra:   Single Dish Observations.

The NEMO commands listed in this document have 3 ways of help:

1. Online man page on https://teuben.github.io/nemo/man_html/index1.html
2. The command "man XXX", if you have NEMO installed
3. THe command "XXX help=h", if you have NEMO installed

## Single Dish Observing

We will take 100 observations, with a line
where the S/N is 0.25, so there is no way the line can be seen in a single spectrum. We do this by 
creating a so-called Waterfall Map, where Time (observation) is along the Y axis, and the Spectrum
dimension (say velocity) along the X axis. We will then stack the spectra along the Y axis, which
in NEMO is taking a moment along the 2nd axis.

First creating the waterfall map:

	sn=0.25
	pos=40
	wid=10
	ccdmath  out=map2.ccd fie="rang(0,1)+$sn*exp(-((%x-$pos)/$wid)**2)" size=100,100 seed=123
	ccdfits  map2.ccd map2.fits
	ccdplot  map2.ccd
	ds9 map2.ccd
	
Now we will take the mean of all spectra, by taking a moment along the
2nd axis. The newly produces map2.mom image is then converted to a
simple ascii spectrum with colum 1 and 2 the "velocity" and
"intensity" resp.:

	ccdmom   map2.ccd map2.mom axis=2 mom=-1
	ccdprint map2.mom x= newline=t label=x > map2.spec
	
It can be plotted, note the 0 baseline that is plotted

	tabplot map2.spec line=1,1 ycoord=0
	
and fitted with a gaussian. We're cheating and giving it the initial conditions that we know should
exist (for low S/N as this spectrum, the initial condition estimator is not very good
	
	tabnllsqfit map2.spec fit=gauss1d par=0,$sn,$pos,$wid
	Fitting a+b*exp(-(x-c)^2/(2*d^2)):
	a= -0.00077811 0.0127007 
	b= 0.259647 0.0336685 
	c= 41.4001 1.10724
	d= 7.67932 1.22789
	rms= 0.0957405
	
If you set sn=2.5 you will see the command

	tabnllsqfit map2.spec fit=gauss1d 
	
has no problem finding the correct solution.

## Rotating Disk Galaxies

The gas in a rotating disk (galaxy, proto-planetary disk, etc.) follows some expected V(R) rotation, to
first order. Programs like **rotcur** can use a tilted ring model and find the V(R) in a set of rings.
The program **ccdzshift** will need a model (e.g. a velocity field) and then shift the spectra so they
are all aligned. They could then be stacked. An example could be weak 13-CO emission, but you do have
strong 12-CO.  This way you might be able to get a more reliably 13CO/12CO.  Another application is
the density profile of the gas. As you align the spectra, using the programm **ccdellint** you can
effectively stack the spectra along an ellipse and increase S/N that way.

## Outflows

Related to rotating disks is the potential outflow along the minor axis. Here we expect the outflow
along the near-side of the minor axis to have receding velocities. The program **rvstack** takes advantage
of this

ccdvel: make a velocity field
velcube: make a cube


nemoinp 0:100:5 > map1.radt
 nemoinp 0:100:5 | tabmath - - "exp(-%1/20)" all > map1.dent
 
 
 
nemoinp 0:100:1 > map1.radt
tabmath map1.radt - "%1/sqrt(40+%1*%1)" all > map1.velt


ccdvel out=map1.den radii=@map1.radt vrot=@map1.dent pa=30 inc=60 amp=t

ccdvel out=map1.vel radii=@map1.radt vrot=@map1.velt pa=30 inc=60 size=256
velcube out=map1.cube invel=map1.vel
