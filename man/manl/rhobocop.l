.TH RhoBoCoP l "10 December 1994" 
.SH NAME 
rhobocop - density bootstrap computing program 
.SH SYNOPSIS
\fBrhobocop <\fP \fIparameter_file\fP

.SH DESCRIPTION 
\fBrhobocop\fP (Density BOotstrap COmputing Program) is a   standard F77
program  to evaluate several quantities of interest in N-body systems.
It gives you the density and velocity profiles in spherical shells
centered on the density maximum of the system considered. For each
measured quantity it gives the associated error due to finite
N-sampling. This is made employing standard bootstrap method (see Efron
& Tibshirani:'An Introduction to the Bootstrap', Chapman & Hall, 1993;
Heyl, Hernquist & Spergel Apj, 427,165).   Multi system data (e.g.
cosmological simulations containing many separate halos) can be
examined.
.PP
Comments and feedback are welcome! (see AUTHOR below)
.PP
Input file must be in SNAP205 ASCII or TIPSY ASCII format.
Output files are in ASCII.
.SH INPUT PARAMETERS
The input parameters can be entered interactively, as the program
reads from standard input. You can also prepare an input parameter file,
and run the program in 'batch' mode as given on the SYNOPSIS line. 
The input parameters are, in order:
.TP 20
\fBInput file:\fP
name of input file
.TP
\fBRoot name:\fP
4 chars identifying the root of output files' names
.TP
\fBFirst & last:\fP
first and last indexes of particles in the input file considered by RhoBoCoP
.TP
\fBmulti system option:\fP
\fB1\fP if you need to translate to a particular region. If you use this,
two new keywords are introduced:

\fBx,y,z\fP: coordinates of the subsystem center

\fBradius:\fP cutoff radius for particles to be examined (i.e. the radius
if the individual system under study
.TP
\fBnmin,facr & radius:\fP
the center of the system under study is found iteratively
removing the external particles starting at "radius" and
obtaining the center of mass at each step reducing the
radius by a factor "facr" (<1 !!) until the number of
particles remaining  in the sphere considered is less then
or equal to "nmin".

nmin best value is strongly dependent from the steepness
of  the density profile, for homogeneous systems is recommanded 
to use nmin near the total number of particles, otherwise
the center of density will be chosen in an almost random way!!
.TP
\fBparticles per bin:\fP
number of particles in each shell (min number is really
problem dependent, but usually > 32)
.TP
\fBNumber of resamplings:\fP
number of reboot cycles, usually a number comparable
to the number of particles for individual shells
.TP
\fBRandom seed:\fP
a (big) integer used to initialize pseudo-random generator.
.PP
An \fBexample paramter file\fP:
.nf
.ta +2i
king       	input file name
1        	input file name type
king      	root file name
1 4096      	first and last
0          	multi-system option
256,0.95,3	nmin, facr & radius
128     	particles per bin
64      	number of resamplings
1234567890	random seed
.fi
.SH OUTPUT FILES
Several (ASCII) output files are created  (xxxx is the root input filename)
.nf
.ta +1i +0.5i +0.5i +0.5i +0.5i +0.5i +0.5i +0.5i +0.5i +0.5i +0.5i 

xxxxrho		r	r^1/4	Log(r)	rho(r)	errrho(r)	m(r)	errm(r)
xxxxrms		r	vrmsr	err	Vrmsth	err	Vrmsph	err
xxxxvel   	r	<vr>	err  	<vtheta>	err	<vphi>	err
xxxxsigma	r	sigmar	err	sigmat	err	sigmaph	err	beta	err


r        	= mean radius of particles belonging to each shell
rho(r)   	= density in each shell
errrho(r)	= 1-sigma error on density of each shell
m(r)     	= cumulative mass function
errm(r)  	= 1-sigma error on m(r)

vrmsr    	= root mean square of radial velocity at radius r
vrmsth   	= root mean square of theta velocity  at radius r
vrmsph   	= root mean square of phi   velocity  at radius r
err      	= 1-sigma error of previous estimated quantity
<vr>     	= mean radial velocity  at radius r
<vtheta> 	= mean theta  velocity  at radius r
<vphi>   	= mean phi    velocity  at radius r

sigmar   	= standard deviation (sigma)  for radial velocity at radius r
sigmat   	= standard deviation (sigma)  for theta velocity at radius r
sigmaph  	= standard deviation (sigma)  for phi velocity at radius r
beta 	 	= anisotropy parameter as in Jeans Equation;
           beta=1-(sigmath**2/sigmar**2)In source file, expressions for alternative 
           formulations of anisotropy parameters are available; since a bootstrap is 
           made, beta is not simply the result of a calculation on other columns in
           file; ( f(<q>) not always equal to <f(q)>!!).
.fi
.SH FILES
.nf
.ta +2i
$NEMO/usr/governato/	(location in NEMO)
 rhobocop.f	the source code
 rhobocop.doc	this document
 king        	King model W0=7 N=4096 (SNAP205 format)
 rbpar    	parameters for analysis of king
.fi
.SH AUTHOR 
F. Leoni (contact: \fIleoni@bach.mi.astro.it\fP) & F. Governato 
.SH HISTORY 
.nf 
.ta +1i +4i
10-dec-1994   	V1.1 released version	FL&FG 
.fi 

