.TH SMOOTH 2.01 "11 OCTOBER 1994" "Stadel Release 2.01" "SMOOTH"
.SH NAME
\fBsmooth\fP \- Calculate smoothed quantities over all particles
.SH SYNOPSIS
.B smooth
[\-s <nSmooth>[dgs]] [\-b <nBucket>] [\-g]
.br
 	[\-o <Output Name>] [\-p <xyzPeriod>] 
.br
	[\-px <xPeriod>] [\-py <yPeriod>] [\-pz <zPeriod>]
.br
	[\-do <MarkFile>]
.br
	[density] [meanvel] [speed] [veldisp] [phase] [mach]
.br
	[all] [null]
.br
.br
Reads particles from stdin in the TIPSY BINARY format.

.SH OPTIONS
.B [\-s <nSmooth>[dgs]]:
Default nSmooth = 64 (dgs)
.br
Specifies the number of neighbors to smooth over. Directly
following the integer any of the letters dgs may follow. A 'd' 
specifies that smooth consider the dark matter in the input
file in its calculations. A 'g' sepecifies that smooth consider 
the gas in the input file and an 's' specifies that stars are to
be considered. For example \-s 64dgs (\-s 64) would smooth over
64 particles considering all the matter in the simulation
(the default). Specifying \-s 32gs would consider only the gas
and the stars and ignore any dark matter in the input file.

.B [\-b <nBucket>]:
Default nBucket = 16
.br
Specifies the number of particles in the leaves of the search
tree. This number can be changed to improve performance, but does not
influence the calculation.

.B [\-g]:
specifies that only the 'gather' kernel be used in calculating
mean quantities. The default is to use the full SPH gather\-scatter
kernel.

.B [\-o <Output Name>]:
Default output name is "smooth".
.br
This argument allows the user to specify an alternate output
name to be used for all the smooth outputs. For example -o dark would
produce dark.den, dark.mvl, dark.dsp, etc.
 
.B [\-p <xyzPeriod>]:
Default non-periodic in x, y and z dimensions
.br
Specifies that smoothing take into account periodicity in the
all three dimensions given by the same period <xyzPeriod>. This is a
shorthand to specifying all three parameters given below.

.B [\-px <xPeriod>]:
Default non-periodic in x dimension
.br
Specifies that smoothing take into account periodicity in the
x\-dimension given by the period <xPeriod>.

.B [\-py <yPeriod>] [\-pz <zPeriod>]:
As above.

.B [\-do <Markfile>]
.br
Specifies that smooth should calculate mean field quantities only
for those particles "marked" in the <MarkFile>. The MarkFile is in 
TIPSY ARRAY format, where a 1 signifies a marked particle and a 
0 an unmarked particle (unmarked particles not computed). Whereas
the (dgs) specification of \-s (see above) ignores COMPLETELY particle
types not selected, this takes into account the EFFECT of unmarked 
particles. The values of quantities for unmarked particles are 
simply zero (or null vector in the case of meanvel). The kernel used
is automatically the GATHER-ONLY kernel (\-g)! This option is useful
when, for example, a fast density calculation is required for a small
subset of particles in a large simulation. Such a subset could be 
found with a group finder (see FOF and DENMAX) or from a minimum
density cut.

.B [density]:
produces a density output file (smooth.den) in tipsy
ascii array format

.B [meanvel]:
produces a mean velocity output file (smooth.mvl) in
tipsy ascii vector format

.B [speed]:
produces a mean speed output file (smooth.spd) in
tipsy array format

.B [veldisp]:
produces a velocity dispersion output file (smooth.dsp)
in tipsy array format

.B [mach]:
produces a mach number output file (smooth.mch) in
tipsy array format

.B [phase]:
produces a phase density output file (smooth.phs) in
tipsy array format

.B [all]:
produces all of the above output files

.B [null]:
produces no output files but searches for nearest neighbors and then
gathers the nearest neighbors (Smoothing and Resmoothing). This is useful for
gauging the performance of smooth.

.SH DESCRIPTION
Smooth can calculate several mean quantities for all particles in an 
N-Body simulation output file. The program produces a file for each
type of output specified on the command line. This output file is in
ascii format with one smoothed quantity for each particle. The program
uses a symmetric SPH (Smoothed Particle Hydrodynamics) smoothing kernel
to find the mean quantities.

	<Ai> = SUMj(Aj*Mj/Dj*W(i,j)) over all j

where <Ai> is the smoothed value for particle i of quantity A,
Mj is the Mass of particle j, Dj is the density of particle j and W(i,j)
is the symmetric smoothing kernel defined by,

	W(i,j) = 1/2(w(Rij,Hi) + w(Rij,Hj)) and,

Rij is the distance from particle i to particle j.
.br
Hi is HALF the distance to the k'th most distant neighbor of particle i.

	w(R,H) = 1/(PI*H^3)*(1 - 0.75*(2 - R/H)*(R/H)^2
.br
		if (R/H) < 1 and,

	w(R,H) = 1/(PI*H^3)*0.25*(2 - R/H)^3
.br   
		if (R/H) >= 1.

The density for particle i is simply,
.br
	Di = SUMj(Mj*W(i,j))
.br

.br
Quantities which smooth will calculate are density, mean velocity, mean
speed, velocity dispersion, mach number and phase density. For simulations
which have periodic boundary conditions, the periods in x, y, and z must
be specified using the \-px, \-py and \-pz arguments for smooth to produce
the correct results at the boundaries. As a shorthand, if the periods in 
the three dimensions are the same (usual) the \-p argument suffices.

.SH EXAMPLES

>
.B smooth density < tipsy.bin

This produces a density output file called smooth.den from the tipsy binary
file tipsy.bin. The program has smoothed density over 64 (default) nearest
neighbors and over all the particles in the tipsy.bin file.

>
.B smooth \-s 32g -o tipsy_gas density veldisp mach < tipsy.bin

This produces a density output file, tipsy_gas.den, a velocity dispersion
output file, tipsy_gas.dsp and a mach number output file, tipsy_gas.mch.
Smoothing has been done over 32 nearest neighbors and taking into account
only the gas particles.

>
.B smooth \-s 64d \-px 256 \-py 256 \-pz 256 density < tipsy.bin
(or)

>
.B smooth \-s 64d \-p 256 density < tipsy.bin

This produces a density output file, smooth.den, which has been smoothed
over a periodic cube of 256 units to a side. Only the dark matter has 
been considered in calculating the density.

.SH NOTE ABOUT OUTPUT

It is important to note that the output arrays always correspond to the
input file. For example if an input file is given with 32000 dark matter,
32000 gas and 1000 star particles the output tipsy array's will have
65000 entries, even if only the gas particles were smoothed over. The 
particles which were ignored have a value of zero in the array or
vector format.

.SH WARNINGS

When using periodic boundary conditions assure that the periodicity is 
identical to that used during the N-Body simulation which generated the
input file being used.

Smooth output files can be quite large since they are in ascii format.
Make sure that sufficient disk space is available. When running smooth
on several input files make sure to specify unique output names for the
files otherwise these files will simply be over\-written.

Smoothing time scales slightly worse than linear with number of particles
and linear with nSmooth (number of nearest neighbors).

.SH BUGS
  
Bugs? Not!

.SH SEE ALSO
.B tipsy(1), denmax(1), fof(1)



