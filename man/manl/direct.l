.TH DIRECT 1.2 "Jan. 18, 1995" "Stadel Release 1.2" "DIRECT"
.SH NAME
.B direct
\- Calculates the gravitational acceleration and potential for all particles
in an N\-body simulation output using the direct, O(N^2), method.
.SH SYNOPSIS
.B direct
[\-e <fSoftRadius>]  [\-spline] [\-plummer] [\-uniform] 
.br
	[\-p <xyzPeriod>] [\-G <fGravConst>] [\-b <iBlockSize>] 
.br 
	[\-dgs] [\-do <MarkFile>] [\-o <FileName>]
.br
	[\-i <iChkptInterval>] [\-restart] [\-v] [\-t]
.br
	

Reads particles from stdin in the tipsy binary format.

.SH OPTIONS
.B [\-e <fSoftRadius>]:
Default: takes softening from input file
.br
This argument sets the softening length used in calculating the
gravitational interactions for all particles.

.B [\-spline]:
This is the default.
.br
Specifies that cubic spline softening is to be used. The equations for
the interaction are those specified by equations (A1) and (A2) from
the paper TREESPH: A UNIFICATION OF SPH WITH THE HIERARCHICAL TREE
METHOD, Lars Hernquist and Neal Katz,
APJ Supplement Series 70:416-446, 1989.

.B [\-plummer]:
.br
Specifies that plummer softening is to be used. The interaction is given
by phi = -m/(|r|^2 + h^2)^(1/2) and a = -m*r/(|r|^2 + h^2)^(3/2), where h
is the softening length and r is the vector from particle j to particle i. 
Note: this type of softening is NOT newtonian beyond some distance as
are the other softening methods. This can put an otherwise equilibrium
model slightly out of equilibrium.

.B [\-uniform]:
.br
Specifies that uniform\-density sphere softening is to be used. The 
interaction is given by phi = -m*|r|^2/h^3, a = -m*r/h^3  for (|r| < h)
and phi = -m/|r| and a = -m*r/|r|^3 for (|r| >= h) where h is the
softening length and r is the vector from particle j to particle i.

.B [\-p <xyzPeriod>]:
Default: non-periodic calculation
.br
This argument specifies that periodic boundary conditions are to be
considered in calculating gravity. The argument value sets the size 
(length of each side) of the periodic box which contains the particles.
The Ewald summation method is used when calculating the gravitational
interactions.

.B [\-G <fGravConst>]:
Default fGravConst = 1.0
.br
This argument allows the user to specify the gravitational constant
used by direct to calculate the gravitational interactions.

.B [\-b <iBlockSize>]:
Default iBlockSize = 256
.br
This argument optimizes direct's use of the primary cache on the machine.
It should be set to reflect the size of your machine's primary data
cache. The default value assumes a primary data cache size of 16k bytes.
If you know the size of your primary data cache then the formula is
.B iBlockSize
= nPrimaryDataBytes / 64. You can also try different
.B iBlockSize
values on a small data\-set, to get the optimum value.

.B [\-dgs]:
Default: all particles considered (-dgs)
.br
This set of flags allows the user to specify the particle types to be
considered for the gravity calculation. If the user specifies
.B \-d,
then only the dark matter particles are considered. With
.B \-dg,
both dark and gas particles are considered but any star particles are 
ignored. The default is, effectively
.B \-dgs.

.B [\-do <MarkFile>]:
Default action: all
.B [\-dgs]
selected particles are calculated

This flag allows only marked particles to be calculated.
The effects of the other selected particles are still taken into
account. The mark file is in TIPSY ARRAY format, where an entry
of NOT\-zero means that the particle is marked for calculation. A zero
entry means the particle's acceleration and potential will not be
calculated, but its effects on any marked particles will be 
included. The accelerations and potentials of all the unmarked
particles are output as zero.

.B [\-o <FileName>]:
Default Output Name = "direct"
.br
This allows the user to specify an alternate name for the 
output files. For example, if the user specifies
.B \-o king
on the command line, then direct will by produce the files
king.acc and king.pot.

.B [\-i <iChkptInterval>]:
Default is no checkpoint files
.br
This argument enables the generation of checkpoint files. The
interval between successive checkpoints is given in number of
blocks by
.B <iChkptInterval>.
This is useful when large data\-sets are being generated to allow
recovery after a system crash or other problem. The 
.B \-restart
flag allows restarting the calculation from the last checkpoint.

.B [\-restart]:
Default is no restart
.br
This flag enables restart mode. All command line arguments are 
ignored other than 
.B \-o <FileName>, \-i <iChkptInterval>
and
.B \-v.
If
.B \-o <FileName>
is specified, direct will restart the calculation from the checkpoint
file <FileName>.chk, otherwise it looks for the file direct.chk and
restarts from that. Note: if you do not specify a new <iChkptInterval>,
then direct will not remove the checkpoint file at the end of the
calculation.

.B [\-v]:
Default is no output
.br
This flag allows the user to enable diagnostic output (on stdout). 
The output produced is the block indecies which direct is currently
processing. The time between outputs should be roughly constant. This
will make you feel better when direct is calculating a lagre data\-set,
where the time between start and any real output generated is large.
It will also output the CPU time used by direct a the end.

.B [\-t]:
Default is test\-mode not enabled
.br
This flag is a developer's flag and simply does a "no\-brainer" N^2
algorithm to test the other features of this tool. Checkpointing is
not enabled when this flag is specified.

.SH DESCRIPTION

.B Direct
is a tool designed to test more complicated gravity codes such as
Barnes\-Hut, Fast Multipole and error\-controlling tree\-codes as well
as FFT codes. It is also the method of choice when dealing with less than 
about 10000 particles. The primary application remains the testing and
error analysis of more sophisticated codes. Typical things looked at 
are the distrbution of relative and absolute errors in the acceleration
as well as the maximum and rms values. Another application when a well
understood mass distribution is calculated (such as a king model), is the
analysis of the discreteness noise as a function of N.

The method used by direct to calculate the gravitational accelerations and
potentials is to process the interaction matrix (an NxN matrix with 
entry_ij being the interaction between particle i and particle j) in blocks.
The interaction matrix is also anti\-symmetric with respect to the
accelerations (a_ij + a_ji = 0, i!=j) and symmetric with respect to 
potentials (p_ij - p_ji = 0, i!=j). This allows only the lower triangle of
blocks to be calculated, with special treatment of the diagonal blocks.
This makes direct quite a bit faster than the simplest implementation,
and when the block size is set properly
(see
.B \-b <iBlockSize>
above) makes it scale exactly as O(N^2) in CPU time.

The softening is now "symmetric", such that the softening used for an
interaction is taken to be the average of the softenings of the two
interacting particles, h_ij = 0.5*(h_i + h_j). For the case where all
particles have the same softening, the results are of course the same. 

.SH OUTPUT FILES

.B direct.acc
: This ASCII file is in TIPSY VECTOR format 
and contains the acceleration vectors of all particles in 
the input file. If certain types of particles were omitted with the
.B -dgs
flag combination, then the acceleration in this file will be zero
for those particles. This file can be read in by tipsy or
any other analysis tool able to read this format.

.B direct.pot
: This ASCII file is in TIPSY ARRAY format 
and contains the potentials of all particles in 
the input file. If certain types of particles were omitted with the
.B -dgs
flag combination, then the potential in this file will be zero
for those particles. This file can be read in by tipsy or
any other analysis tool able to read this format.

.SH EXAMPLES

> direct \-e 0.1 \-o dark < dark.bin

This calculates the gravitational acceleration and potential for all
particles in the tipsy binary input file dark.bin. The softening radius
is set to 0.1 and the masses for the particles are taken from the 
input file. The gravitational constant is the default G = 1.0.
The files dark.acc and dark.pot are produced.

> direct \-p 1 \-o dark_p < dark.bin

This calculates the gravitational acceleration and potential for all
particles in the tipsy binary input file dark.bin. The softening radii
and masses for the particles are taken from the input file. The calculation
will take into consideration periodic boundary conditions with a periodic
box of length 1.0. The Ewald summation technique is used to include the
effects of periodic boundary conditions.
The gravitational constant is the default G = 1.0.
The files dark_p.acc and dark_p.pot are produced.

> direct \-e 0.025 \-v \-d < gasrun.bin

This calculates the gravitational acceleration and potential for only
the dark matter particles of the input file gasrun.bin. The accelerations
and potentials of the gas particles will be set to zero in the output 
files. Verbose output is requested, and will cause direct to output 
diagnostic output on stdout. Softening radius of
0.025 is used and the gravitational constant is again G = 1.0.
The default files direct.acc and direct.pot are produced.

> direct \-o big \-v \-i 10000 \-restart

This restarts a previous calculation from the checkpoint file, big.chk.
Verbose output is requested and a new checkpoint interval of 10000
is specified (this is likely the same as the previously used value).

.SH WARNINGS

1) Since
.B direct
is an O(N^2) method (and scales precisely so), it can take VERY long
to calculate large problems. For a million particles, a week to several
weeks of CPU time may be required.

2) If you are using periodic boundary conditions EXPECT the calculation
to take about 100 times longer! This increase in time is due to the 
very expensive Ewald interactions. Note: faster, more approximate 
techniques exist for handling periodic BCs, but this code was designed
to provide a correct answer.

3) Do NOT modify the doubles in the code to floats. Although this 
provides significant speedup, the results may not be useful due to 
round\-off errors.

4) Round\-off and truncation errors are an aspect where a high
precision multipole calculation (such as FMM code) may well produce
more accurate results. Be careful in using direct to test such codes!

.SH BUGS

Please report any!
  
.SH SEE ALSO

tipsy(1), smooth(1)















