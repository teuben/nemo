.TH FOF 1 "July 28, 1994" "Mumey Release 1.0" "AFOF"
.SH NAME
\fBafof\fP \- Find groups in N\-body simulations using an approximate
friends\-of\-friends method.
.SH SYNOPSIS
.B fof
\-e <Linking Length> \-a <Approximate Factor> [\-m <nMinMembers>] [-dgs] [-v]
.br
	[\-o <Output Name>] [\-p <xyzPeriod>]
.br
	[\-px <xPeriod>] [\-py <yPeriod>] [\-pz <zPeriod>]
.br

Reads particles from stdin in the tipsy binary format.

.SH OPTIONS
.B \-e <Linking Length>:
No Default tau value!
.br
This argument specifies the linking length used by the
friends\-of\-friends method.

.B \-a <Approximation>
No Default epsilon value!
.br
This argument specifies amount of error, epsilon, tolerated.
That is, points up to ( 1 + epsilon ) tau apart may be linked.

.B [\-m <nMinMembers>]:
Default nMinMembers = 8
.br
This argument allows the user to specify what the minimum
number of members (particles) per group should be. All groups 
with less than this value will be dicarded.

.B [\-dgs]:
Default: all particles considered (-dgs)
.br
This set of flags allows the user to specify the particle
types to be considered for grouping. If the user specifies \-d, then
only the dark matter particles are considered. With \-dg, both dark
and star particles are considered but any star particles are 
ignored. The default is, effectively \-dgs.

.B [\-v]:
Default is no output
.br
This flag allows the user to enable diagnostic output
(on stdout). This includes the number of groups	found by the
friends\-of\-friends method, the number of groups meeting the
minimum members criterion and also the CPU time needed to find
the groups.

.B [\-o <Output Name>]:
Default Output Name = "afof".
.br
This allows the user to specify an alternate name for the 
output files. For example, if the user specifies -o sim5 on the
command line, then fof will by default produce the file sim5.grp.

.B [\-px <xPeriod>]:
Default non-periodic in x dimension
.br
Specifies that afof take into account periodicity in the
x\-dimension given by the period <xPeriod>.

.B [\-py <yPeriod>] [\-pz <zPeriod>]:
As above.

.B [\-p <xyzPeriod>]:
Default non-periodic in x,y and z dimensions
.br
Specifies that afof take into account periodicity in all
three dimension given by the period <xyzPeriod>. This can be used
in place of the above as a shorthand.

.SH DESCRIPTION

The afof program finds groups in N\-body simulations using an approximate
friends\-of\-friends method. It rejects any groups found which
have less than nMinMembers number of particles. A
friends\-of\-friends group is one in
which every particle in the group has
a "friend" particle within a distance less than or equal to 
some specified "linking length". The linking length is set
with the \-e <Linking Length> argument.

.SH OUTPUT FILE

.B afof.grp
: This ASCII file is in TIPSY ARRAY format 
and contains the group number to which each particle in 
the input file belongs. Group number zero means this 
particle was not grouped. This file can be read in by tipsy or
any other analysis tool able to read this format.

.SH EXAMPLES

> afof \-e 1.0e\-4 \-a 0.1 \-p 1 \-o dark < dark.bin

This example groups all particles in the file dark.bin, with at most a 10%
error on linking lengths.
The simulation is periodic with period length of 1.0 in
each dimension. The minimum group members is taken to be the
default value of 8 and the name of the output file is 
specified to be dark.grp. 

> afof \-e 1.0 \-a 0.0 \-v \-d \-m 4 \< gasrun.bin

This example groups only the dark matter particles in the input
file gasrun.bin.  A epsilon setting of 0.0 means that the exact solution
will be found.  The simulation state is treated as non-periodic,
since there are no periods specified. Verbose (\-v) ouput is
requested. The number of minimum group members has been reduced
to 4 and the name of the output file will be fof.grp.

.SH BUGS

Please report any!
  
.SH SEE ALSO

fof(1), tipsy(1), smooth(1), denmax(1)















