.TH FOF 1.1 "Oct 26, 1994" "Stadel Release 1.1" "FOF"
.SH NAME
\fBfof\fP \- Find groups in N\-body simulations using the friends\-of\-friends method.
.SH SYNOPSIS
.B fof
\-e <Linking Length> [\-m <nMinMembers>] [-dgs] [-v]
.br
	[\-o <Output Name>] [\-p <xyzPeriod>]
.br
	[\-px <xPeriod>] [\-py <yPeriod>] [\-pz <zPeriod>]
.br

Reads particles from stdin in the tipsy binary format.

.SH OPTIONS
.B \-e <Linking Length>:
No Default fEpsilon value!
.br
This argument specifies the linking length used by the
friends\-of\-friends method.

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
Default Output Name = "fof".
.br
This allows the user to specify an alternate name for the 
output files. For example, if the user specifies -o sim5 on the
command line, then fof will by default produce the file sim5.grp.

.B [\-px <xPeriod>]:
Default non-periodic in x dimension
.br
Specifies that fof take into account periodicity in the
x\-dimension given by the period <xPeriod>.

.B [\-py <yPeriod>] [\-pz <zPeriod>]:
As above.

.B [\-p <xyzPeriod>]:
Default non-periodic in x,y and z dimensions
.br
Specifies that fof take into account periodicity in all
three dimension given by the period <xyzPeriod>. This can be used
in place of the above as a shorthand.

.SH DESCRIPTION

The fof program finds groups in N\-body simulations using the
friends\-of\-friends method. It rejects any groups found which
have less than nMinMembers number of particles. A
friends\-of\-friends group is one in
which every particle in the group has
a "friend" particle within a distance less than or equal to 
some specified "linking length". The linking length is set
with the \-e <Linking Length> argument.

.SH OUTPUT FILE

.B fof.grp
: This ASCII file is in TIPSY ARRAY format 
and contains the group number to which each particle in 
the input file belongs. Group number zero means this 
particle was not grouped. This file can be read in by tipsy or
any other analysis tool able to read this format.

.SH EXAMPLES

> fof \-e 1.0e\-4 \-p 1 \-o dark < dark.bin

This example groups all particles in the file dark.bin.
The simulation is periodic with period length of 1.0 in
each dimension. The minimum group members is taken to be the
default value of 8 and the name of the output file is 
specified to be dark.grp. 

> fof \-e 1.0 \-v \-d \-m 4 \< gasrun.bin

This example groups only the dark matter particles in the input
file gasrun.bin. The simulation state is treated as non-periodic,
since there are no periods specified. Verbose (\-v) ouput is
requested. The number of minimum group members has been reduced
to 4 and the name of the output file will be fof.grp.

.SH WARNINGS

If the linking length is set unreasonably large or if a large number
of particles are found within a ball of radius the linking length, then
the fof program will take VERY long to complete. The reason for this
is that in such cases the algorithm used here degenerates to O(N^2).
For typically desired linking lengths the fof program completes in a
reasonable amount of time. Use AFOF (approximate friends\-of\-friends)
when group finding with large linking lengths. The algorithm for
AFOF does not suffer from this degeneracy and is an O(N) (linear scaling) 
method when some user defined fuzziness in the linking length can be
tolerated.

.SH BUGS

Please report any!
  
.SH SEE ALSO

afof(1), denmax(1), tipsy(1)















