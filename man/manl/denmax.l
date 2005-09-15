.TH DENMAX 1.1 "May 28 1994" "Stadel Release 1.1" "DENMAX"
.SH NAME
\fBdenmax\fP \- Find gravitationally bound groups in N-body simulations.
.SH SYNOPSIS
.B denmax
\-e <fEpsilon> [\-s <nSmooth>]
.br
	[\-d <fMinDensity>] [\-t <fMaxTemperature>]
.br
	[\-z <fRedShift>] [\-G <fGravConst>]
.br
	[\-O <fOmega>] [\-H <fHubble>] [\-m <nMinMembers>]
.br
	[\-p <xyzPeriod>]
.br
	[\-px <xPeriod>] [\-py <yPeriod>] [\-pz <zPeriod>]
.br
	[\-c <xyzCenter>]
.br
	[\-cx <xCenter>] [\-cy <yCenter>] [\-cz <zCenter>]
.br
	[\-r] [-nu] [\-o <Output Name>]
.br

Reads particles from stdin in the tipsy binary format.

.SH OPTIONS
.B \-e <fEpsilon>:
No Default fEpsilon value!
.br
This parameter should be set to the length
scale of interest for the groups. For example in a 
gas simulation it might be set to the rough size of
galaxies. This parameter also sets the resolution
limit of the simulation and thus the resolution of
the density gradients used for sliding the particles.
For more information of this parameter's signifigance 
see the discussion below.

.B [\-s <nSmooth>]:
Default nSmooth = 64
.br
Specifies the number of neighbors to smooth over when finding 
densities.

.B [\-d <fMinDensity>]:
Default fMinDensity = 0.0
.br
This sets the minimum density for dark matter
or gas particles which are considered for grouping.

.B [\-t <fMaxTemperature>]:
Default fMaxTemperature = HUGE
.br
This sets the maximum temperature for gas 
particles which are considered for grouping.

.B [\-z <fRedShift>]:
Default fRedShift = 0.0
.br
The redshift of the input file. This is used to
determine physical coordinates from comoving ones. The
default of 0.0 specifies that the input file is the state
of the simulation at our (the current) epoch.

.B [\-G <fGravConst>]:
Default fGravConst = 1.0
.br
The relevant gravitational constant in the
units of the simulation. Most N-body simulations have
their units set such that this constant is 1.0.

.B [\-O <fOmega>]:
Default fOmega = 1.0
.br
This sets the Omega parameter for the simulation in
question. The default value of 1.0 specifies a flat universe.

.B [\-H <fHubble>]:
Default fHubble = 0.0
.br
Sets the Hubble constant for the simulation. If the
simulation was not performed in comoving coordinates then this
can be left at zero. Note that this must be set in the units
of the simulation.

.B [\-m <nMinMembers>]:
Default nMinMembers = 4
.br
This argument allows the user to specify what the minimum
number of members (particles) per group should be. All groups 
with less than this value will be dicarded.

.B [\-px <xPeriod>]:
Default non-periodic in x dimension
.br
Specifies that denmax take into account periodicity in the
x\-dimension given by the period <xPeriod>.

.B [\-py <yPeriod>] [\-pz <zPeriod>]:
As above.

.B [\-p <xyzPeriod>]:
Default non-periodic in x,y and z dimensions
.br
Specifies that denmax take into account periodicity in all
three dimension given by the period <xyzPeriod>. This can be used
in place of the above as a shorthand.

.B [\-cx <xCenter>]:
Default box center is 0.0 in x dimension
.br
Specifies the center of the periodic box when denmax
takes into account periodicity in the x\-dimension. This arguement
can be ignored for a non-periodic simulation.

.B [\-cy <yCenter>] [\-cz <zCenter>]:
As above for y and z dimensions.

.B [\-c <xyzCenter>]:
Default box center is 0.0 in x,y and z dimensions
.br
For periodic simulations specifies the center of the box in
all three dimensions. This is a shorthand to the above.

.B [\-r]:
Causes the denmax.ray output file to be produced.
.br
By default denmax does not produce the denmax.ray output
file. With this argument present it will produce this (tipsy vector
format) file.

.B [\-nu]:
Disables Unbinding (No Unbinding).
.br
When this flag is specified denmax will bypass the unbinding
procedure but still output all the normal group information. This is 
useful if special unbinding methods are required for a simulation or
if unbinding is completely undesired. The user can write a special
unbinding program which would read the denmax.grp file.

.B [\-o <Output Name>]:
Default Output Name = "denmax".
.br
This allows the user to specify an alternate name for the 
output files. For example, if the user specifies -o sim5 on the
command line, then denmax will by default produce the files sim5.grp
and sim5.gtp (and if requested sim5.ray).

.SH DESCRIPTION

Denmax finds gravitationally bound groups in N-body simulations. The 
general procedure is first to decide which particles should be grouped,
these are called the "moving" particles.
The particles considered may be of a certain type (there are three types
in the tipsy format, dark matter, gas and stars). They may have to 
satisfy a minimum density criterion and may also (in the case of gas) 
need to satisfy a maximum temperature criterion. This choice of particles
is affected by the type of input file and the users settings for the 
minimum density and maximum temperature. There are seven possible types
of tipsy input files and here is how each is handled.

.B Dark matter only:
density and density gradients are calculated from
all the particles (dark), and the moving particles are those meeting
the minimum density criterion.

.B Gas only:
density and density gradients are calculated from all the 
particles (gas), and the moving particles are those meeting both the 
minimum density and maximum temperature criteria.

.B Stars only:
density gradients are calculated from all the particles
(stars), and all the particles are moved.

.B Gas and dark matter:
density and density gradients are calculated from
the gas particles only, and the moving particles are the gas particles 
which meet the density and temperature criteria.

.B Stars and gas:
density gradients are calculated from the 
star particles only, and the moving particles are all the star
particles.

.B Stars and dark matter:
Same as stars and gas case.

.B Stars and gas and dark matter:
density and density gradients are
calculated from all the star particles and all the gas particles.
The moving particles are then the gas particles meeting the density
and temperature criteria and all the star particles.

Once the moving particles have been found they are then moved
along the initial density gradients toward regions of higher density.
The particles are successively stepped a distance of fEpsilon/4.0 in
this manner until they stay within a distance of fEpsilon over 5 steps.
This mean that they are oscillating in a local high density region.
We then group all the particles in each high density region together using
the friends-of-friends method with a linking length of fEpsilon.
Note that this also links together high density regions seperated by 
less than fEpsilon. Once all particles are localized in this way,
we proceed to "squish" them together even more to find the high density
center. Lastly we discard groups with less than nMinMembers number of 
particles.

The next stage (if the user has not specified \-nu) is to
remove particles from a group which are not gravitationally bound 
to it. We call this phase unbinding and the initial positions are
considered here, not the moved positions. There are cases for 
unbinding, one for dark matter only or star particle only inputs
(case I unbinding) and one for the other input types (case II 
unbinding).

.B Case I unbinding:

1. Calculate the potential energies of all the particles in the 
group taking into account the redshift of the simulation to get
physical distances (redshift zero applies also to a simulation
in physical coordinates).

2. The center of mass and center of mass velocity for the group
is found. The velocity relative to the center of mass velocity 
if found for all the particles in the group. This relative velocity
is converted to a physical velocity by taking into account the 
redshift and adding a term to include the Hubble flow. For a 
simulation in physical coordinates fHubble must be 0.0.
(Note that this is the default!) From this the kinetic energy with
respect to the center of mass is found for all particles in the 
group.

3. The least bound particle is found. If this particle is bound then
all particles are bound and we have finished unbinding for this
group. Otherwise we must remove it from the group, adjusting the 
potential energies and kinetic energies (the center of mass changes)
of the remaining particles and return to step 3.

.B Case II unbinding:

Here we include some of the group's environment to the potential.
For example a galaxy depends on its dark matter halo in order to 
stay together.

1. Calculate the potential energies as per case I, step 1.

1a. Include the potential contribution of all non-grouped mass 
within a 2*fEpsilon ball about the center of the high density
region of the group. In otherwords we want to include some of
the mass about the density maximum. 

2. Calculate the kinetic energies as per case I, step 2.

3. The least bound particle is found. If this particle is bound
then again we are finished with the group. Otherwise we remove it
from the group, adjusting the kinetic energies as in case I. Then 
we go back to step 3. Note that this time we do not subtract the 
potential energy contribution of this particle from the others in
the group. The potential energies remain fixed as any removed 
particle is still considered to be part of the groups environment.

After unbinding we check once more that all groups have 
at least nMinMembers number of particles. We discard those that 
don't and output the group information.

.SH OUTPUT FILES

Denmax produces three different output files, by default a 
denmax.grp file, and a denmax.gtp file. Optionally also a 
denmax.ray file.

.B denmax.grp:
This ASCII file is in TIPSY ARRAY format 
and contains the group number to which each particle in 
the input file belongs. Group number zero means this 
particle was not grouped (could have been removed by 
unbinding). This file can be read in by tipsy or
any other analysis tool able to read this format.

.B denmax.gtp:
This file is in TIPSY BINARY format (same as
the input file) and contains one star particle
to represent each group that denmax finds. Each particle
in this file has as its position, the position of the 
group's density maximum. As its velocity, the center of 
mass velocity of the group. Each star particle's mass is 
the mass of the entire group and the star_particle.eps
field contains the radial size of the group. (The time
stamp of the input file is also stored in each 
star_particle.tform field.) This file can also be read
in by tipsy or any tool accepting TIPSY BINARY format.

.B denmax.ray:
This ASCII file is in TIPSY VECTOR format
and contains a vector pointing from the initial position
of each particle to the localized (moved and squished) 
position of that particle. For particles not moved the
vector stored is null. This file can be read in by tipsy
and used to analyse the actions of denmax.

.SH EXAMPLES

> denmax \-e 9e\-4 \-H 2.8944 \-d 170 \-p 1 \-o dark < dark.bin

This example groups dark matter particles in the file
dark.bin. The simulation is periodic with period length
of 1.0 in each dimension and box centered on (0,0,0).
The Hubble constant in system units in this case is 
2.8944 (sqrt(8*Pi/3)) and the redshift defaults to 
0 (current epoch). The density cut is made at a density 
of 170 (mean density is 1 in these units). The files
dark.grp and dark.gtp are produced.

> denmax \-e 9e\-4 \-H 2.8944 \-d 170 \-t 30000 \-z 1 \-p 1 < all.bin

Here we suppose that dark matter, gas and stars are in 
the input file. A cut in the gas is made at a density
minimum of 170 and temperature maximum of 30000 
(cold dense gas) and groups of this gas and all the
stars are formed at a redshift of 1. The files denmax.grp
and denmax.gtp are created.

> denmax \-e 0.5 \-d 100 -m 8 \-r \-nu \-o phys < phys.bin

Here denmax processes a dark matter file in physical 
coordinates which is non-periodic. It makes a density
cut at 100 and finds the groups, but does not try 
to unbind any particles. It does remove groups with
less than 8 members. The output files phys.grp, phys.gtp
and phys.ray (\-r) are produced. 

.SH WARNINGS

1. Make sure you set fEpsilon to a reasonable size for your
simulation. If too small it can take a long time to converge
and may not give reasonable groups do to the limited resolution
of objects in the simulation. If too large, it will miss the 
smaller scale objects and may group together many smaller
objects to form large associations. For example the desired
groups could be galaxies in a gas and dark matter simulation,
epsilon should reflect this fact.

2. Do not run denmax with an nSmooth of less than 64 when 
finding groups in a dark matter only input file. For gas
simulations an nSmooth of 32 is acceptable.

3. The unbinding procedure is O(Ngroup^2), this means that for
very large simulations and specifically situations where groups
larger than a couple of thousand members are found the unbinding
procedure could be prohibitive. So far this has not been a real
problem, even in a 2.1 million particle simulation, however,
future versions of denmax may need to use a tree code for the 
group pootential energies. (The friends-of-friends is actually
also O(Ngroup^2) but this is not such a big problem)

.SH BUGS

Please report any!
  
.SH SEE ALSO
.B tipsy(1), smooth(1), fof(1)
.nf
astro-ph/0405445; Fast Identification of Bound Structures in Large N-body Simulations (Weller et al)
.fi
