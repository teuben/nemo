/*
 * ORBLIST:  list an orbit path
 *
 *	16-jul-87  V1.1  minor expansion
 *	28-jul-87  V2.0  new orbit(5) structure
 *	 2-jun-87  V2.1  small mod, alsao new filestruct	PJT
 *	22-may-90  V2.2  minor improvement for new getparam()	PJT
 *	24-jul-92  V2.3  new NEMO - etc. PJT
 */

#include <stdinc.h>		/* also gets <stdio.h	*/
#include <getparam.h>
#include <vectmath.h>		/* otherwise NDIM undefined */
#include <filestruct.h>
#include <potential.h>
#include <orbit.h>

string defv[] = {
	"in=???\n           Input filename (an orbit)",
	"n=1\n              stride in time through the orbit",
        "maxsteps=10000\n   Maximum number of steps allowed",
	"VERSION=2.3\n      24-jul-92 PJT",
	NULL,
};

#define HPI  1.5702
#define RPD (3.1415/360.0)
#ifndef HUGE
# define HUGE 1.0e20
#endif
#define UNDEF HUGE
#define EPS  0.00001

string	infile;			/* file names */
stream  instr;			/* file streams */

orbitptr optr;

int    n;
double trange[2];	/* future expansion will allow timerange */
int maxsteps;

nemo_main ()
{
	int ndim;

	infile = getparam("in");
	n = getiparam("n");
        maxsteps = getiparam("maxsteps");

	instr = stropen (infile,"r");

	optr = NULL;
#if 0
        allocate_orbit(&optr,ndim,maxsteps);
#endif

	while (read_orbit(instr,&optr)) {
	    list_orbit(optr, -HUGE,  HUGE, n);
            Nsteps(optr) = maxsteps;        /* reset */
        }
	strclose(instr);

}

