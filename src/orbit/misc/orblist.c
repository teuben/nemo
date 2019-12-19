/*
 * ORBLIST:  list an orbit path
 *
 *	16-jul-87  V1.1  minor expansion
 *	28-jul-87  V2.0  new orbit(5) structure
 *	 2-jun-87  V2.1  small mod, alsao new filestruct	PJT
 *	22-may-90  V2.2  minor improvement for new getparam()	PJT
 *	24-jul-92  V2.3  new NEMO - etc. PJT
 *      15-feb-03  V2.4  added format=				PJT
 *    10-dec-2019  V2.5  deal with PHI/ACC                      PJT
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
  "format=%g\n        Format for output",
  "mode=t,p,v\n       Output:  o,t,p,v (not used)",
  "VERSION=2.5\n      10-dec-2019 PJT",
  NULL,
};

string usage = "list an orbit path";

#ifndef HUGE
# define HUGE 1.0e20
#endif
#define EPS  0.00001

local string	infile;			/* file names */
local stream  instr;			/* file streams */

local orbitptr optr;

local int    n;
local double trange[2];   /* time range to plot */
local int maxsteps;

void nemo_main(void)
{
	int ndim;
	string format = getparam("format");

	infile = getparam("in");
	n = getiparam("n");
        maxsteps = getiparam("maxsteps");

	instr = stropen (infile,"r");

	optr = NULL;
#if 0
        allocate_orbit(&optr,ndim,maxsteps);
#endif

	while (read_orbit(instr,&optr)) {
	    list_orbit(optr, -HUGE,  HUGE, n, format);
            Nsteps(optr) = maxsteps;        /* reset */
        }
	strclose(instr);
}

