/*
 * ORBMAX:	find places where orbit has maximum radial extent
 *
 *      29-jul-09 1.0   created
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <orbit.h>

string defv[] = {
    "in=???\n		Input orbit",
    "mode=1\n           1: show max    2: show min     3: show both",
    "first=false\n      Also try and figure out if first point was min or max",
    "VERSION=1.1\n  	30-jul-09 PJT",
    NULL,
};

string usage="Tabulate radial maxima (or minima) along an orbit";

string cvsid="$Id$";

string	infile;			/* file name */
stream  instr;			/* file stream */

orbitptr optr;

void max_orbit(orbitptr o, int mode);

nemo_main ()
{
    int mode;

    infile = getparam("in");
    instr = stropen (infile,"r");
    mode = getiparam("mode");

    optr=NULL;
    while (read_orbit (instr,&optr)) 
      max_orbit(optr, mode);
    strclose(instr);
}



void max_orbit(orbitptr o, int mode)
{
  real x,y,s,t, x0,y0,s0,t0, x1,y1,t1, r1;
  int i;
  permanent bool first = TRUE;
    
  for (i=0; i<Nsteps(o); i++) {
    x = Xorb(o,i);
    y = Yorb(o,i);
    s = x*Uorb(o,i) + y*Vorb(o,i);
    t = Torb(o,i);
    if (!first) {
      if (s*s0 < 0) {    /* sign changed, we have an extremum */
	if (  (s0>0 && mode&0x01) || (s0<0 && mode&0x02) ) {
	x1 = x0 + (x-x0)*s0/(s-s0);
	y1 = y0 + (y-y0)*s0/(s-s0);
	t1 = t0 + (t-t0)*s0/(s-s0);
	r1 = sqrt(x1*x1+y1*y1);
	printf("%g  %g    %g %g %g\n",t1,atan2(y1,x1),x1,y1,r1);
	}
      }
    } else
      first = FALSE;
    x0 = x;
    y0 = y;
    s0 = s;
    t0 = t;
  }
}
