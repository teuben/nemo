/*
 * ORBSTAT:	 tabulate orbit and some of their statistics
 *
 *	23-mar-95 V1	Created					pjt
 *	14-apr-01       header now in dprintf()
 *      21-feb-03       gather statistics on checking how well ellipsoidal orbit is
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <orbit.h>
#include <moment.h>

string defv[] = {
    "in=???\n		Input orbit",
    "ellipse=f\n        Check on how well ellipsoidal",
    "VERSION=1.1a\n  	2-feb-05 PJT",
    NULL,
};

string usage="Tabulate some statistics of orbit(s)";

string cvsid="$Id$";

string	infile;			/* file name */
stream  instr;			/* file stream */

orbitptr optr;

void stat_orbit(orbitptr o, bool Qell);

nemo_main ()
{
    string mode;
    bool pabs,vabs,Qell;

    infile = getparam("in");
    instr = stropen (infile,"r");
    Qell = getbparam("ellipse");

    optr=NULL;
    while (read_orbit (instr,&optr)) 
        stat_orbit(optr,Qell);
    strclose(instr);
}



void stat_orbit(orbitptr o, bool Qell)
{
  Moment  xm, ym, um, vm, jm, r2m, v2m;
  real jz, t, e, xmax, ymax, umax, vmax, r2, v2;
  int i;
  permanent bool first = TRUE;
    
  ini_moment(&xm,-1,0);
  ini_moment(&ym,-1,0);
  ini_moment(&um,-1,0);
  ini_moment(&vm,-1,0);
  ini_moment(&jm,2,0);
  ini_moment(&r2m,2,0);
  ini_moment(&v2m,2,0);
  
  for (i=0; i<Nsteps(o); i++) {
    accum_moment(&xm,Xorb(o,i),1.0);
    accum_moment(&ym,Yorb(o,i),1.0);
    accum_moment(&um,Uorb(o,i),1.0);
    accum_moment(&vm,Vorb(o,i),1.0);
    jz = Xorb(o,i)*Vorb(o,i) - Yorb(o,i)*Uorb(o,i);
    accum_moment(&jm,jz,1.0);
  }
  t = Torb(o,Nsteps(o)-1);
  e = I1(o);
  if (first) {
    if (Qell)
      dprintf(0,"# T\tE\tr2_min\tr2_max\tv2_min\tv2_max\tr2_mean\tr2_sigma\tv2_mean\tv2_sigma\tLmin/Lmax\n");
    else
      dprintf(0,"# T\tE\tx_max\ty_max\tu_max\tv_max\tj_mean\tj_sigma\n");
    first=FALSE;
  }
  if (Qell) {
    xmax = max_moment(&xm);
    ymax = max_moment(&ym);
    umax = max_moment(&um);
    vmax = max_moment(&vm);
    for (i=0; i<Nsteps(o); i++) {
      r2 = sqr(Xorb(o,i)/xmax) + sqr(Yorb(o,i)/ymax);
      v2 = sqr(Uorb(o,i)/umax) + sqr(Vorb(o,i)/vmax);
      accum_moment(&r2m,r2,1.0);
      accum_moment(&v2m,v2,1.0);
    }
    printf("%g %g  %g %g %g %g  %g %g %g %g  %g\n",
	   t,e,
	   min_moment(&r2m), max_moment(&r2m), min_moment(&v2m), max_moment(&v2m), 
	   mean_moment(&r2m), sigma_moment(&r2m), mean_moment(&v2m), sigma_moment(&v2m),
	   ymax*umax/(xmax*vmax));
  } else {
    printf("%g %g %g %g %g %g %g %g\n",
	   t,e,
	   max_moment(&xm), max_moment(&ym),
	   max_moment(&um), max_moment(&vm),
	   mean_moment(&jm), sigma_moment(&jm));
  }
}
