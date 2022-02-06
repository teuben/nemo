/*
 * ORBSTAT:	 tabulate orbit and some of their statistics
 *
 *	23-mar-95 V1	Created					pjt
 *	14-apr-01       header now in dprintf()
 *      21-feb-03       gather statistics on checking how well ellipsoidal orbit is
 *      27-apr-05 2.0   add energy check
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <orbit.h>
#include <potential.h>
#include <moment.h>

string defv[] = {
    "in=???\n		Input orbit",
    "ellipse=f\n        Check on how well ellipsoidal",
    "relative=f\n       Energy check relative to recorded version",
    "tab=\n             Output table of quantities **not implemented**",
    "pot=t\n            Use potential stored with orbit?",
    "VERSION=2.1\n  	6-feb-2022 PJT",
    NULL,
};

string usage="Tabulate some statistics of orbit(s)";

string cvsid="$Id$";

string	infile;			/* file name */
stream  instr;			/* file stream */

orbitptr optr;
bool Qpot;

void stat_orbit(orbitptr o, bool Qell, bool Qrel, stream tabstr);

void nemo_main ()
{
    string mode;
    stream tabstr;
    bool pabs,vabs,Qell,Qrel;

    infile = getparam("in");
    instr = stropen (infile,"r");
    Qell = getbparam("ellipse");
    Qrel = getbparam("relative");
    Qpot = getbparam("pot");
    if (hasvalue("tab"))
      tabstr = stropen(getparam("tab"),"w");
    else 
      tabstr = 0;

    optr=NULL;
    while (read_orbit (instr,&optr)) 
      stat_orbit(optr,Qell,Qrel,tabstr);
    strclose(instr);
}



void stat_orbit(orbitptr o, bool Qell, bool Qrel, stream tabstr)
{
  Moment  xm, ym, um, vm, em, jm, r2m, v2m;
  real jz, t, e, xmax, ymax, umax, vmax, r2, v2;
  real etot, pos[NDIM], acc[NDIM], time;
  int i, ndim=3;
  permanent bool first = TRUE;
  proc pot;
    
  ini_moment(&xm,-1,0);
  ini_moment(&ym,-1,0);
  ini_moment(&um,-1,0);
  ini_moment(&vm,-1,0);
  ini_moment(&em,2,0);
  ini_moment(&jm,2,0);
  ini_moment(&r2m,2,0);
  ini_moment(&v2m,2,0);
  pot=get_potential(PotName(o), PotPars(o), PotFile(o));
  if (pot == NULL && Qpot==TRUE) {
    warning("No potential found, using internal values");
    Qpot = FALSE;
  }
  
  t = Torb(o,Nsteps(o)-1);
  e = I1(o);
  dprintf(0,"E=%g\n",e);
  for (i=0; i<Nsteps(o); i++) {
    accum_moment(&xm,Xorb(o,i),1.0);
    accum_moment(&ym,Yorb(o,i),1.0);
    accum_moment(&um,Uorb(o,i),1.0);
    accum_moment(&vm,Vorb(o,i),1.0);
    jz = Xorb(o,i)*Vorb(o,i) - Yorb(o,i)*Uorb(o,i);
    accum_moment(&jm,jz,1.0);

    pos[0] = Xorb(o,i);
    pos[1] = Yorb(o,i);
    pos[2] = Zorb(o,i);
    time = Torb(o,i);
    if (Qpot)
      (*pot)(&ndim,pos,acc,&etot,&time);
    else
      etot = Porb(o,i);
    etot += 0.5*(sqr(Uorb(o,i)) + sqr(Vorb(o,i)) + sqr(Worb(o,i)));
    if (Qrel) etot -= e;
    accum_moment(&em,etot,1.0);
  }
  if (first) {
    if (Qell)
      dprintf(0,"# T\tE\tr2_min\tr2_max\tv2_min\tv2_max\tr2_mean\tr2_sigma\tv2_mean\tv2_sigma\tLmin/Lmax\n");
    else
      dprintf(0,"# T\tE\tx_max\ty_max\tu_max\tv_max\tj_mean\tj_sigma\te_mean\te_sigma %s\n",
	      Qrel ? "[rel]" : "[abs]");
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
    printf("%g %g %g %g %g %g %g %g %g %g\n",
	   t,e,
	   max_moment(&xm), max_moment(&ym),
	   max_moment(&um), max_moment(&vm),
	   mean_moment(&jm), sigma_moment(&jm),
	   mean_moment(&em), sigma_moment(&em));
  }
}
