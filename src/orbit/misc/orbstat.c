/*
 * ORBSTAT:	 tabulate orbit and some of their statistics
 *
 *	23-mar-95 V1	Created					pjt
 *	14-apr-01       header now in dprintf()
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <orbit.h>
#include <moment.h>

string defv[] = {
    "in=???\n		Input orbit",
    "VERSION=1.0a\n  	14-apr-01 PJT",
    NULL,
};

string usage="Tabulate some statistics of orbit(s)";


string	infile;			/* file name */
stream  instr;			/* file stream */

orbitptr optr;

nemo_main ()
{
    string mode;
    bool pabs,vabs;

    infile = getparam("in");
    instr = stropen (infile,"r");

    optr=NULL;
    while (read_orbit (instr,&optr)) 
        stat_orbit(optr);
    strclose(instr);
}



stat_orbit(orbitptr o)
{
    Moment  xm, ym, um, vm, jm;
    real jz;
    real t, e;
    int i;
    permanent bool first = TRUE;
    
    ini_moment(&xm,-1);
    ini_moment(&ym,-1);
    ini_moment(&um,-1);
    ini_moment(&vm,-1);
    ini_moment(&jm,2);
    
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
    	dprintf(0,"# T\tE\tx_max\ty_max\tu_max\tv_max\tj_mean\tj_sigma\n");
    	first=FALSE;
    }
    printf("%g %g %g %g %g %g %g %g\n",t,e,
    	max_moment(&xm), max_moment(&ym),
    	max_moment(&um), max_moment(&vm),
    	mean_moment(&jm), sigma_moment(&jm));
}
