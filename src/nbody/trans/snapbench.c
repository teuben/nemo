/*
 *  SNAPBENCH:  snapshot benchmark to scale masses
 *
 *    mkplummer p6 10000000 massname='n(m)' massrange=1,2
 *    /usr/bin/time snapbench p6 'mass=3.1415'   bodytrans=f iter=1000    6.6"
 *    /usr/bin/time snapbench p6 'mass=3.1415'   bodytrans=t iter=1000   13.6
 *    /usr/bin/time snapbench p6 'mass=3.1415*m' bodytrans=t iter=1000    6.8
 *  
 *     13-mar-05  Created after Walter's comment at Vegas05            PJT
 *     25-feb-2024   ansi cleanup + documented odd behavior            PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <timers.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>
#include <bodytrans.h>

string defv[] = {
    "in=???\n		      Input (snapshot) file",
    "mass=1\n	              Expression for new masses",
    "bodytrans=t\n            Use bodytrans",
    "iter=10\n                Number of iterations to test",
    "body=t\n                 Keep Body() structure, or simple array?",
    "out=\n                   output (snapshot) file, if needed",
    "VERSION=1.2\n            27-feb-2024 PJT",
    NULL,
};

string usage="benchmark (re)assign masses to a snapshot";

// using an inline made no difference
inline real mult(real a, real b) { return a*b; }


void nemo_main()
{
    stream instr, outstr;
    real   tsnap, mscale, *mass;
    Body  *btab = NULL, *bp;
    int i, j, nbody, bits;
    rproc  bfunc;
    bool Qtrans = getbparam("bodytrans");
    bool Qbody = getbparam("body");
    int niter = getiparam("iter");

    if (!Qbody) Qtrans=FALSE;

    init_timers2(niter,1);
    stamp_timers(0);

    instr = stropen(getparam("in"), "r");
    outstr = hasvalue("out") ? stropen(getparam("out"),"w") : NULL;
    
    get_history(instr);
    if (!get_tag_ok(instr, SnapShotTag)) 
      error("not a snapshot");
    get_snap(instr, &btab, &nbody, &tsnap, &bits);
    
    if (Qtrans) {
      dprintf(0,"bodytrans scaling, iter=%d\n",niter);
      bfunc = btrtrans(getparam("mass"));     /* use bodytrans expression */
      stamp_timers(1);
      for (j=0; j<niter; j++)
	for (bp=btab, i=0; i<nbody; bp++,i++)
	  Mass(bp) = bfunc(bp, tsnap, i);
      dprintf(2,"final mass:  %g\n", Mass(btab));
    } else {
      mscale = getdparam("mass");           // use it as adding parameter
      if (Qbody) {
	dprintf(0,"simple inline scaling, iter=%d mscale=%g\n",niter,mscale);
	stamp_timers(1);	
	for (j=0; j<niter; j++)
	  for (bp=btab, i=0; i<nbody; bp++,i++)
	    Mass(bp) += mscale;       // or use mult(mscale, Mass(bp));
	dprintf(2,"final mass:  %g\n", Mass(btab));
      } else {
	dprintf(0,"simple array scaling, iter=%d mscale=%g\n",niter,mscale);
	mass = (real *) allocate(nbody * sizeof(real));
	for (bp=btab, i=0; i<nbody; bp++,i++)
	  mass[i] = Mass(bp);
	stamp_timers(1);	
	for (j=0; j<niter; j++)
	  for (i=0; i<nbody; i++)
	    mass[i] += mscale;
	dprintf(2,"final mass: %g\n",mass[0]);
      }
    }
    stamp_timers(2);

    strclose(instr);
    if (outstr) { 
      put_history(outstr);
      put_snap(outstr, &btab, &nbody, &tsnap, &bits);
      strclose(outstr);
    }
    stamp_timers(3);
    dprintf(0,"%Ld %Ld %Ld ticks\n",
	    diff_timers(0,1),
	    diff_timers(1,2),
	    diff_timers(2,3));
	    
    dprintf(0,"%g %g %g sec\n",
	    diff_timers2(0,1),
	    diff_timers2(1,2),
	    diff_timers2(2,3));
	    
}
