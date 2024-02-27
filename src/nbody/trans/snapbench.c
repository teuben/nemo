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
    "in=???\n		      input (snapshot) file",
    "mass=1\n		      expression for new masses",
    "bodytrans=t\n            Use bodytrans",
    "iter=10\n                Number of iterations to test",
    "out=\n                   output (snapshot) file, if needed",
    "VERSION=1.1\n            25-feb-2024 PJT",
    NULL,
};

string usage="benchmark (re)assign masses to a snapshot";

// using an inline made no difference
inline real mult(real a, real b) { return a*b; }


void nemo_main()
{
    stream instr, outstr;
    real   tsnap, mscale;
    Body  *btab = NULL, *bp;
    int i, j, nbody, bits;
    rproc  bfunc;
    real t0,t1,t2;
    bool Qtrans = getbparam("bodytrans");
    int niter = getiparam("iter");

    init_timers(niter);
    stamp_timers(0);

    instr = stropen(getparam("in"), "r");
    outstr = hasvalue("out") ? stropen(getparam("out"),"w") : NULL;
    

    get_history(instr);
    if (!get_tag_ok(instr, SnapShotTag)) 
      error("not a snapshot");
    get_snap(instr, &btab, &nbody, &tsnap, &bits);
    stamp_timers(1);
    
    if (Qtrans) {
      dprintf(0,"bodytrans scaling, iter=%d\n",niter);
      bfunc = btrtrans(getparam("mass"));     /* use bodytrans expression */

      for (j=0; j<niter; j++)
	for (bp=btab, i=0; i<nbody; bp++,i++)
	  Mass(bp) = bfunc(bp, tsnap, i);
    } else {
      dprintf(0,"simple inline scaling, iter=%d\n",niter);
      mscale = getdparam("mass");            // use it as scaling parameter
      for (j=0; j<niter; j++)
	for (bp=btab, i=0; i<nbody; bp++,i++)
	  Mass(bp) = mscale*Mass(bp);    // or use mult(mscale, Mass(bp));
    }
    stamp_timers(2);

    strclose(instr);
    if (outstr) { 
      put_history(outstr);
      put_snap(outstr, &btab, &nbody, &tsnap, &bits);
      strclose(outstr);
    }
    stamp_timers(3);
    dprintf(0,"%Ld %Ld %Ld\n",
	    diff_timers(0,1),
	    diff_timers(1,2),
	    diff_timers(2,3));
	    
}
