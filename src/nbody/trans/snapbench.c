/*
 *  SNAPBENCH:  snapshot benchmark to scale masses (see also the broken testio.c)
 *              also benchmarks vectors of bodies (get_snap) vs. body of vectors (get_snapshot)
 *
 *  
 *     13-mar-2005  Created after Walter's comment at Vegas05         PJT
 *                   (what comment?)
 *     25-feb-2024  ansi cleanup + documented odd behavior            PJT
 *     27-feb-2024  overhaul, just using mode=0,1,2,3                 PJT
 *
 * mode=0     snapshot I/O but keeping linear arrays
 *      1     body with a bodytrans function
 *      2     body with a Mass(bp)
 *      3     body copy to mass[], and scaling that
 *
 *  0 and 3 are same speed, super fast;   0 saves on I/O
 *  1 is slowest 
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <timers.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>
#include <bodytrans.h>

#include <snapshot/get_snapshot.c>

string defv[] = {
    "in=???\n		      Input (snapshot) file",
    "mass=1\n	              Expression for new masses",
    "iter=10\n                Number of iterations to test",
    "mode=1\n                 reading input mode (1=body 2=arrays)",
    "out=\n                   output (snapshot) file, if needed",
    "VERSION=2.0a\n           10-mar-2024 PJT",
    NULL,
};

string usage="benchmark (re)assigning masses to a snapshot";

// using an inline made no difference
inline real mult(real a, real b) { return a*b; }

void ini_snapshot(SS *ssp) 
{
    ssp->nbody = 0;
    ssp->bits = 0;
    ssp->time = 0.0;
    ssp->mass = NULL;
    ssp->phase = NULL;
}



void nemo_main()
{
    stream instr, outstr;
    real   tsnap, mscale, *mass = NULL;
    Body  *btab = NULL, *bp;
    int i, j, nbody, bits;
    rproc  bfunc;
    int mode = getiparam("mode");
    int niter = getiparam("iter");
    
    dprintf(0,"mode=%d\n",mode);
    
    instr = stropen(getparam("in"), "r");
    outstr = hasvalue("out") ? stropen(getparam("out"),"w") : NULL;
    get_history(instr);

    init_timers2(niter,1);
    stamp_timers(0);
    if (mode==0) {
      SS ss;
      ini_snapshot(&ss);
      if (get_snapshot(instr, &ss) == 0) return;
      dprintf(2,"initial mass: %g\n",ss.mass[0]);
      mscale = getdparam("mass");           // use it as adding parameter
      nbody = ss.nbody;
      stamp_timers(1);	
      for (j=0; j<niter; j++)
	for (i=0; i<nbody; i++)
	  ss.mass[i] += mscale;
      dprintf(2,"final mass: %g\n",ss.mass[0]);
    } else {
      get_snap(instr, &btab, &nbody, &tsnap, &bits);
      if (mode == 1) {
	dprintf(0,"bodytrans scaling, iter=%d\n",niter);
	bfunc = btrtrans(getparam("mass"));     /* use bodytrans expression */
	stamp_timers(1);
	for (j=0; j<niter; j++)
	  for (bp=btab, i=0; i<nbody; bp++,i++)
	    Mass(bp) = bfunc(bp, tsnap, i);
	dprintf(2,"final mass:  %g\n", Mass(btab));
      } else {
	mscale = getdparam("mass");           // use it as adding parameter
	if (mode == 2) {
	  dprintf(0,"simple inline scaling, iter=%d mscale=%g\n",niter,mscale);
	  stamp_timers(1);	
	  for (j=0; j<niter; j++)
	    for (bp=btab, i=0; i<nbody; bp++,i++)
	      Mass(bp) += mscale;       // or use mult(mscale, Mass(bp));
	  dprintf(2,"final mass:  %g\n", Mass(btab));
	} else if (mode == 3) {
	  dprintf(0,"simple array scaling, iter=%d mscale=%g\n",niter,mscale);
	  mass = (real *) allocate(nbody * sizeof(real));
	  for (bp=btab, i=0; i<nbody; bp++,i++)
	    mass[i] = Mass(bp);
	  stamp_timers(1);	
	  for (j=0; j<niter; j++)
	    for (i=0; i<nbody; i++)
	      mass[i] += mscale;
	  dprintf(2,"final mass: %g\n",mass[0]);
	} else {
	  error("mode=%d not implemented",mode);
	} /* 2,3, non */
      } /* 1 vs. non-1 */
    } /* 0 vs. non-0 */
    stamp_timers(2);

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
