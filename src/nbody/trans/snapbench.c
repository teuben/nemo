/*
 *  SNAPBENCH:  snapshot benchmark to scale masses
 *
 *    mkplummer p6 1000000 massname='n(m)' massrange=1,2
 *    time snapbench p6 'mass=3.1415'   bodytrans=f iter=10
 *    time snapbench p6 'mass=3.1415*m' bodytrans=f iter=10  
 *  
 *     13-mar-05  Created after Walter's comment at Vegas05          PJT
 */
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {
    "in=???\n		      input (snapshot) file",
    "mass=1\n		      expression for new masses",
    "bodytrans=t\n            Use bodytrans",
    "iter=10\n                Number of iterations to test",
    "out=\n                   output (snapshot) file, if needed",
    "VERSION=1.0\n            13-mar-05 PJT",
    NULL,
};

string usage="(re)assign masses to a snapshot";

#define TIMEFUZZ	0.0001	/* tolerance in time comparisons */

extern double frandom(double, double, rproc);

nemo_main()
{
    stream instr, outstr;
    real   tsnap, mscale;
    Body  *btab = NULL, *bp;
    int i, j, n, nbody, nbodymass, bits, bitsmass, seed;
    rproc  bfunc, btrtrans();

    instr = stropen(getparam("in"), "r");
    outstr = hasvalue("out") ? stropen(getparam("out"),"w") : NULL;

    n = getiparam("iter");

    get_history(instr);
    if (!get_tag_ok(instr, SnapShotTag)) 
      error("not a snapshot");
    get_snap(instr, &btab, &nbody, &tsnap, &bits);

    if (getbparam("bodytrans")) {
      dprintf(0,"bodytrans scaling, iter=%d\n",n);
      bfunc = btrtrans(getparam("mass"));     /* use bodytrans expression */

      for (j=0; j<n; j++)
	for (bp=btab, i=0; i<nbody; bp++,i++)
	  Mass(bp) = bfunc(bp, tsnap, i);
    } else {
      dprintf(0,"simple inline scaling, iter=%d\n",n);
      mscale = getdparam("mass");
      for (j=0; j<n; j++)
	for (bp=btab, i=0; i<nbody; bp++,i++)
	  Mass(bp) = mscale*Mass(bp);
    }

    strclose(instr);
    if (outstr) { 
      put_history(outstr);
      put_snap(outstr, &btab, &nbody, &tsnap, &bits);
      strclose(outstr);
    }
}
