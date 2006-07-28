/*
 * SNAPSQUASH.C: squash a set of particles with some rotation,
 *               following the recipe in Bosma's (1978) thesis (ch. 5.4)
 *
 *     28-jul-2006    V1.0    Created       - Peter Teuben & Albert Bosma
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
    "in=???\n	    Input file name",
    "out=???\n	    Output file name",
    "f=1.0\n        Factor to squash by",
    "omega=0.0\n    Pattern Speed",
    "times=all\n    Times to select snapshots from",
    "VERSION=1.0\n  28-jul-06 PJT",
    NULL,
};

string usage="squash a snapshot using Bosma's thesis recipe";

string cvsid="$Id$";

#define TIMEFUZZ	0.000001	/* tolerance in time comparisons */

bool uscalar(real x)
{
    return x==1.0;
}

bool uvector(vector v)
{
    register int i;

    for (i=0; i<NDIM; i++)
        if(v[i] != 1.0) return FALSE;
    return TRUE;
}


void nemo_main()
{
  stream instr, outstr;
  real   f, omega, fac1, fac2, fac3, fac4,x,y,u,v,rp,vc,tsnap;
  string times;
  int i, nbody, bits;
  Body *btab = NULL, *bp;
  
  f = getdparam("f");
  omega = getdparam("omega");
  times = getparam("times");
  fac1 = (f*f-1)*omega;
  fac2 = f*f;
  fac3 = (1.0-1.0/fac2)*omega;
  fac4 = 1.0/fac2;

  instr = stropen(getparam("in"), "r");   /* open files */
  outstr = stropen(getparam("out"), "w");

  get_history(instr);
  put_history(outstr);		
  for (;;) {
    get_history(instr);		/* skip over stuff we can forget */
    if (!get_tag_ok(instr, SnapShotTag))
      break;			/* done with work in loop */
    get_snap(instr, &btab, &nbody, &tsnap, &bits);
    if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0) {
      continue;       /* just skip it's probably a diagnostics */
    }
    
    if ((bits & TimeBit) == 0)
      tsnap = 0.0;
    else if (!streq(times,"all") && !within(tsnap, times, TIMEFUZZ))
      continue;
    dprintf (1,"Scaling snapshot at time= %f bits=0x%x\n",tsnap,bits);
    
    for (bp = btab; bp < btab+nbody; bp++) {
      x = Pos(bp)[0];
      y = Pos(bp)[1];
      Pos(bp)[0] = x*f;
      Pos(bp)[1] = y/f;
      u = Vel(bp)[0];
      v = Vel(bp)[1];
      rp = sqrt(x*x+y*y);
      vc = sqrt(u*u+v*v);     /* wrong if we don't have circular orbits */
      Vel(bp)[0] = y*(fac1-fac2*vc/rp);
      Vel(bp)[1] = x*(fac3+fac4*vc/rp);
    }
    put_snap(outstr, &btab, &nbody, &tsnap, &bits);
  }
}

