/*
 * SNAPSYM.C:    symmeterize a snapshot
 *
 *     23-oct-2021    V0.1   Created       - Peter Teuben
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
    "mode=\n        (TBD) mode of symmetry",
    "VERSION=0.1\n  23-oct-2021 PJT",
    NULL,
};

string usage="symmetrize a snapshot";

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
  real   x,y,z,u,v,w,tsnap;
  string times;
  int i, nbody, bits, nbody8, bits8;
  Body *btab = NULL, *bpi, *bpo;
  Body *btab8 = NULL;
  
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
    dprintf (1,"Symmetrizing snapshot at time= %f bits=0x%x\n",tsnap,bits);
    bits8  = (MassBit | PhaseSpaceBit);
    nbody8 = 8*nbody;
    btab8  = (Body *) allocate(nbody8*sizeof(Body));


    bpo = btab8;
    for (bpi = btab; bpi < btab+nbody; bpi++, bpo++) {              // Q1 +z
      x = Pos(bpi)[0];    y = Pos(bpi)[1];   z = Pos(bpi)[2];
      u = Vel(bpi)[0];    v = Vel(bpi)[1];   w = Vel(bpi)[2];
      Pos(bpo)[0] = +x;   Pos(bpo)[1] = +y;  Pos(bpo)[2] = +z;
      Vel(bpo)[0] = +u;   Vel(bpo)[1] = +v;  Vel(bpo)[2] = +w;
      Mass(bpo) = Mass(bpi);
    }

    for (bpi = btab; bpi < btab+nbody; bpi++, bpo++) {              // Q2 +z
      x = Pos(bpi)[0];    y = Pos(bpi)[1];   z = Pos(bpi)[2];
      u = Vel(bpi)[0];    v = Vel(bpi)[1];   w = Vel(bpi)[2];
      Pos(bpo)[0] = -x;   Pos(bpo)[1] = +y;  Pos(bpo)[2] = +z;
      Vel(bpo)[0] = -v;   Vel(bpo)[1] = +u;  Vel(bpo)[2] = +w;
      Mass(bpo) = Mass(bpi);
    }

    for (bpi = btab; bpi < btab+nbody; bpi++, bpo++) {              // Q3 +z
      x = Pos(bpi)[0];    y = Pos(bpi)[1];   z = Pos(bpi)[2];
      u = Vel(bpi)[0];    v = Vel(bpi)[1];   w = Vel(bpi)[2];
      Pos(bpo)[0] = -x;   Pos(bpo)[1] = -y;  Pos(bpo)[2] = +z;
      Vel(bpo)[0] = -u;   Vel(bpo)[1] = -v;  Vel(bpo)[2] = +w;
      Mass(bpo) = Mass(bpi);
    }

    for (bpi = btab; bpi < btab+nbody; bpi++, bpo++) {              // Q4 +z
      x = Pos(bpi)[0];    y = Pos(bpi)[1];   z = Pos(bpi)[2];
      u = Vel(bpi)[0];    v = Vel(bpi)[1];   w = Vel(bpi)[2];
      Pos(bpo)[0] = +x;   Pos(bpo)[1] = -y;  Pos(bpo)[2] = +z;
      Vel(bpo)[0] = +v;   Vel(bpo)[1] = -u;  Vel(bpo)[2] = +w;
      Mass(bpo) = Mass(bpi);
    }
    
    for (bpi = btab; bpi < btab+nbody; bpi++, bpo++) {              // Q1 -z
      x = Pos(bpi)[0];    y = Pos(bpi)[1];   z = Pos(bpi)[2];
      u = Vel(bpi)[0];    v = Vel(bpi)[1];   w = Vel(bpi)[2];
      Pos(bpo)[0] = +x;   Pos(bpo)[1] = +y;  Pos(bpo)[2] = -z;
      Vel(bpo)[0] = +u;   Vel(bpo)[1] = +v;  Vel(bpo)[2] = -w;
      Mass(bpo) = Mass(bpi);
    }

    for (bpi = btab; bpi < btab+nbody; bpi++, bpo++) {              // Q2 -z
      x = Pos(bpi)[0];    y = Pos(bpi)[1];   z = Pos(bpi)[2];
      u = Vel(bpi)[0];    v = Vel(bpi)[1];   w = Vel(bpi)[2];
      Pos(bpo)[0] = -x;   Pos(bpo)[1] = +y;  Pos(bpo)[2] = -z;
      Vel(bpo)[0] = -v;   Vel(bpo)[1] = +u;  Vel(bpo)[2] = -w;
      Mass(bpo) = Mass(bpi);
    }

    for (bpi = btab; bpi < btab+nbody; bpi++, bpo++) {              // Q3 -z
      x = Pos(bpi)[0];    y = Pos(bpi)[1];   z = Pos(bpi)[2];
      u = Vel(bpi)[0];    v = Vel(bpi)[1];   w = Vel(bpi)[2];
      Pos(bpo)[0] = -x;   Pos(bpo)[1] = -y;  Pos(bpo)[2] = -z;
      Vel(bpo)[0] = -u;   Vel(bpo)[1] = -v;  Vel(bpo)[2] = -w;
      Mass(bpo) = Mass(bpi);
    }

    for (bpi = btab; bpi < btab+nbody; bpi++, bpo++) {              // Q4 -z
      x = Pos(bpi)[0];    y = Pos(bpi)[1];   z = Pos(bpi)[2];
      u = Vel(bpi)[0];    v = Vel(bpi)[1];   w = Vel(bpi)[2];
      Pos(bpo)[0] = +x;   Pos(bpo)[1] = -y;  Pos(bpo)[2] = -z;
      Vel(bpo)[0] = +v;   Vel(bpo)[1] = -u;  Vel(bpo)[2] = -w;
      Mass(bpo) = Mass(bpi);
    }
    
    put_snap(outstr, &btab8, &nbody8, &tsnap, &bits8);
  }
}

