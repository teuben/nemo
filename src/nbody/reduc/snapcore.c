/*
 * SNAPCORE   - sample code for feedback into potname's plummerv
 *
 *   cloned off snapmradii - hence the sometimes odd looking code
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

string defv[] = {
    "in=???\n                   Input file name (snapshot)",
    "potfile=???\n              Output plummerv file",
    "fraction=0.1\n             Fractional masses to print radii of",
    "VERSION=0.1\n              18-jun-03 PJT",
    NULL,
};

string usage="Example";


#define MFRACT 256

local void snapsort(Body *, int , real);
local int rank_aux(Body *, Body *);
local void write_tmr(string potfile, real  t1, real m1, real r1);


void nemo_main()
{
  stream instr;
  string outfile = getparam("potfile");
  real   tsnap, mf[MFRACT], tmass, cmass, fmass, mold, rold;
  real   msnap, rsnap;
  int    k, nbody, bits, nfract;
  Body *btab = NULL, *bp;
  
  nfract = 1;
  mf[0] = getdparam("fraction");

  instr = stropen(getparam("in"), "r");           /* open input file */
  get_history(instr);			    /* accumulate data history */
  for(;;) {				 	 /* loop for all times */
        get_history(instr);                         /* for paranoidici */
        if (!get_tag_ok(instr, SnapShotTag))          /* check if done */
	    break;
        get_snap(instr, &btab, &nbody, &tsnap, &bits);      /* get one */
        if ((bits & PhaseSpaceBit) == 0)
            continue;                       /* if no positions -  skip */ 
        snapsort(btab,nbody,tsnap);        /* sort; hide radius in aux */
        for (bp=btab, tmass=0.0; bp<btab+nbody; bp++)
            tmass += Mass(bp);
        if (tmass == 0.0) error("No masses available in this snapshot");
        k=0;  
        fmass=mf[0]*tmass;
        mold=0.0;
        rold=0.0;
        cmass=0.0;
        for (bp=btab; bp<btab+nbody; bp++) {    /* loop over bodies */
            cmass += Mass(bp);
            if(cmass>=fmass) {
	      msnap = cmass/mf[0];
	      rsnap = absv(Pos(bp));
	      dprintf(0,"snapcore: %g %g %g\n",tsnap,msnap,rsnap);
	      write_tmr(outfile,tsnap,msnap,rsnap);
	      break;
	    }
        }
#if 0
	/* yuck, if you free, the new one may not have masses ... */
        free(btab);
	btab = NULL; /* free old to make sure newly allocated .... (sic)  */
#endif
    }   /* for(;;) */
} /* nemo_main() */

void snapsort(Body *btab, int nbody, real tsnap)
{
    int i;
    Body *b;

    for (i = 0, b = btab; i < nbody; i++, b++)
	Aux(b) = absv(Pos(b));
    qsort(btab, nbody, sizeof(Body), rank_aux);
}

int rank_aux(Body *a, Body *b)
{
    return (Aux(a) < Aux(b) ? -1 : Aux(a) > Aux(b) ? 1 : 0);
}

void write_tmr(string potfile, real  t1, real m1, real r1)
{
  stream f = stropen(potfile,"w!");
  put_set(f,"plummerv");
  put_data(f,"Time",RealType,&t1,0);
  put_data(f,"Mass",RealType,&m1,0);
  put_data(f,"Radius",RealType,&r1,0);
  put_tes(f,"plummerv");
  strclose(f);
}
