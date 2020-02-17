/*
 * SNAPKSTAT.C: show some statistics of the keys in a snapshot
 *              Only looks at first snapshot
 *
 *      If sorting is done, the histogram is proper, and masses
 *      could in principle be tagges in 'Key' if <body.h> is
 *      used.
 *
 *      17-feb-2010  V0.1   quick hack, from snapmstate     PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>  
#include <snapshot/get_snap.c>

string defv[] = {
    "in=???\n       Input file name",
    "sort=f\n       Sort masses before processing?",    
    "species=100\n  Maximum number of species",
    "VERSION=0.2\n  17-feb-2020 PJT",
    NULL,
};

string usage="report some statistics of the keys in a snapshot";


local void snapsort(Body *, int);
local int rank_key(Body *, Body *);

void nemo_main(void)
{
    stream instr;
    real   mold, tsnap, mtot, mcum;
    int    i, iold, icum, nbody, bits, nspecies, maxspecies, scount;
    bool   Qsort;
    Body *btab = NULL, *bp;
    int *nsp;

    Qsort = getbparam("sort");
    maxspecies = getiparam("species");
    nsp = (int *) allocate(maxspecies * sizeof(int));

    instr = stropen(getparam("in"), "r");
    get_history(instr);
    if (!get_tag_ok(instr, SnapShotTag))
	error("Not a snapshot");
    get_snap(instr, &btab, &nbody, &tsnap, &bits);
    if ((bits & KeyBit) == 0)
	error("No keys");
    if (Qsort)
        snapsort(btab,nbody);

    mold = Key(btab);      /* set first key */
    mtot = mold;
    mcum = mold;
    iold = 0;
    scount = -1;
    for (i=0, bp = btab; bp < btab+nbody; i++, bp++) {
        if (Key(bp) != mold) {
	   dprintf(1,"%d %d:%d  = %d Key= %g TotKey= %g CumKey= %g\n",
                scount+1, iold, i-1, i-iold,mold, mtot, mcum);
	   nsp[scount+1] = i-iold;
           scount++;
	   if (scount == maxspecies) error("Too many key species, use sort=t ?");
           mold = Key(bp);     /* remember old key */
           iold = i;           /* remember where we were */
           mtot = 0.0;         /* reset total mass */
        }
        if (i>0) {
            mtot += Key(bp);
            mcum += Key(bp);
        }
    }
    dprintf(1,"%d %d:%d = %d Key= %g TotKey= %g CumKey= %g\n",
                scount+1, iold, i-1, i-iold, mold, mtot, mcum);
    nsp[scount+1] = i-iold;
    scount++;
    if (scount == maxspecies) error("Too many key species, use sort=t ?");
    scount++;

    if (!Qsort) {
      dprintf(0,"# Found %d key species:\n",scount);
      printf("select=");
      
      iold = 0;
      icum = 0;
      for (i=0; i<scount; i++) {
	icum = icum + nsp[i];
	if (i==scount-1)
	  printf("%d:%d ",iold,icum-1);
	else
	  printf("%d:%d,",iold,icum-1);
	iold = icum ;
      }
      printf("\n");
    }
}


local void snapsort(Body *btab, int nbody)
{
  qsort(btab, nbody, sizeof(Body), rank_key);
}

local int rank_key(Body *a, Body *b)
{
  return (Key(a) < Key(b) ? -1 : Key(a) > Key(b) ? 1 : 0);
}

