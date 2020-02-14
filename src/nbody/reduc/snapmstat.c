/*
 * SNAPMSTAT.C: show some statistics of the masses in a snapshot
 *              Only looks at first snapshot
 *
 *      If sorting is done, the histogram is proper, and masses
 *      could in principle be tagges in 'Key' if <body.h> is
 *      used.
 *
 *	18-sep-90  V1.0  created - no sorting done yet		PJT
 *	12-may-92  V1.1  default sorted by mass
 *	15-mar-95  V1.2  default not sorted by mass, output format	PJT
 *       6-feb-2010  V2.0  also report a select= type for glnemo2       PJT
 *                         default for sort=f
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>	
#include <snapshot/barebody.h>  /* not a full body needed */
#include <snapshot/get_snap.c>

string defv[] = {
    "in=???\n       Input file name",
    "sort=f\n       Sort masses before processing?",    
    "species=100\n  Maximum number of species",
    "VERSION=2.0\n  6-feb-2020 PJT",
    NULL,
};

string usage="report some statistics of the masses in a snapshot";


local void snapsort(Body *, int);
local int rank_mass(Body *, Body *);

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
    if ((bits & MassBit) == 0)
	error("No masses");
    if (Qsort)
        snapsort(btab,nbody);

    mold = Mass(btab);      /* set first mass */
    mtot = mold;
    mcum = mold;
    iold = 0;
    scount = -1;
    for (i=0, bp = btab; bp < btab+nbody; i++, bp++) {
        if (Mass(bp) != mold) {
           printf("%d %d:%d  = %d Mass= %g TotMas= %g CumMas= %g\n",
                scount+1, iold, i-1, i-iold,mold, mtot, mcum);
	   nsp[scount+1] = i-iold;
           scount++;
	   if (scount == maxspecies) error("Too many mass species, use sort=t ?");
           mold = Mass(bp);    /* remember old mass */
           iold = i;           /* remember where we were */
           mtot = 0.0;         /* reset total mass */
        }
        if (i>0) {
            mtot += Mass(bp);
            mcum += Mass(bp);
        }
    }
    printf("%d %d:%d = %d Mass= %g TotMas= %g CumMas= %g\n",
                scount+1, iold, i-1, i-iold, mold, mtot, mcum);
    nsp[scount+1] = i-iold;
    scount++;
    if (scount == maxspecies) error("Too many mass species, use sort=t ?");
    scount++;

    if (!Qsort) {
      printf("# Found %d species:\n",scount);
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
    qsort(btab, nbody, sizeof(Body), rank_mass);
}

local int rank_mass(Body *a, Body *b)
{
    return (Mass(a) < Mass(b) ? -1 : Mass(a) > Mass(b) ? 1 : 0);
}

