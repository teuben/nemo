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
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/barebody.h>  /* no a full body needed */
#include <snapshot/get_snap.c>

string defv[] = {
    "in=???\n       Input file name",
    "sort=t\n       Sort masses before processing?",    
    "species=100\n  Maximum number of species",
    "VERSION=1.2\n  15-mar-95 PJT",
    NULL,
};

string usage="report some statistics of the masses in a snapshot";


local void snapsort(Body *, int);
local int rank_mass(Body *, Body *);

nemo_main()
{
    stream instr;
    real   mold, tsnap, mtot, mcum;
    int    i, iold, nbody, bits, nspecies, maxspecies, scount;
    bool   msort;
    Body *btab = NULL, *bp;

    msort = getbparam("sort");
    maxspecies = getiparam("species");

    instr = stropen(getparam("in"), "r");
    get_history(instr);
    if (!get_tag_ok(instr, SnapShotTag))
	error("Not a snapshot");
    get_snap(instr, &btab, &nbody, &tsnap, &bits);
    if ((bits & MassBit) == 0)
	error("No masses");
    if (msort)
        snapsort(btab,nbody);

    mold = Mass(btab);      /* set first mass */
    mtot = mold;
    mcum = mold;
    iold = 0;
    scount = -1;
    for (i=0, bp = btab; bp < btab+nbody; i++, bp++) {  /* first pass */
        if (Mass(bp) != mold) {
           printf("%d %d:%d  = %d Mass= %g TotMas= %g CumMas= %g\n",
                scount+1, iold, i-1, i-iold,mold, mtot, mcum);
           scount++;
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
}


local void snapsort(Body *btab, int nbody)
{
    qsort(btab, nbody, sizeof(Body), rank_mass);
}

local int rank_mass(Body *a, Body *b)
{
    return (Mass(a) < Mass(b) ? -1 : Mass(a) > Mass(b) ? 1 : 0);
}

