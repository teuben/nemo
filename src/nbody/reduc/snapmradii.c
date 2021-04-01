/*
 * SNAPMRADII.C: mass radii
 *
 *	12-may-92  V1.0  created	PJT
 *      19-nov-93  V1.1  don't free(btab) all the time;    PJT
 *	10-jun-95     a  qsort from <stdlib.h> (linux)     pjt
 *	10-aug-95  V1.2  added tab= keyword for full table output	pjt
 *      15-aug-96  V1.3  fixed bug in which total mass=1 was assumed    pjt
 *	25-mar-97     a  fixed for SINGLEPREC				pjt
 *      10-mar-04  V1.4  add log=                                       pjt
 *      27-jul-05   1.5  added sort=                                    pjt
 *       1-apr-21   1.6  deal with no masses in snapshot for Tjeerd     pjt
 *
 *  Bug: if the massfractions are too close such that there
 *       are bins withouth mass, this algorithm fails
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <bodytransc.h>

string defv[] = {
    "in=???\n                   Input file name (snapshot)",
    "fraction=0.1:0.9:0.1\n     Fractional masses to print radii of",
    "tab=f\n			Full table of r,m(r) ? ",
    "log=f\n                    Print radii in log10() ? ",
    "sort=r\n                   Observerble to sort masses by",
    "VERSION=1.6\n              1-apr-2021 PJT",
    NULL,
};

string usage="Langrangian mass-fraction radii of a snapshot";

string cvsid="$Id$";


#define MFRACT 256

local void snapsort(Body *, int , real, rproc_body);
local int rank_aux(Body *, Body *);


void nemo_main()
{
    stream instr;
    real   tsnap, mf[MFRACT], tmass, cmass, fmass, mold, rold, rlag;
    int    k, nbody, bits, nfract;
    bool   Qtab = getbparam("tab");
    bool   Qlog = getbparam("log");
    Body *btab = NULL, *bp;
    rproc_body sortptr;

    sortptr = btrtrans(getparam("sort"));
    
    nfract = nemoinpr(getparam("fraction"),mf,MFRACT);
    if (nfract<1) error("Illegal or bad parsed fraction=%s",
                        getparam("fraction"));
    for (k=0; k<nfract; k++) {    /* check for being sorted */
        if (mf[k]<0.0 || mf[k]>1.0) error("Bad fraction[%d]=%g (must be (0,1)",
                                           k+1,mf[k]);
        if (k>0 && mf[k-1]>=mf[k]) error("fraction= must be sorted");
    }

    instr = stropen(getparam("in"), "r");           /* open input file */
    get_history(instr);			    /* accumulate data history */
    for(;;) {				 	 /* loop for all times */
        get_history(instr);                         /* for paranoidici */
        if (!get_tag_ok(instr, SnapShotTag))          /* check if done */
	    break;
        get_snap(instr, &btab, &nbody, &tsnap, &bits);      /* get one */
        if ((bits & PhaseSpaceBit) == 0)
            continue;                       /* if no positions -  skip */
        if (!Qtab) printf("%g",tsnap);
        snapsort(btab,nbody,tsnap,sortptr);/* sort; hides radius in aux */
        for (bp=btab, tmass=0.0; bp<btab+nbody; bp++)
            tmass += Mass(bp);
        if (tmass == 0.0) {
	  warning("No masses available in this snapshot- using equal masses");
	  tmass = 1.0;
	  for (bp=btab;  bp<btab+nbody; bp++)
	    Mass(bp) = 1.0/nbody;
	}
        k=0;  
        fmass=mf[0]*tmass;
        mold=0.0;
        rold=0.0;
        cmass=0.0;
        for (bp=btab; bp<btab+nbody; bp++) {    /* loop over bodies */
            cmass += Mass(bp);
            if(cmass>=fmass) {
            	if (Qtab) printf("%g", mf[k]);
		rlag = rold + (tmass*mf[k]-mold)*(Aux(bp)-rold)/(cmass-mold);
		if (Qlog) rlag = log10(rlag);
                printf(" %g", rlag);
                if (Qtab) printf("\n");
                k++;
                if (k >= nfract) break;
                fmass = mf[k] * tmass;
            }
            rold = Aux(bp);
            mold = cmass;
        }
        if (!Qtab) printf("\n");
#if 0
	/* yuck, if you free, the new one may not have masses ... */
        free(btab);
	btab = NULL; /* free old to make sure newly allocated .... (sic)  */
#endif
    }   /* for(;;) */
} /* nemo_main() */

void snapsort(Body *btab, int nbody, real tsnap, rproc_body sortptr) 
{
    int i;
    Body *b;

    for (i = 0, b = btab; i < nbody; i++, b++)
      Aux(b) = sortptr(b,tsnap,i);
    qsort(btab, nbody, sizeof(Body), rank_aux);
}

int rank_aux(Body *a, Body *b)
{
    return (Aux(a) < Aux(b) ? -1 : Aux(a) > Aux(b) ? 1 : 0);
}

