/*
 * SNAPRSTAT.C: statistics 
 *
 *	22-oct-90  V1.0  replace some of snapstat's functionality	PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/barebody.h>  /* not a full body needed */
#include <snapshot/get_snap.c>

string defv[] = {
    "in=???\n        Input file name (snapshot)",
    "VERSION=1.1\n   21-sep-2011 PJT",
    NULL,
};

string usage="snapshot statistics";



nemo_main()
{
    stream instr;
    real   tsnap, fact;
    vector pmin, pmax, vmin, vmax;
    vector tmp, sump1, sump2, sumv1, sumv2, meanp, meanv, sigp, sigv;
    int    i, nbody, bits;
    Body *btab = NULL, *bp;

    instr = stropen(getparam("in"), "r");           /* open input file */
    get_history(instr);			    /* accumulate data history */
    for(;;) {				 	 /* loop for all times */
        get_history(instr);                         /* for paranoidici */
        if (!get_tag_ok(instr, SnapShotTag))          /* check if done */
	    break;
        get_snap(instr, &btab, &nbody, &tsnap, &bits);      /* get one */
        if ((bits & PhaseSpaceBit) == 0)
            continue;                       /* if no positions -  skip */
        CLRV(sump1);				 /* clear some vectors */
        CLRV(sump2);
        CLRV(sumv1);
        CLRV(sumv2);
	SETV(pmin,Pos(btab));
	SETV(pmax,Pos(btab));
	SETV(vmin,Vel(btab));
	SETV(vmax,Vel(btab));
        for (bp = btab; bp < btab+nbody; bp++) {    /* loop all bodies */
            SADDV(sump1,Pos(bp));
            MULVV(tmp,Pos(bp),Pos(bp));
            SADDV(sump2,tmp);
            SADDV(sumv1,Vel(bp));
            MULVV(tmp,Vel(bp),Vel(bp));
            SADDV(sumv2,tmp);
	    SMINV(pmin,Pos(bp))
	    SMAXV(pmax,Pos(bp))
	    SMINV(vmin,Vel(bp))
	    SMAXV(vmax,Vel(bp))
        }   /* for (bp) */
        fact = 1.0 / nbody;
        SMULVS(sump1,fact);                     /* position statistics */
        SETV(meanp,sump1);
        MULVV(sump1,sump1,sump1);
        SMULVS(sump2,fact);
        SUBV(sigp,sump2,sump1);
        SQRTV(sigp);

        SMULVS(sumv1,fact);                     /* velocity statistics */
        SETV(meanv,sumv1);
	MULVV(sumv1,sumv1,sumv1);
        SMULVS(sumv2,fact);
        SUBV(sigv,sumv2,sumv1);
        SQRTV(sigv);
                                                  /* Output of results */
        printf("time: %f\n",tsnap);

        printf("pos: ");
        for (i=0; i<NDIM; i++)
            printf(" %lf +/- %lf ",meanp[i],sigp[i]);
        printf("\n");

        printf("vel: ");
        for (i=0; i<NDIM; i++)
            printf(" %lf +/- %lf ",meanv[i],sigv[i]);
        printf("\n");

	printf("mnmx pos: ");
        for (i=0; i<NDIM; i++)
            printf(" %lf %lf ",pmin[i],pmax[i]);
        printf("\n");

	printf("mnmx vel: ");
        for (i=0; i<NDIM; i++)
            printf(" %lf %lf ",vmin[i],vmax[i]);
        printf("\n");

	btab = NULL; /* free old to make sure newly allocated .... (sic)  */
    }   /* for(;;) */
} /* nemo_main() */
