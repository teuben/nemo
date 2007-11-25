/****************************************************************************/
/* TREETEST.C: routines to make test data for tree-sph code.                */
/* Copyright (c) 2000 by Jin Koda, Tokyo, JAPAN.                            */
/****************************************************************************/

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#define global                                  /* don't default to extern  */
#include "treecode.h"

/* Local routines and variables to set initial condition. */

local void inithknl(bodyptr, int);
local void threadinit(nodeptr, nodeptr, real);

#define Size(x)   (((nodeptr) (x))->mass)
/*
 * TESTDATA: generate initial disk.
 */

#define MFRAC  0.999                            /* cut off 1-MFRAC of mass  */

void testdata(void)
{
    real rsc, r, v, ar;
    vector rcm, vcm;
    bodyptr p;

    bodytab = (bodyptr) allocate(nbody * sizeof(body));
                                                /* alloc space for bodies   */
    rsc = rinit / rscale;                       /* init. radius in sys unit */
    CLRV(rcm);                                  /* zero out cm position     */
    for (p = bodytab; p < bodytab+nbody; p++) { /* loop over bodies         */
        Type(p) = BODY;                         /* tag as a body            */
        Mass(p) = fgas * massdk / nbody;        /* set masses equal         */
	pickball(Pos(p), NDIM, rsc);
	ABSV(r, Pos(p));
        ADDMULVS(rcm, Pos(p), 1.0 / nbody);     /* accumulate cm position   */
	Csnd(p) =  9.291746e-11 * tscale / rscale * rsqrt(1.0e4);
                                                /* sound velocity           */
    }
    for (p = bodytab; p < bodytab+nbody; p++) { /* loop over bodies again   */
        SUBV(Pos(p), Pos(p), rcm);              /* subtract cm position     */
	CLRV(Acc(p));                           /* clear acceleration       */
    }
    if (selfgrav) {                             /* if compute self-gravity  */
	planttree(bodytab, nbody);              /* construct tree structure */
	managetree();                           /* calc com and branch link */
	gravcalc();                             /* and compute grav. potent.*/
    }
    extforce(bodytab, nbody, tnow);             /* compute external forces  */

    CLRV(vcm);                                  /* zero out cm velocity     */
    for (p = bodytab; p < bodytab+nbody; p++) { /* loop over bodies         */
        ABSV(r, Pos(p));                        /* calc radius              */
	DOTVP(ar, Pos(p), Acc(p));              /* dot product of pos & acc */
	v = rsqrt(ABS(ar));                     /* calc tangential velocity,*/
	Vel(p)[0] = v * (-Pos(p)[1]/r);         /* and set velocity vector  */
	Vel(p)[1] = v * ( Pos(p)[0]/r);
        ADDMULVS(vcm, Vel(p), 1.0 / nbody);     /* accumulate cm velocity   */
    }
    for (p = bodytab; p < bodytab+nbody; p++)   /* loop over bodies again   */
        SUBV(Vel(p), Vel(p), vcm);              /* subtract cm velocity     */
    inithknl(bodytab, nbody);                   /* initialize hknl          */
}

/*
 * INITHKNL: initialize the value of hknl.
 */

local void inithknl(bodyptr btab, int nbody)
{
    int actlen;
    nodeptr *active;
    bodyptr *ngbrptr, p;

    planttree(btab, nbody);                     /* construct tree structure */
    threadinit((nodeptr) root, NULL, rsize);    /* thread tree for search   */
    actlen = nbody;                             /* set length of active list*/
    active = (nodeptr *) allocate(actlen * sizeof(nodeptr));
    active[0] = (nodeptr) root;                 /* search from root         */
    ngbrlen = nbody;                            /* set length of ngbr list  */
    ngbrlist = (bodyptr *) allocate(ngbrlen * sizeof(bodyptr));
    for (p = btab; p < btab+nbody; p++) {       /* loop over all body       */
	Hknl(p) = rsqrt((real) nnbr / (real) nbody);
                                                /* initial guess of rad.    */
	do {
	    Hknl(p) *= 1.2;                     /* make hknl larger         */
	    Ngbrs(p) = ngbrlist;
	    Ngbrf(p) = ngbrlist + nnbr;
	    scanngbr(active, active+1, active+actlen-1, active+actlen-1,
		     (nodeptr) p, ngbrlist);    /* scan neighbors           */
    	} while ((Ngbrf(p) - Ngbrs(p)) < nnbr); /* loop if ngbrs not enof   */
	Hknl(p) = 0.5 * hcorrect(p, Ngbrs(p), Ngbrf(p), nnbr);
                                                /* set half of search rad.  */
    }
    free(active);                               /* release memory           */
    free(ngbrlist);
}

/*
 * THREADINIT: the same as function "threadngbr" in treengbr.c
 */

local void threadinit(nodeptr p, nodeptr n, real psize)
{
    int ndesc, i;
    nodeptr desc[NSUB+1];

    Next(p) = n;                                /* set link to next node    */
    if (Type(p) == CELL) {                      /* if descendents to thread */
        ndesc = 0;                              /* start counting them      */
        Size(p) = ((real) 0.5) * psize;         /* set half size of cell    */
        for (i = 0; i < NSUB; i++)              /* loop over all subcells   */
            if (Subp(p)[i] != NULL)             /* if this one is occupied  */
                desc[ndesc++] = Subp(p)[i];     /* then store it in table   */
        More(p) = desc[0];                      /* set link to 1st one      */
        desc[ndesc] = n;                        /* thread last one to next  */
        for (i = 0; i < ndesc; i++)             /* loop over descendents    */
            threadinit(desc[i], desc[i+1], psize/2.0);
                                                /* and thread them together */
    }
}
