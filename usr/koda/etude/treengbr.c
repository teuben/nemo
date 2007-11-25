/****************************************************************************/
/* TREENGBR.C: routines to make neigbor list.                               */
/* Public routines: makengbrlist()                                          */
/* Copyright (c) 2000 by Jin Koda, Tokyo, JAPAN.                            */
/****************************************************************************/

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "treecode.h"

/*
 * Local routines and variables to perform neighbor search.
 */

local void setrsrch(cellptr, real, int);        /* set critical search size */
local void threadngbr(nodeptr, nodeptr, real);  /* set next and more links  */
local void walktreengbr(nodeptr *, nodeptr *, nodeptr *, nodeptr *, nodeptr, int);
                                                /* walk tree to search nbr. */
local void pickngbr(nodeptr *, nodeptr *, nodeptr, bodyptr *);
                                                /* pick up real neighbors   */
local bool rejectcell(nodeptr, nodeptr);        /* cell rejection criterion */
local bool rejectbody(nodeptr, nodeptr);        /* body rejection criterion */

#define MAXLEVEL  32                            /* max height of tree       */

#define Maxh(x)   (((nodeptr) (x))->mass)
#define Size(x)   (((nodeptr) (x))->mass)
#define Rsrch2(x) (((nodeptr) (x))->radius)

#if !defined(FACTIVE)
# define FACTIVE 2.0                            /* active list fudge factor */
#endif

local int cellhist[MAXLEVEL];                   /* count cells by level     */
local int subnhist[MAXLEVEL];                   /* count subnodes by level  */
local int actlen;                               /* length as allocated      */ 
local nodeptr *active;                          /* list of nodes tested     */
local bodyptr *ngbrptr;                         /* pointer to ngbr list     */

/*
 * MAKENGBRLIST: make neighbor list for SPH calculation with tree structure.
 */

void makengbrlist(bodyptr btab, int nbody)
{
    double cpustart;
    bodyptr p;
    int i;

    cpustart = cputime();                       /* record time at start     */
    tdepth = 0;                                 /* init count of levels     */
    for (i = 0; i < MAXLEVEL; i++)              /* and init tree histograms */
        cellhist[i] = subnhist[i] = 0;
    for (p=btab; p<btab+nbody; p++)             /* search radius is twice   */
	Hknl(p) *= 2.0;                         /* of smoothing length      */
    setrsrch(root, rsize, 0);                   /* set search radii for SPH */
    threadngbr((nodeptr) root, NULL, rsize);    /* add next and more links  */
    ngbrptr = ngbrlist;                         /* neighbor list            */
    actlen = FACTIVE * 512 * tdepth;            /* estimate list length     */
    active = (nodeptr *) allocate(actlen * sizeof(nodeptr));
    actmax = actbmax = actcmax = 0;             /* zero cumulative counters */
    active[0] = (nodeptr) root;                 /* initialize active list   */
    walktreengbr(active, active+1, active+actlen-1, active+actlen-1,
		 (nodeptr) root, 0);            /* perform neighbor search  */
    for (p=btab; p<btab+nbody; p++)             /* search radius is twice   */
     	Hknl(p) *= 0.5;                         /* of smoothing length      */
    free(active);                               /* release memory           */
    cpungbr = cputime() - cpustart;             /* store elapsed CPU time   */
}

/*
 * SETRSRCH: assign neighbor search box size for cell p; descending
 * tree to find maximum hknl values in cells, and seting critical
 * search size, if appropriate.
 */

local void setrsrch(cellptr p, real psize, int lev)
{
    int i;
    nodeptr q;
    real maxhknl;

    tdepth = MAX(tdepth, lev);                  /* remember maximum level   */
    cellhist[lev]++;                            /* count cells by level     */
    Maxh(p) = 0.0;                              /* init cell's total mass   */
    for (i = 0; i < NSUB; i++)                  /* loop over the subnodes   */
        if ((q = Subp(p)[i]) != NULL) {         /* skipping the NULLs       */
            subnhist[lev]++;                    /* count existing subnodes  */
            if (Type(q) == CELL) {              /* and if node is a cell    */
                setrsrch((cellptr) q, psize/2, lev+1);
                                                /* then do the same for it  */
		Maxh(p) = MAX(Maxh(p), Maxh(q));/* and find max. of hknl    */
	    } else {                            /* if node is a body        */
		Maxh(p) = MAX(Maxh(p), Hknl(q));/* find max. of hknl        */
	    }
            Update(p) |= Update(q);             /* propagate update request */
        }
    Rsrch2(p) = ((real) 0.5) * psize + Maxh(p); /* set critical size        */
}

/*
 * THREADNGBR: do a recursive treewalk starting from node p,
 * with next stop n, installing Next and More links, and
 * set half the sizes of cells.
 */

local void threadngbr(nodeptr p, nodeptr n, real psize)
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
            threadngbr(desc[i], desc[i+1], psize/2.0);
                                                /* and thread them together */
    }
}


/*
 * WALKTREENGBR: search neighbors with a complete walk of the tree.
 */

local void walktreengbr(nodeptr *aptr, nodeptr *nptr, nodeptr *bptr,
			nodeptr *mptr, nodeptr p, int lev)
{
    nodeptr *ap, *np, *bp, *mp, q;
    int actsafe;

    if (Update(p)) {                            /* new ngbr list needed?    */
	mp = mptr;                              /* start active body list   */
	for (bp = bptr; bp > mptr; bp--)        /* loop over body list      */
	    if (!rejectbody(*bp, p))            /* if body is not rejected  */
		*mp-- = *bp;                    /* forward in body list     */
	np = nptr;                              /* start active node list   */
	actsafe = NSUB - 1;                     /* leave room for NSUB more */
	for (ap = aptr; ap < nptr; ap++) {      /* loop over active nodes   */
	    if (!rejectcell(*ap, p)) {          /* if the cell not rejected */
		if (mp-np < actsafe)            /* check list has room      */
		    error("walktreengbr: active list overflow\n");
		for (q = More(*ap); q != Next(*ap); q = Next(q))
                                                /* loop over sub nodes      */
		    if (Type(q) == CELL) {      /* if node is a cell        */
			if (!rejectcell(q, p))  /* and unless rejected      */
			    *np++ = q;          /* put on active node list  */
		    } else {                    /* else if node is a body   */
			if (!rejectbody(q, p))  /* unless rejected          */
			    *mp-- = q;          /* put on active body list  */ 
		    }
	    }
	}
	actcmax = MAX(actcmax, np-active);      /* keep track of maximum    */
	actbmax = MAX(actbmax,active+actlen-mp);/* active cell and body     */
	if (Type(p) == CELL)                    /* if node is a cell        */
	    for (q = More(p); q != Next(p); q = Next(q))
		walktreengbr(nptr, np, mptr, mp, q, lev+1);
                                                /* decend tree to next lev. */
	else {                                  /* else node is a body      */
	    scanngbr(nptr, np, mptr, mp, p, ngbrptr);
                                                /* scan ngbr from actives   */
	    ngbrptr = checkhknl(p, np, mp);     /* check if h is in range   */
	}
    }
}

/*
 * SCANNGBR: iterative routine to scan neighbors for a particle "p".
 */

void scanngbr(nodeptr *aptr, nodeptr *nptr,
	      nodeptr *bptr, nodeptr *mptr, nodeptr p, bodyptr *ngbrs)
{
    nodeptr *ap, *mp, q;

    mp = mptr;                              /* start list for new active*/
    for (ap = aptr; ap < nptr; ap++) {      /* loop over active nodes   */
	q = *ap;                            /* start from top of active */
	while (q != Next(*ap))              /* loop over all decendants */
	    if (Type(q) == CELL) {          /* if it is a cell          */
		if (rejectcell(q, p))       /* and rejected             */  
		    q = Next(q);            /* then step to next node   */
		else                        /* else not rejected        */
		    q = More(q);            /* then decend to children  */
	    } else {                        /* if it is body            */
		if (!rejectbody(q, p)) {    /* and not rejected         */
		    if (mp-nptr < 1)        /* check list has room      */
			error("scanngbr: active list overflow\n");
		    *mp-- = q;              /* put on active body list  */
		}
		q = Next(q);                /* step to next node        */
	    }
    }
    pickngbr(bptr, mp, p, ngbrs);           /* pick up real ngbr        */
}

/*
 * PICKNGBR: Pick up real neighbors from potential neighbor list.
 */

local void pickngbr(nodeptr *bptr, nodeptr *mptr, nodeptr p, bodyptr *ngbrs)
{
    nodeptr *q;
    real dr2, hknl2;
    vector dr;
    int safelen;

    safelen = ngbrlen - (int) (ngbrs-ngbrlist); /* leave room for new list  */
    if ((int) (bptr - mptr) > safelen)          /* check ngbr list has room */
	error("pickngbr: ngbr list overflow\n");
    hknl2 = rsqr( Hknl(p) );                    /* compute square of h      */
    ngbrptr = ngbrs;                            /* start neighbor list      */
    for (q = bptr; q > mptr; q--)               /* loop over potential ngbr */
	if (p != *q) {
	    DOTPSUBV(dr2, dr, Pos(*q), Pos(p)); /* compute separation       */
	    if (dr2 < hknl2)                    /* if body is in h          */
		*ngbrptr++ =  (bodyptr) *q;     /* put it on ngbr list      */
	}
    Ngbrs(p) = ngbrs;                           /* set start and            */
    Ngbrf(p) = ngbrptr;                         /* end of neghbor list      */
}

/*
 * REJECTCELL: criterion rejects any cell not touching searching box of node p.
 */

local bool rejectcell(nodeptr c, nodeptr p)
{
    real dk;
    int k;

    for (k=0; k < NDIM; k++) {
	dk = ABS(Pos(c)[k] -Pos(p)[k])-Size(c); /* find distance to midpnt  */
	if (dk > Rsrch2(p))                     /* if c is out of boundary  */
	    return(TRUE);                       /* then reject it           */
    }
    return(FALSE);
}

/*
 * REJECTBODY: criterion rejects any body not touching searching box of node p.
 */

local bool rejectbody(nodeptr c, nodeptr p)
{
    real dk;
    int k;

    for (k=0; k < NDIM; k++) {
	dk = ABS(Pos(c)[k] -Pos(p)[k]);         /* find distance to midpnt  */
	if (dk > Rsrch2(p))                     /* if c is out of boundary  */
	    return(TRUE);                       /* then reject it           */
    }
    return(FALSE);
}
