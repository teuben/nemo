/****************************************************************************/
/* TREELOAD.C: routines to create tree.  Public routines: maketree().       */
/* Copyright (c) 1999 by Joshua E. Barnes, Tokyo, JAPAN.                    */
/****************************************************************************/

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "treedefs.h"

/*
 * Local routines and variables to perform tree construction.
 */

local void newtree(void);                       /* flush existing tree      */
local cellptr makecell(void);                   /* create an empty cell     */
local void expandbox(bodyptr, int);             /* set size of root cell    */
local void loadbody(bodyptr);                   /* load body into tree      */
local int subindex(bodyptr, cellptr);           /* compute subcell index    */
local void hackcofm(cellptr, real, int);        /* find centers of mass     */
local void setrcrit(cellptr, vector, real);     /* set cell's crit. radius  */
local void threadtree(nodeptr, nodeptr);        /* set next and more links  */
local void hackquad(cellptr);                   /* compute quad moments     */

local bool bh86, sw94;                          /* use alternate criteria   */
local nodeptr freecell = NULL;                  /* list of free cells       */

#define MAXLEVEL  32                            /* max height of tree       */

local int cellhist[MAXLEVEL];                   /* count cells by level     */
local int subnhist[MAXLEVEL];                   /* count subnodes by level  */

/*
 * MAKETREE: initialize tree structure for hierarchical force calculation.
 */

void maketree(bodyptr btab, int nbody)
{
    double cpustart;
    bodyptr p;
    int i;

    cpustart = cputime();                       /* record time at start     */
    newtree();                                  /* flush existing tree, etc */
    root = makecell();                          /* allocate the root cell   */
    CLRV(Pos(root));                            /* initialize the midpoint  */
    expandbox(btab, nbody);                     /* and expand cell to fit   */
    for (p = btab; p < btab+nbody; p++)         /* loop over all bodies     */
        loadbody(p);                            /* insert each into tree    */
    bh86 = scanopt(options, "bh86");            /* set flags for alternate  */
    sw94 = scanopt(options, "sw94");            /* ...cell opening criteria */
    if (bh86 && sw94)                           /* can't have both at once  */
        error("maketree: incompatible options bh86 and sw94\n");
    tdepth = 0;                                 /* init count of levels     */
    for (i = 0; i < MAXLEVEL; i++)              /* and init tree histograms */
        cellhist[i] = subnhist[i] = 0;
    hackcofm(root, rsize, 0);                   /* find c-of-m coords, etc  */
    threadtree((nodeptr) root, NULL);           /* add next and more links  */
    if (usequad)                                /* if including quad terms  */
        hackquad(root);                         /* find quadrupole moments  */
    cputree = cputime() - cpustart;             /* store elapsed CPU time   */
}

/*
 * NEWTREE: reclaim cells in tree, prepare to build new one.
 */

local void newtree(void)
{
    static bool firstcall = TRUE;
    nodeptr p;

    if (! firstcall) {                          /* if cells to reclaim      */
        p = (nodeptr) root;                     /* start with the root      */
        while (p != NULL)                       /* loop scanning tree       */
            if (Type(p) == CELL) {              /* if we found a cell to    */
                Next(p) = freecell;             /* then save existing list  */
                freecell = p;                   /* and add it to the front  */
                p = More(p);                    /* then scan down tree      */
            } else                              /* else, skip over bodies   */
                p = Next(p);                    /* by going on to the next  */
    } else                                      /* else nothing to reclaim  */
        firstcall = FALSE;                      /* so just note it          */
    root = NULL;                                /* flush existing tree      */
    ncell = 0;                                  /* reset cell count         */
}

/*
 * MAKECELL: return pointer to free cell.
 */

local cellptr makecell(void)
{
    cellptr c;
    int i;

    if (freecell == NULL)                       /* if no free cells left    */
        c = (cellptr) allocate(sizeof(cell));   /* then allocate a new one  */
    else {                                      /* else use existing cell   */
        c = (cellptr) freecell;                 /* take the one in front    */
        freecell = Next(c);                     /* and go on to next one    */
    }
    Type(c) = CELL;                             /* initialize node type     */
    Update(c) = FALSE;                          /* and force update flag    */
    for (i = 0; i < NSUB; i++)                  /* loop over subcells       */
        Subp(c)[i] = NULL;                      /* and empty each one       */
    ncell++;                                    /* count one more cell      */
    return (c);                                 /* return pointer to cell   */
}

/*
 * EXPANDBOX: find range of coordinate values (with respect to root)
 * and expand root cell to fit.  The size is doubled at each step to
 * take advantage of exact representation of powers of two.
 */

local void expandbox(bodyptr btab, int nbody)
{
    real dmax, d;
    bodyptr p;
    int k;

    dmax = 0.0;                                 /* keep track of max value  */
    for (p = btab; p < btab+nbody; p++)         /* loop over all bodies     */
        for (k = 0; k < NDIM; k++) {            /* and over all dimensions  */
            d = rabs(Pos(p)[k] - Pos(root)[k]); /* find distance to midpnt  */
            if (d > dmax)                       /* if bigger than old one   */
                dmax = d;                       /* store new max value      */
        }
    while (rsize < 2 * dmax)                    /* loop until value fits    */
        rsize = 2 * rsize;                      /* doubling box each time   */
}

/*
 * LOADBODY: descend tree and insert body p in appropriate place.
 */

local void loadbody(bodyptr p)
{
    cellptr q, c;
    int qind, k;
    real qsize, dist2;
    vector distv;

    q = root;                                   /* start with tree root     */
    qind = subindex(p, q);                      /* get index of subcell     */
    qsize = rsize;                              /* keep track of cell size  */
    while (Subp(q)[qind] != NULL) {             /* loop descending tree     */
        if (Type(Subp(q)[qind]) == BODY) {      /* if another body reached  */
            DOTPSUBV(dist2, distv, Pos(p), Pos(Subp(q)[qind]));
            if (dist2 == 0.0)                   /* check positions differ   */
                error("loadbody: two bodies have same position\n");
            c = makecell();                     /* then allocate new cell   */
            for (k = 0; k < NDIM; k++)          /* and initialize midpoint  */
                Pos(c)[k] = Pos(q)[k] +         /* offset from parent       */
                    (Pos(p)[k] < Pos(q)[k] ? - qsize : qsize) / 4;
            Subp(c)[subindex((bodyptr) Subp(q)[qind], c)] = Subp(q)[qind];
                                                /* put body in cell         */
            Subp(q)[qind] = (nodeptr) c;        /* link cell in tree        */
        }
        q = (cellptr) Subp(q)[qind];            /* advance to next level    */
        qind = subindex(p, q);                  /* get index to examine     */
        qsize = qsize / 2;                      /* shrink current cell      */
    }
    Subp(q)[qind] = (nodeptr) p;                /* found place, store p     */
}

/*
 * SUBINDEX: compute subcell index for body p in cell q.
 */

local int subindex(bodyptr p, cellptr q)
{
    int ind, k;

    ind = 0;                                    /* accumulate subcell index */
    for (k = 0; k < NDIM; k++)                  /* loop over dimensions     */
        if (Pos(q)[k] <= Pos(p)[k])             /* if beyond midpoint       */
            ind += NSUB >> (k + 1);             /* then skip over subcells  */
    return (ind);
}

/*
 * HACKCOFM: descend tree finding center-of-mass coordinates;
 * also sets critical cell radii, if appropriate.
 */

local void hackcofm(cellptr p, real psize, int lev)
{
    vector cmpos, tmpv;
    int i, k;
    nodeptr q;
    real dpq;

    tdepth = MAX(tdepth, lev);                  /* remember maximum level   */
    cellhist[lev]++;                            /* count cells by level     */
    Mass(p) = 0.0;                              /* init cell's total mass   */
    CLRV(cmpos);                                /* and center of mass pos   */
    for (i = 0; i < NSUB; i++)                  /* loop over the subnodes   */
        if ((q = Subp(p)[i]) != NULL) {         /* skipping the NULLs       */
            subnhist[lev]++;                    /* count existing subnodes  */
            if (Type(q) == CELL)                /* and if node is a cell    */
                hackcofm((cellptr) q, psize/2, lev+1);
                                                /* then do the same for it  */
            Update(p) |= Update(q);             /* propagate update request */
            Mass(p) += Mass(q);                 /* accumulate total mass    */
            MULVS(tmpv, Pos(q), Mass(q));       /* weight position by mass  */
            ADDV(cmpos, cmpos, tmpv);           /* and sum c-of-m position  */
        }
    if (Mass(p) > 0.0) {                        /* usually, cell has mass   */
        DIVVS(cmpos, cmpos, Mass(p));           /* so find c-of-m position  */
    } else {                                    /* but if no mass inside    */
        SETV(cmpos, Pos(p));                    /* use geo. center for now  */
    }
    for (k = 0; k < NDIM; k++)                  /* check c-of-m of cell     */
        if (cmpos[k] < Pos(p)[k] - psize/2 ||   /* if actually outside cell */
              Pos(p)[k] + psize/2 <= cmpos[k])
            error("hackcofm: tree structure error\n");
#if !defined(QUICKSCAN)
    setrcrit(p, cmpos, psize);                  /* set critical radius      */
#endif
    SETV(Pos(p), cmpos);                        /* and center-of-mass pos   */
}

#if !defined(QUICKSCAN)

/*
 * SETRCRIT: assign critical radius for cell p, using center-of-mass
 * position cmpos and cell size psize.
 */

local void setrcrit(cellptr p, vector cmpos, real psize)
{
    real bmax2, d;
    int k;

    if (theta == 0.0)                           /* if exact calculation     */
        Rcrit2(p) = rsqr(2 * rsize);            /* then always open cells   */
    else if (sw94) {                            /* if using S&W's criterion */
        bmax2 = 0.0;                            /* compute max distance^2   */
        for (k = 0; k < NDIM; k++) {            /* loop over dimensions     */
            d = cmpos[k] - Pos(p)[k] + psize/2; /* get dist from corner     */
            bmax2 += rsqr(MAX(d, psize - d));   /* and sum max distance^2   */
        }
        Rcrit2(p) = bmax2 / rsqr(theta);        /* use max dist from cm     */
    } else if (bh86)                            /* if using old criterion   */
        Rcrit2(p) = rsqr(psize / theta);        /* then use size of cell    */
    else {                                      /* else use new criterion   */
        DISTV(d, cmpos, Pos(p));                /* find offset from center  */
        Rcrit2(p) = rsqr(psize / theta + d);    /* use size plus offset     */
    }
}

#endif

/*
 * THREADTREE: do a recursive treewalk starting from node p,
 * with next stop n, installing Next and More links.
 */

local void threadtree(nodeptr p, nodeptr n)
{
    int ndesc, i;
    nodeptr desc[NSUB+1];

    Next(p) = n;                                /* set link to next node    */
    if (Type(p) == CELL) {                      /* if descendents to thread */
        ndesc = 0;                              /* start counting them      */
        for (i = 0; i < NSUB; i++)              /* loop over all subcells   */
            if (Subp(p)[i] != NULL)             /* if this one is occupied  */
                desc[ndesc++] = Subp(p)[i];     /* then store it in table   */
        More(p) = desc[0];                      /* set link to 1st one      */
        desc[ndesc] = n;                        /* thread last one to next  */
        for (i = 0; i < ndesc; i++)             /* loop over descendents    */
            threadtree(desc[i], desc[i+1]);     /* and thread them together */
    }
}

/*
 * HACKQUAD: descend tree, evaluating quadrupole moments.  Note that this
 * routine is coded so that the Subp() and Quad() components of a cell can
 * share the same memory locations.
 */

local void hackquad(cellptr p)
{
    int ndesc, i;
    nodeptr desc[NSUB], q;
    vector dr;
    real drsq;
    matrix drdr, Idrsq, tmpm;

    ndesc = 0;                                  /* count occupied subnodes  */
    for (i = 0; i < NSUB; i++)                  /* loop over all subnodes   */
        if (Subp(p)[i] != NULL)                 /* if this one's occupied   */
            desc[ndesc++] = Subp(p)[i];         /* copy it to safety        */
    CLRM(Quad(p));                              /* init quadrupole moment   */
    for (i = 0; i < ndesc; i++) {               /* loop over real subnodes  */
        q = desc[i];                            /* access ech one in turn   */
        if (Type(q) == CELL)                    /* if it's also a cell      */
            hackquad((cellptr) q);              /* then process it first    */
        SUBV(dr, Pos(q), Pos(p));               /* find displacement vect.  */
        OUTVP(drdr, dr, dr);                    /* form outer prod. of dr   */
        DOTVP(drsq, dr, dr);                    /* and dot prod. dr * dr    */
        SETMI(Idrsq);                           /* init unit matrix         */
        MULMS(Idrsq, Idrsq, drsq);              /* and scale by dr * dr     */
        MULMS(tmpm, drdr, 3.0);                 /* scale drdr by 3          */
        SUBM(tmpm, tmpm, Idrsq);                /* now form quad. moment    */
        MULMS(tmpm, tmpm, Mass(q));             /* from cm of subnode       */
        if (Type(q) == CELL)                    /* if subnode is cell       */
            ADDM(tmpm, tmpm, Quad(q));          /* then include its moment  */
        ADDM(Quad(p), Quad(p), tmpm);           /* increment moment of cell */
    }
}
