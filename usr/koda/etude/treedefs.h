/****************************************************************************/
/* TREEDEFS.H: include file for hierarchical force calculation routines.    */
/* Copyright (c) 1999 by Joshua E. Barnes, Tokyo, JAPAN,                    */
/*    revised for SPH calculation by Jin Koda, Tokyo, JAPAN. 2000.11        */
/****************************************************************************/

#ifndef _treedefs_h
#define _treedefs_h

/* fix for NEMO, once we've sorted out the PRECISION flag stuff ZENO uses */
#if 0
#define rsqrt sqrt
#define rabs  abs
#define rpow  pow
#define rsqr  sqr
#endif


/* in NEMO there is a function 'root', so root needs to be redefined */

#define root treecode_root

/*
 * NODE: data common to BODY and CELL structures.
 */

typedef struct _node {
    short type;                 /* code for node type */
    bool update;                /* status in force calc */
    real mass;                  /* total mass of node */
    vector pos;                 /* position of node */
    real radius;                /* some kinds of radius */
    struct _node *next;         /* link to next force calc */
} node, *nodeptr;

#define Type(x)   (((nodeptr) (x))->type)
#define Update(x) (((nodeptr) (x))->update)
#define Mass(x)   (((nodeptr) (x))->mass)
#define Pos(x)    (((nodeptr) (x))->pos)
#define Next(x)   (((nodeptr) (x))->next)

#define BODY 01                 /* type code for bodies */
#define CELL 02                 /* type code for cells */

/*
 * BODY: data structure used to represent SPH particles.
 */

typedef struct _body {
    node bodynode;              /* data common to all nodes */
    vector vel;                 /* velocity of body */
    vector velm;                /* velocity of body at midpnt of step */
    vector acc;                 /* acceleration of body */
    real phi;                   /* potential at body */
    real rho;                   /* density of body */
    real csnd;                  /* internal energy */
    real divv;                  /* divergence of velocity */
#if defined(THREEDIM)
    vector rotv;                /* rotation of velocity */
#else
    real rotv;                  /* rotation of velocity */
#endif
    real vsprs;                 /* viscousity suppresing factor */
    struct _body **ngbrstart;   /* pointer to top of neighbor list */
    struct _body **ngbrfinish;  /* pointer to end of neighbor list */
} body, *bodyptr;

#define Vel(x)     (((bodyptr) (x))->vel)
#define Velm(x)    (((bodyptr) (x))->velm)
#define Acc(x)     (((bodyptr) (x))->acc)
#define Phi(x)     (((bodyptr) (x))->phi)
#define Hknl(x)    (((nodeptr) (x))->radius)
#define Rho(x)     (((bodyptr) (x))->rho)
#define Csnd(x)    (((bodyptr) (x))->csnd)
#define Divv(x)    (((bodyptr) (x))->divv)
#define Rotv(x)    (((bodyptr) (x))->rotv)
#define Vsprs(x)   (((bodyptr) (x))->vsprs)
#define Ngbrs(x)   (((bodyptr) (x))->ngbrstart)
#define Ngbrf(x)   (((bodyptr) (x))->ngbrfinish)

/*
 * CELL: structure used to represent internal nodes of tree.
 */

#define NSUB (1 << NDIM)        /* subcells per cell */

typedef struct {
    node cellnode;              /* data common to all nodes */
    nodeptr more;               /* link to first descendent */
    union {
        nodeptr subp[NSUB];     /* descendents of cell */
        matrix quad;            /* quad. moment of cell */
    } sorq;
} cell, *cellptr;

#if !defined(QUICKSCAN)
#define Rcrit2(x) (((nodeptr) (x))->radius)
#endif
#define More(x)   (((cellptr) (x))->more)
#define Subp(x)   (((cellptr) (x))->sorq.subp)
#define Quad(x)   (((cellptr) (x))->sorq.quad)

/*
 * GLOBAL: pseudo-keyword for storage class.
 */

#if !defined(global)
#  define global extern
#endif

/*
 * Parameters for tree construction and force calculation.
 */

#if !defined(QUICKSCAN)
global real theta;                      /* force accuracy parameter         */
#endif
global string options;                  /* various option keywords          */
global bool selfgrav;                   /* compute self-gravity             */
global bool usequad;                    /* use quadrupole corrections       */
global real eps;                        /* density smoothing parameter      */


/*
 * Tree construction.
 */

void planttree(bodyptr, int);           /* construct tree structure         */
void manegetree(void);                  /* comp c-of-m and make branch link */

global cellptr root;                    /* pointer to root cell             */
global real rsize;                      /* side-length of root cell         */
global int ncell;                       /* count of cells in tree           */
global int tdepth;                      /* count of levels in tree          */
global real cputree;                    /* CPU time to build tree           */

/*
 * Neighbor search.
 */

void makengbrlist(bodyptr, int);        /* construct neighbor list for SPH  */
void scanngbr(nodeptr *, nodeptr *, nodeptr *, nodeptr *, nodeptr, bodyptr *);
                                        /* scan neighbors for a particle    */

global int actcmax;                     /* max length of active cell list   */
global int actbmax;                     /* max length of active body list   */
global real cpungbr;                    /* CPU time to search neighbors     */
global bodyptr *ngbrlist;               /* neighbor list                    */
global int ngbrlen;                     /* length of neighbor list          */

/*
 * SPH calculation.
 */

void sphcalc(bodyptr, int);             /* update SPH force                 */
void stephknl(bodyptr, int, real);      /* update SPH smoothing length      */
bodyptr *checkhknl(nodeptr p, nodeptr *nptr, nodeptr *mptr);
                                        /* check if suitable ngbrs in 2h    */
real hcorrect(bodyptr, bodyptr *, bodyptr *, int);
                                        /* correct h if not in range        */

global real fgas;                       /* fraction of gas in mass          */
global real alpha;                      /* coefficient of viscousity        */
global real beta;                       /* coefficient of viscousity        */
global int nnbr;                        /* number of ngbrs requested        */
global int nmax;                        /* maximum num. of ngbrs acceped    */
global int nmin;                        /* minimum num. of ngbrs acceped    */
global real cpudens;                    /* CPU time for SPH density calc    */
global real cpusphf;                    /* CPU time for SPH force calc      */

/*
 * Gravity calculation.
 */

void gravcalc(void);                    /* update gravi. force on bodies    */

global int actmax;                      /* maximum length of active list    */
global int nbbcalc;                     /* total body-body interactions     */
global int nbccalc;                     /* total body-cell interactions     */
global real cpugrav;                    /* CPU time for force calc          */

/*
 * External force calculation.
 */

void initextf(void);                    /* initialize external force param. */
void extforce(bodyptr, int, real);      /* compute external force           */
real massext(real);                     /* total mass of ext. pot. within r */

global real mode;                       /* spiral mode; m=2 for two arms    */
global real pitch;                      /* pitch angle of stellar spiral    */


global real rinit;                      /* radius of initial disk           */
global real rcore;                      /* core radius of disk potential    */
global real vmax;                       /* maximum velocity of rotation     */
global real omgb;                       /* patarn speed of bar              */
global real fbar;                       /* strength of bar potential        */

global real rcorebh;                    /* core radius of central black hole*/
global real fbh;                        /* potential fraction of central BH */

/*
 * Make testdata.
 */

void testdata(void);                    /* generate test data               */

/*
 * Plot data.
 */

void plot(bodyptr, int);

#endif /* ! _treedefs_h */
