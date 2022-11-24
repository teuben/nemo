/****************************************************************************/
/* DEFS.H: include file for hierarchical force calculation routines.  The   */
/* definitions in this file are needed for load.c and grav.c; this file     */
/* does not provide definitions for other parts of the N-body code.         */
/*                                                                          */
/* Copyright (c) 1993 by Joshua E. Barnes, Honolulu, HI.                    */
/* It's free because it's yours.                                            */
/****************************************************************************/
 
#include "stdinc.h"
#include "real.h"
#include "vectmath.h"
 
/*
 * Body and cell data structures are used to represent the tree.  During
 * tree construction, descendent pointers are stored in the subp arrays:
 *
 *          +-------------------------------------------------------------+
 * root --> | CELL: mass, pos, next, rcrit2, more, subp:[/,o,/,/,/,/,o,/] |
 *          +----------------------------------------------|---------|----+
 *                                                         |         |
 *     +---------------------------------------------------+         |
 *     |                                                             |
 *     |    +--------------------------------------+                 |
 *     +--> | BODY: mass, pos, next, vel, acc, phi |                 |
 *          +--------------------------------------+                 |
 *                                                                   |
 *     +-------------------------------------------------------------+
 *     |
 *     |    +-------------------------------------------------------------+
 *     +--> | CELL: mass, pos, next, rcrit2, more, subp:[o,/,/,o,/,/,o,/] |
 *          +--------------------------------------------|-----|-----|----+
 *                                                      etc   etc   etc
 *
 * After the tree is complete, it is threaded to permit linear force
 * calculation, using the next and more pointers.  The storage used for
 * the subp arrays may be reused to store quadrupole moments.
 *
 *          +-----------------------------------------------+
 * root --> | CELL: mass, pos, next:/, rcrit2, more:o, quad |
 *          +---------------------------------------|-------+
 *                                                  |
 *     +--------------------------------------------+
 *     |                             
 *     |    +----------------------------------------+
 *     +--> | BODY: mass, pos, next:o, vel, acc, phi |
 *          +-----------------------|----------------+
 *                                  |
 *     +----------------------------+
 *     |
 *     |    +-----------------------------------------------+
 *     +--> | CELL: mass, pos, next:/, rcrit2, more:o, quad |
 *          +---------------------------------------|-------+
 *                                                 etc
 */
 
/*
 * NODE: data common to BODY and CELL structures.
 */
 
typedef struct _node {
    short type;                 /* code for node type */
    real mass;                  /* total mass of node */
    vector pos;                 /* position of node */
    struct _node *next;		/* link to next force-calc */   
} node, *nodeptr;
 
#define Type(x) (((nodeptr) (x))->type)
#define Mass(x) (((nodeptr) (x))->mass)
#define Pos(x)  (((nodeptr) (x))->pos)
#define Next(x) (((nodeptr) (x))->next)

/*
 * BODY: data structure used to represent particles.
 */
 
#define BODY 01                 /* type code for bodies */
 
typedef struct {
    node bodynode;              /* data common to all nodes */
    vector vel;                 /* velocity of body */
    vector acc;                 /* acceleration of body */
    real phi;                   /* potential at body */
} body, *bodyptr;
 
#define Body    body
 
#define Vel(x)  (((bodyptr) (x))->vel)
#define Acc(x)  (((bodyptr) (x))->acc)
#define Phi(x)  (((bodyptr) (x))->phi)

/*
 * CELL: structure used to represent internal nodes of tree.
 */
 
#define CELL 02                 /* type code for cells */
 
#define NSUB (1 << NDIM)        /* subcells per cell */
 
typedef struct {
    node cellnode;              /* data common to all nodes */
    real rcrit2;                /* critical c-of-m radius^2 */
    nodeptr more;		/* link to first descendent */   
    union {			/* shared storage for... */
	nodeptr subp[NSUB];     /* descendents of cell */
	matrix quad;            /* quad. moment of cell */
    } stuff;
} cell, *cellptr;
 
#define Rcrit2(x) (((cellptr) (x))->rcrit2)
#define More(x)   (((cellptr) (x))->more)
#define Subp(x)   (((cellptr) (x))->stuff.subp)
#define Quad(x)   (((cellptr) (x))->stuff.quad)

/*
 * Variables used in tree construction.
 */
 
global cellptr root;                    /* pointer to root cell             */
 
global real rsize;                      /* side-length of root cell         */
 
global int cellused;			/* count of cells in tree           */

global int maxlevel;			/* count of levels in tree          */

/*
 * Parameters and results for gravitational calculation.
 */
 
global string options;                  /* various option keywords          */
 
global real theta;                      /* accuracy parameter: 0.0 => exact */
 
global bool usequad;		        /* use quadrupole corrections       */

global real eps;                        /* potential softening parameter    */
 
global int n2bterm;                     /* number 2-body of terms evaluated */
 
global int nbcterm;                     /* num of body-cell terms evaluated */

/*
 * Utility routines used in load.c and grav.c.  These are defined in
 * util.c, which must be compiled with the same choice of precision.
 */

bool scanopt(string, string);			/* find option in string    */
real cputime(void);				/* return elapsed CPU time  */
void *allocate(int);				/* allocate and zero memory */
real distv(vector, vector);			/* distance between vectors */
void error(string, ...);			/* report error and exit    */
void eprintf(string, ...);			/* printf to error stream   */
