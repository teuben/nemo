/*
 * DEFS: include file for hack programs.
 */

#include <stdinc.h>
#include <vectmath.h>
#include <assert.h>

/*
 * GLOBAL: pseudo-keyword for storage class.
 */
 
#if !defined(global)
#  define global extern
#endif


/*
 * BODY and CELL data structures are used to represent the tree:
 *
 *         +-----------------------------------------------+
 * root--> | CELL: mass, pos, quad, /, o, /, /, /, /, o, / |
 *         +---------------------------|--------------|----+
 *                                     |              |
 *    +--------------------------------+              |
 *    |                                               |
 *    |    +--------------------------------+         |
 *    +--> | BODY: mass, pos, vel, acc, phi |         |
 *         +--------------------------------+         |
 *                                                    |
 *    +-----------------------------------------------+
 *    |
 *    |    +-----------------------------------------------+
 *    +--> | CELL: mass, pos, quad, o, /, /, o, /, /, o, / |
 *         +------------------------|--------|--------|----+
 *                                 etc      etc      etc
 */

/*
 * NODE: data common to BODY and CELL structures.
 */

typedef struct {
    short type;                 /* code for node type */
    real mass;                  /* total mass of node */
    vector pos;			/* position of node */
} node, *nodeptr;

#define Type(x) (((nodeptr) (x))->type)
#define Mass(x) (((nodeptr) (x))->mass)
#define Pos(x)  (((nodeptr) (x))->pos)

/*
 * BODY: data structure used to represent particles.
 */

#define BODY 01                 /* type code for bodies */

typedef struct {
    short type;
    real mass;                  /* mass of body */
    vector pos;                 /* position of body */
    vector vel;                 /* velocity of body */
    vector acc;			/* acceleration of body */
    real phi;			/* potential at body */
} body, *bodyptr;

#define Body    body
#define Vel(x)  (((bodyptr) (x))->vel)
#define Acc(x)  (((bodyptr) (x))->acc)
#define Phi(x)  (((bodyptr) (x))->phi)

/*
 * PHASEBODY: alternate definition introduced for I/O.
 */

typedef struct {
    short type;
    real mass;
    vector phase[2];            /* position, velocity of body */
    vector acc;
    real phi;
} phasebody;

#define Phase(x)  (((phasebody *) (x))->phase)

/*
 * CELL: structure used to represent internal nodes of tree.
 */

#define CELL 02                 /* type code for cells */

#define NSUB (1 << NDIM)        /* subcells per cell */

typedef struct {
    short type;
    real mass;                  /* total mass of cell */
    vector pos;                 /* cm. position of cell */
#ifdef QUADPOLE
    matrix quad;		/* quad. moment of cell */
#endif
    nodeptr subp[NSUB];         /* descendents of cell */
} cell, *cellptr;

#ifdef QUADPOLE
#define Quad(x) (((cellptr) (x))->quad)
#endif

#define Subp(x) (((cellptr) (x))->subp)


#if defined(cray)
#define IMAX (1 << 30)
#else
#define IMAX (1 << (8 * sizeof(int) - 2))       /* highest bit */
#endif

