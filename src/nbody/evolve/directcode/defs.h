/*
 * DEFS: include file for hack programs.
 */

#include <stdinc.h>
#include <vectmath.h>

/*
 * GLOBAL: pseudo-keyword for storage class.
 */
  
#if !defined(global)
#  define global extern
#endif
 

/*
 * BODY: data structure used to represent particles.
 */

#define BODY 01                 /* type code for bodies */

typedef struct {
    real mass;                  /* mass of body */
    vector pos;                 /* position of body */
    vector vel;                 /* velocity of body */
    vector acc;			/* acceleration of body */
    real phi;			/* potential at body */
} body, *bodyptr;

#define Body    body
#define Mass(x) (((bodyptr) (x))->mass)
#define Pos(x)  (((bodyptr) (x))->pos)
#define Vel(x)  (((bodyptr) (x))->vel)
#define Acc(x)  (((bodyptr) (x))->acc)
#define Phi(x)  (((bodyptr) (x))->phi)

/*
 * PHASEBODY: alternate definition introduced for I/O.
 */

typedef struct {
    real mass;
    vector phase[2];            /* position, velocity of body */
    vector acc;
    real phi;
} phasebody;

#define Phase(x)  (((phasebody *) (x))->phase)

