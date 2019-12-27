/*
 * BODY.H: Structure and accessor macro definitions for a vanilla-flavored
 * body structure to be used in conjunction with SnapShot binary files.
 */
#ifndef _body_h
#define _body_h


#define _body_h_dens

typedef struct {
    real   mass;		/* mass of body			    */
    vector phase[2];		/* phase-space coordinates	    */
    real   phi;			/* gravitational potential	    */
    vector acc;			/* gravitational acceleration	    */
    real   aux;			/* misc. real value assoc. w. body  */
    int    key;			/* misc. int. value assoc. w. body  */
#ifdef _body_h_dens
    real   dens;		/* density associated w. body       */
    real   eps;                 /* softening length w. body         */
#endif
} body;

typedef int  (*btiproc)(body *, real, int);
typedef real (*btrproc)(body *, real, int);


#define Body     body

#define Mass(b)  ((b)->mass)
#define Phase(b) ((b)->phase)
#define Pos(b)   ((b)->phase[0])
#define Vel(b)   ((b)->phase[1])
#define Phi(b)   ((b)->phi)
#define Acc(b)   ((b)->acc)
#define Aux(b)   ((b)->aux)
#define Key(b)   ((b)->key)
#ifdef _body_h_dens
#define Dens(b)  ((b)->dens)
#define Eps(b)   ((b)->eps)
#endif

#endif
