/*
 * BODY.H: Structure and accessor macro definitions for a vanilla-flavored
 * body structure to be used in conjunction with SnapShot binary files.
 */
#ifndef _body_h
#define _body_h


#define _body_h_dens

typedef struct {
    real   bodymass;			/* mass of body			    */
    vector bodyphase[2];		/* phase-space coordinates	    */
    real   bodyphi;			/* gravitational potential	    */
    vector bodyacc;			/* gravitational acceleration	    */
    real   bodyaux;			/* misc. real value assoc. w. body  */
    int    bodykey;			/* misc. int. value assoc. w. body  */
#ifdef _body_h_dens
    real   bodydens;			/* density associated w. body       */
    real   bodyeps;                     /* softening length w. body         */
#endif
} body;

typedef int  (*btiproc)(body *, real, int);
typedef real (*btrproc)(body *, real, int);


#define Body     body

#define Mass(b)  ((b)->bodymass)
#define Phase(b) ((b)->bodyphase)
#define Pos(b)   ((b)->bodyphase[0])
#define Vel(b)   ((b)->bodyphase[1])
#define Phi(b)   ((b)->bodyphi)
#define Acc(b)   ((b)->bodyacc)
#define Aux(b)   ((b)->bodyaux)
#define Key(b)   ((b)->bodykey)
#ifdef _body_h_dens
#define Dens(b)  ((b)->bodydens)
#define Eps(b)   ((b)->bodyeps)
#endif

#endif
