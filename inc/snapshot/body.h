/*
 * BODY.H: Structure and accessor macro definitions for a vanilla-flavored
 * body structure to be used in conjunction with SnapShot binary files.
 */
#ifndef _body_h
#define _body_h


#define _body_h_dens

typedef struct {
    real   _mass;		/* mass of body			    */
    vector _phase[2];		/* phase-space coordinates	    */
    real   _phi;		/* gravitational potential	    */
    vector _acc;	 	/* gravitational acceleration	    */
    real   _aux;		/* misc. real value assoc. w. body  */
    int    _key;		/* misc. int. value assoc. w. body  */
#ifdef _body_h_dens
    real   _dens;		/* density associated w. body       */
    real   _eps;                /* softening length w. body         */
#endif
} body;

typedef int  (*btiproc)(body *, real, int);
typedef real (*btrproc)(body *, real, int);


#define Body     body

#define Mass(b)  ((b)->_mass)
#define Phase(b) ((b)->_phase)
#define Pos(b)   ((b)->_phase[0])
#define Vel(b)   ((b)->_phase[1])
#define Phi(b)   ((b)->_phi)
#define Acc(b)   ((b)->_acc)
#define Aux(b)   ((b)->_aux)
#define Key(b)   ((b)->_key)
#ifdef _body_h_dens
#define Dens(b)  ((b)->_dens)
#define Eps(b)   ((b)->_eps)
#endif

#endif
