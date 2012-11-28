/*
 * SPHBODY.H: Structure and accessor macro definitions for an SPH-flavored
 * body structure to be used in conjunction with SPH SnapShot binary files.
 */

typedef struct {
    real   bodymass;			/* mass of body			    */
    vector bodyphase[2];		/* phase-space coordinates	    */
    real   bodyphi;			/* gravitational potential	    */
    vector bodyacc;			/* gravitational acceleration	    */
    real   bodydens;			/* density at body position         */
    real   bodytemp;                    /* temperature                      */
    real   bodysl;                      /* smoothing lenght                 */
    real   bodyaux;			/* misc. real value assoc. w. body  */
    int    bodykey;			/* misc. int. value assoc. w. body  */
} body;

#define Body     body

#define Mass(b)  ((b)->bodymass)
#define Phase(b) ((b)->bodyphase)
#define Pos(b)   ((b)->bodyphase[0])
#define Vel(b)   ((b)->bodyphase[1])
#define Phi(b)   ((b)->bodyphi)
#define Acc(b)   ((b)->bodyacc)
#define Dens(b)  ((b)->bodydens)
#define Temp(b)  ((b)->bodytemp)
#define Sl(b)    ((b)->bodysl)
#define Aux(b)   ((b)->bodyaux)
#define Key(b)   ((b)->bodykey)
