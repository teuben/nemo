/*
 * MYBODY.H: Structure and accessor macro definitions for a vanilla-flavored
 * body structure to be used in conjunction with SnapShot binary files.
 * 
 * This is a standard body with stories added to it
 */

typedef struct {
    real   bodymass;			/* mass of body			    */
    vector bodyphase[2];		/* phase-space coordinates	    */
    real   bodyphi;			/* gravitational potential	    */
    vector bodyacc;			/* gravitational acceleration	    */
    real   bodyaux;			/* misc. real value assoc. w. body  */
    int    bodykey;			/* misc. int. value assoc. w. body  */
    char  *bodystory;			/* misc. textual info assoc. w. body*/
} body;

#define Body     body

#define Mass(b)  ((b)->bodymass)
#define Phase(b) ((b)->bodyphase)
#define Pos(b)   ((b)->bodyphase[0])
#define Vel(b)   ((b)->bodyphase[1])
#define Phi(b)   ((b)->bodyphi)
#define Acc(b)   ((b)->bodyacc)
#define Aux(b)   ((b)->bodyaux)
#define Key(b)   ((b)->bodykey)
#define Story(b) ((b)->bodystory)
