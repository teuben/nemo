/*
 * BAREBODY.H: Structure and accessor macro definitions for a minimal
 * body structure to be used in conjunction with SnapShot binary files.
 */

typedef struct {
    real   bodymass;			/* mass of body			    */
    vector bodyphase[2];		/* phase-space coordinates	    */
} body;

#define Body     body

#define Mass(b)  ((b)->bodymass)
#define Phase(b) ((b)->bodyphase)
#define Pos(b)   ((b)->bodyphase[0])
#define Vel(b)   ((b)->bodyphase[1])
