/*
 * STDBODY.H: Structure and accessor macro definitions for a vanilla-flavored
 * body structure to be used in conjunction with SnapShot binary files.
 * An additional <extbody.h> is needed to glue the extbody....
 *
 */

typedef struct {
    real   bodymass;			/* mass of body			    */
    vector bodypos;		        /* positions            	    */
    real   bodyvel;			/* gravitational potential	    */
    char  *bodystory;                   /* body story                       */
    char  *extbody;                     /* pointer to extended body         */
} body, *bodyptr;

#define Body     body

#define Mass(b)    ((b)->bodymass)
#define Pos(b)     ((b)->bodypos)
#define Vel(b)     ((b)->bodyvel)
#define Story(b)   ((b)->bodystory)
#define Extbody(b) ((b)->extbody)

