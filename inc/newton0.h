/* newton0.h - NEARBY_PQ0h, coincide, earlier, earliest, later, latest, 
               too_short */

/*
 *  newton0.h: for a gravitational many-body code with collective time steps
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
#ifdef REGULARIZATION
#  define  NDIM  4
#endif

#include  <math.h>
#include  <stdinc.h>
#include  <vectmath.h>
#include  "state.h"

/*-----------------------------------------------------------------------------
 *  coincide, earlier, later  --  compare two points in time, which are
 *                                represented as floating point numbers
 *                                with finite accuracy. These are fuzzy
 *                                equivalents of  = , < , >  in the 
 *                                direction of integration (i.e.,
 *                                if t_end < t_begin, the roles of < and >
 *                                are reversed).
 *          earliest, latest  --  returns the earliest and latest of two
 *                                points in time, in the direction of 
 *                                integration (no fuzziness here).
 *                 too_short  --  shorthand for (coincide ... ,0.0).
 *  BEWARE: these macros use the variable forwards, which has to be defined 
 *          locally before using the macros. This is a bit of a nuisance,
 *          but it is a consequence of the decision not to use any shared
 *          global (i.e. external) variables in newton0. It does lead to larger
 *          freedom, by allowing e.g. temporary time reversals for part of the
 *          system (after time overshoot of a regularized sub-system, say).
 *          Typically the value of forwards is set locally through the macro
 *           Forwards()  acting on the control block of the appropriate state.
 *  NOTE: instead of using  NEARBY , which may be accidentally redefined in a
 *        some part of a program, the "nonsense" addition "_PQ0h" is added to
 *        minimize the chance of collisions with other macros.
 *-----------------------------------------------------------------------------
 */
#define   NEARBY_PQ0h          0.000000001

#define   earlier(t1, t2)         (forwards ? ((t1) < ((t2) - NEARBY_PQ0h)) : \
                                              ((t1) > ((t2) + NEARBY_PQ0h)))
#define   later(t1, t2)           (forwards ? ((t1) > ((t2) + NEARBY_PQ0h)) : \
                                              ((t1) < ((t2) - NEARBY_PQ0h)))
#define   earliest(t1, t2)        (forwards ? (MIN(t1, t2)) : (MAX(t1, t2)))
#define   latest(t1, t2)          (forwards ? (MAX(t1, t2)) : (MIN(t1, t2)))
#define   coincide(t1, t2)        (ABS((t1)-(t2)) < NEARBY_PQ0h)
#define   too_short(delta_t)      (ABS(delta_t) < NEARBY_PQ0h)

/* endof: newton0.h */
