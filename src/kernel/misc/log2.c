/*
 * LOG2: log to the base two.
 * WHY AINT THIS STANDARD? WHATS *W*R*O*N*G* WITH PEOPLE?
 *
 *	22-jun-01	added ZENO mixed precision code		PJT
 *      28-mar-05       log2 now under configure control        PJT
 * 
 *   Warning: gcc3 has log2
 *            cygwin uses a macro, and doesn't detect it via configure
 */

#include <stdinc.h>
#include <math.h>
#include <mathfns.h>

#if !defined(HAVE_LOG2)
double log2(double x)
{

    return (log(x) / 0.693147180559945309417);	/* normalize by ln(2)       */
                   /*  M_LOG2_E  */
}
#endif

#if !defined(DOUBLEPREC)
real rlog2(real x)
{
    return (rlog(x) / M_LN2);
}
#endif
