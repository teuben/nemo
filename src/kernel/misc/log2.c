/*
 * LOG2: log to the base two.
 * WHY AINT THIS STANDARD? WHATS *W*R*O*N*G* WITH PEOPLE?
 *
 */

#include <math.h>

double log2(double x)
{

    return (log(x) / 0.693147180559945);	/* normalize by ln(2)       */
}
