/*
 * XRAND: floating-point random number routine.
 *
 * Arguments: xl, xh: lower, upper bounds on random number.
 * Returns: uniformly distributed random number in [xl..xh]
 * It uses the generic unix rand() routine, which is known
 * to be not so great. See xrandom() for a  slightly better one.
 *
 *     4-mar-94  ansi, added -DNEED_RAND section - added stdin.h
 */

#include <stdinc.h>
#include <stdlib.h>  /* this is where rand() normally resides */

#ifndef RAND_MAX
#define RAND_MAX  2147483647   /* some machines (DOS) may need 32767 */
#endif

double xrand(double xl, double xh)
{
    return xl + (xh - xl) * ((double) rand()) / RAND_MAX;
}

#if defined(NEED_RAND)
/* use -DNEED_RAND if you don't have or need the system's rand() and srand() */
/* this one we owe to: Kernigan & Ritchie; p.46 */
/* returns a pseudo-random number 0..32767 */

static unsigned long int next = 1;
int rand(void)
{
    next = next * 1103515245 + 12345;
    return (unsigned int)(next/65536) % 32768;
}

void srand(unsigned int seed)
{
    next = seed;
}
#endif
