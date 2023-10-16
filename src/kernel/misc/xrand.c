/*
 * XRAND: floating-point random number routine.
 *
 * Arguments: xl, xh: lower, upper bounds on random number.
 * Returns: uniformly distributed random number in [xl..xh]
 * It uses the generic unix rand() routine, which is known
 * to be not so great. See xrandom() for a  slightly better one.
 *
 *     4-mar-94    ansi, added -DNEED_RAND section - added stdin.h
 *    18-apr-2021  von Hoerner's 1957 pretty good algorithm
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


static double seed_svh = -1.0;

/*
 * ran_svh57:    returns a number between 0 and 1
 */

double ran_svh57(double seed)
{
  double x2, xnew;
  if (seed_svh < 0) {
    if (seed < 0.57 || seed > 0.91) error("SvH seed needs to be between 0.57 and 0.91");
    dprintf(0,"SvH seed=%g\n",seed);
    seed_svh = seed;
  }
  x2 = seed_svh*seed_svh;
  if (x2 < 0.57)
    seed_svh = x2 + 0.32;
  else
    seed_svh = x2;
  xnew = (seed_svh*1000) - (int)(seed_svh*1000);
  return xnew;
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
