/*
 * FORCE_MATH: this curious function is never called.  However, its mere
 * existence is sufficient to fool ld into loading the math functions.
 * Useful for routines which use the dynamic object loader.
 */

#include <stdinc.h>
              
static void force_math()
{
    error("mathlinker should never be called");
    
    (void) cbrt(1.0);
    (void) sqrt(1.0);
    (void) qbe(1.0);
    (void) sqr(1.0);
    (void) sin(1.0);
    (void) cos(1.0);
    (void) asin(1.0);
    (void) acos(1.0);
    (void) tan(1.0);
    (void) atan(1.0);
    (void) atan2(1.0, 1.0);
    (void) exp(1.0);
    (void) dex(1.0);
    (void) log(1.0);
    (void) log10(1.0);
    (void) pow(1.0, 1.0);
    (void) fabs(1.0);
    (void) floor(1.0);
    (void) ceil(1.0);
    (void) rint(1.0);
}
