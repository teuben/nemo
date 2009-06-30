/*
 *  in situ data conversion routines:
 *      convert_d2f     double to float, loss of precision of course
 *      convert_f2d     float to double, random info added of course
 *
 *  work done IN SITU; either starting from the bottom upwards, or
 *  top downwards.
 *
 *
 *     25-may-91  written for some new code?     PJT
 *     25-feb-92  amazing, had to make gcc2.0 happy PJT 
 *     19-aug-92  added illegal address protection, <nemoinc>
 *                added convert_f2d but never tested...      PJT
 *     20-jun-01  gcc3
 */

#include <stdinc.h>

int convert_d2f(int n,double *from, float *to)
{
    if (from==NULL) error("convert_d2f: illegal from=NULL address");
    if (to==NULL)   error("convert_d2f: illegal to=NULL address");
    if (n<1) return 0;

    while(n--)                     /* can be done in situ */
	*to++ = *from++;
    return 1;
}

int convert_f2d(int n,float *from,double *to)
{
    if (from==NULL) error("convert_f2d: illegal from=NULL address");
    if (to==NULL)   error("convert_f2d: illegal to=NULL address");
    if (n<1) return 0;

    from = from+n-1;               /* find the address at top */
    to = to+n-1;

    while(n--)                     /* can be done in situ */
	*to-- = *from--;
    return 1;
}

