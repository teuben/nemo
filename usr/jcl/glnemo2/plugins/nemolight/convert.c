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
 *     11-dec-09  half-precision code added   PJT
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
      *to-- = *from--;             /* assumer user didn't cheat */
    return 1;
}

/* 
 * needs: ieeehalfprecision.c 
 * http://www.mathworks.com/matlabcentral/fileexchange/23173
 */

extern int singles2halfp(void *target, void *source, int numel);
extern int doubles2halfp(void *target, void *source, int numel);
extern int halfp2singles(void *target, void *source, int numel);
extern int halfp2doubles(void *target, void *source, int numel);

int convert_h2f(int n,halfp *from,float *to)
{
  return halfp2singles(to,from,n);
}

int convert_h2d(int n,halfp *from,double *to)
{
  return halfp2doubles(to,from,n);
}

int convert_d2h(int n,double *from,halfp *to)
{
  return doubles2halfp(to,from,n);  
}

int convert_f2h(int n,float *from,halfp *to)
{
  return singles2halfp(to,from,n);  
}
