/*
 * SQR, QBE: compute x*x and x*x*x.
 * DEX: compute inverse log^10: 10**x
 *	dark ages	created			JEB
 *	22-jan-95	ansi proto		PJT
 */

#include <math.h>

double sqr(double x)
{
    return x*x;
}

double qbe(double x)
{
    return x*x*x;
}

double dex(double x)
{
    /* extern double pow(double,double); */
    return pow(10.0, x);
}
