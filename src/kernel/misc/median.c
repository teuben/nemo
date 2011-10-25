/*
 * Median: find the median of an array
 *
 *	20-jun-01	gcc3
 */

#include <stdinc.h>

static int *ix=NULL;
static int nix=0;

static int last_n = 0;
static real *last_x = NULL;

extern     void sortptr(real *x, int *ix, int n);

real median(int n, real *x)
{
    if (n <= 0) error("median: n=%d",n);
    if (x == 0) error("median: x=NULL");

    last_n = n;
    last_x = x;

    if (nix==0) {                       /* first time allocate */
        nix = n;
        ix = (int *) allocate(sizeof(real)*nix);
    } else if (nix < n) {               /* re-allocation */
        nix = n;
        ix = (int *) reallocate(ix, sizeof(real)*nix);
    }
    sortptr(x,ix,n);
    if (n % 2)
        return x[ix[(n-1)/2]];
    else
        return  0.5*(x[ix[n/2]] + x[ix[n/2-1]]);
}

/* these quartiles depend on a previous call to median, we
   are reusing the ix[] pointing array
   this code will thus not run in MP mode
*/

real median_q1(int n, real *x)
{
  int n1;
  if (last_n==0) error("No previous median()");
  if (last_n != n || last_x != x) error("Bad previous median() call");
  if (n>4) {
    n1=(n+1)/4;
    return(x[ix[n1]]);
  } else
    error("median_q1: too few points");

}
real median_q3(int n, real *x)
{
  int n3;
  if (last_n==0) error("No previous median()");
  if (last_n != n || last_x != x) error("Bad previous median() call");
  if (n>4) {
    n3=((n+1)*3)/4;
    return(x[ix[n3]]);
  } else
    error("median_q1: too few points");

}


 
void init_median(int size)
{
    nix = size;
    ix = (int *) allocate(sizeof(real)*nix);
}

void finis_median(void)
{
    free(ix);
    ix = 0;
    nix = 0;
    last_n = 0;
    last_x = NULL;
}
