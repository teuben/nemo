/*
 * Median: find the median of an array
 *
 *	20-jun-01	gcc3
 */

#include <stdinc.h>

static int *ix=NULL;
static int nix=0;

extern     void sortptr(real *x, int *ix, int n);

real median(int n, real *x)
{
    if (n <= 0) error("median: n=%d",n);
    if (x == 0) error("median: x=NULL");

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
}
