/*
 * MinMax: find Min and Max of an array.
 *
 */

#include <stdinc.h>

int minmax(int n, real *array, real *amin, real *amax)
{
    int i;
    real datamin, datamax, *a = array;

    if (n<1) return 0;
    
    datamin = *a;
    datamax = *a;

    a++;    /* point to first element */
    
    for (i=1; i<n; i++, a++)
        if (*a > datamax)
            datamax = *a;
        else if (*a < datamin)
            datamin = *a;

    *amin = datamin;
    *amax = datamax;
    return 0;
}
 
