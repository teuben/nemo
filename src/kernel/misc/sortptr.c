/*
 *  SORTPTR:  shell sort of an array of doubles into a int index-array
 *              See K&R pp108
 *
 *      input:  x[]   array of 'double' values
 *              n     number of elements in x to sort
 *
 *      output: idx[] 'pointer' array, such that x[idx[i-1]]<x[idx[i]] for
 *                    i=1..n
 */
        
#include <stdinc.h>

void sortptr (real *x ,int *idx, int n)
{
    int    gap, i, j;
    int    tmp;
#ifdef TESTBED
    int k;
#endif          
    for (i=0; i<n; i++)
        idx[i]=i;               /*  one-to-one */
                
    for (gap=n/2; gap>0; gap /= 2)
        for (i=gap; i<n; i++)
            for (j=i-gap; j>=0; j -= gap) {
                if (x[idx[j]] <= x[idx[j+gap]])
                    break;          /* in order */
                tmp = idx[j];
                idx[j]    = idx[j+gap];
                idx[j+gap]= tmp;
#ifdef TESTBED
                printf ("xchg @ %d %d %d | ",gap,i,j);
                for (k=0; k<10; k++)
                   printf ("%2.0f ",x[idx[k]]);
                printf ("\n");
#endif
            }
}


#ifdef TESTBED

static real   x[10]={9, 5, 4, 7, 6, 2, 1, 3, 8, 0};
static int    idx[10];

main()
{
    int i;
        
    sortptr (x, idx, 10);
        
    for (i=0; i<10; i++) 
        printf ("%d=> %f  & %d sorted %f  \n",
                 i,   x[i], idx[i],   x[idx[i]]);
}
#endif

