/*
 *  SORT:  shell sort of an array of pointers to 'real's'  
 *              See K&R pp108
 *     20-jun-01	gcc3 			pjt
 */

#include <stdinc.h>
        
void sort (real *x[] , int n)
{
    int    gap, i, j;
    real  *temp;
#ifdef TESTBED
    int k;
#endif  
    for (gap=n/2; gap>0; gap /= 2)
        for (i=gap; i<n; i++)
            for (j=i-gap; j>=0; j -= gap) {
                if (*x[j] <= *x[j+gap])
                    break;
                temp = x[j];
                x[j] = x[j+gap];
                x[j+gap] = temp;
#ifdef TESTBED
                printf ("xchg @ %d %d %d | ",gap,i,j);
                for (k=0; k<10; k++)
                   printf ("%2.0f ",*x[k]);
                printf ("\n");
#endif
            }
}


#ifdef TESTBED

double xxx[10]={9, 5, 4, 7, 6, 2, 1, 3, 8, 0};

main()
{
    real *x[10], *xold[10];
    int i;
        
    for (i=0; i<10; i++) 
        x[i] = xold[i] = &xxx[i];

    sort (x, 10);
        
    for (i=0; i<10; i++) 
        printf ("&%d=> %f  & %d=> %f  \n",xold[i],*xold[i],x[i],*x[i]);
}
#endif

                
        
