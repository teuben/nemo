/* NEMOfied  - added local shell sort */

#include <nemo.h>

extern real probks(real alam);

void my_sort(real *arr, int n);		/* nemo's sort, not NR */

void kstwo(real data1[], unsigned long n1, real data2[], unsigned long n2,
        real *d, real *prob)
{
        unsigned long j1=0,j2=0;
        float d1,d2,dt,en1,en2,en,fn1=0.0,fn2=0.0;

        my_sort(data1,n1);
        my_sort(data2,n2);
        en1=n1;
        en2=n2;
        *d=0.0;
        while (j1 < n1 && j2 < n2) {
                if ((d1=data1[j1]) <= (d2=data2[j2])) fn1=(++j1)/en1;
                if (d2 <= d1) fn2=(++j2)/en2;
                if ((dt=fabs(fn2-fn1)) > *d) *d=dt;
        }
        en=sqrt(en1*en2/(en1+en2));
        *prob=probks((en+0.12+0.11/en)*(*d));
}
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

/*
 *  SORT:  shell sort of an array of pointers to 'real's'  
 *              See K&R pp108
 */

#include <stdinc.h>
        
void my_sort (real *x , int n)
{
    int    gap, i, j;
    real   temp;

    for (gap=n/2; gap>0; gap /= 2)
        for (i=gap; i<n; i++)
            for (j=i-gap; j>=0; j -= gap) {
                if (x[j] <= x[j+gap])
                    break;
                temp = x[j];
                x[j] = x[j+gap];
                x[j+gap] = temp;
            }
}


