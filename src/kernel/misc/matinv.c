/*
     MATINV is a routine for inverting a matrix. The algorithm used is
     the Gauss-Jordan algorithm described in Stoer, Numerische Matema-
     tik,1 Teil p     

     Author:  T. Oosterloo (sheltran)
              P. Teuben (c)

     Date:    09-feb-1982

     Updates: 20-apr-1983  calculation of determinant. (K. Begeman)
              22-apr-1983  routine implemented.
               5-may-1983  permutation error repaired (K. Begeman)
               6-aug-1983  no calculation of real determinat  to avoid
                           floating overflow (K. Begeman)
              13-feb-1987  limit extended to 512 (A. Leene)
              30-sep-1990  written in C for Nemo's lsq_ routines            PJT
               9-jul-1991  replaced '=-' with '= -' (appease IRIS compiler) PJT
	      25-feb-92    happy gcc2.0                                     PJT
              20-jun-01    gcc 3.0                                          pjt
*/

#include <stdinc.h>

/*	maximum size of matrix */
#define  LIMIT 32

void matinv(
    real *matrix,		/* a(n,m) us used in here */
    int size,		/* actual size of matrix */
    int sizemp,		/* declared size of matrix by caller */
    real *det)		/* returned determinant; 0=singular */
{
    int  i,j,row,k,evin,per[LIMIT];
    real max,even,hv[LIMIT],mjk;

    *det=1.0;    /* in case of normal end */
    
    for (i=0; i<size; i++)    /* set permutation array */
        per[i]=i;

    for (j=0; j<size; j++) { /* in j-th column, set row with largest element */
        max=ABS(matrix[j*sizemp+j]);
        row=j;
        i=j+1;
        while (i < size) {
            if (ABS(matrix[i+j*sizemp]) > max) {
                max=ABS(matrix[i+j*sizemp]);
                row=i;
            }
            i++;
        }
        if (matrix[row+j*sizemp] == 0.0)  {   /* determinant zero ? */
            *det=0.0;        /* no solution possible */
            return;
        }
       
        if (row > j) {  /* if largest element not on diagonal: */
                        /* then permute rows */
            for (k=0; k<size; k++) {  /*    permutation loop */
                even=matrix[j+k*sizemp];
                matrix[j+k*sizemp]=matrix[row+k*sizemp];
                matrix[row+k*sizemp]=even;
            }
            evin=per[j];        /* keep track of permutation */
            per[j]=per[row];
            per[row]=evin;
        }
        even=1.0/matrix[j+j*sizemp];       /* modify column */
        for(i=0; i<size; i++) 
            matrix[i+j*sizemp]=even*matrix[i+j*sizemp];

        matrix[j+j*sizemp]=even;
        k=0;            /* modify rest of matrix */
        while (k < j) {
            mjk=matrix[j+k*sizemp];
            i=0;
            while (i < j) {
                matrix[i+k*sizemp]=matrix[i+k*sizemp]-matrix[i+j*sizemp]*mjk;
                i++;
            }
            i=j+1;
            while (i < size) {
                matrix[i+k*sizemp]=matrix[i+k*sizemp]-matrix[i+j*sizemp]*mjk;
                i++;
            }
            matrix[j+k*sizemp] = -even*mjk;
            k++;
        }
        k=j+1;
        while (k < size) {
            mjk=matrix[j+k*sizemp];
            i=0;
            while (i < j) {
                matrix[i+k*sizemp]=matrix[i+k*sizemp]-matrix[i+j*sizemp]*mjk;
                i++;
            }
            i=j+1;
            while (i < size) {
                matrix[i+k*sizemp]=matrix[i+k*sizemp]-matrix[i+j*sizemp]*mjk;
                i++;
            }
            matrix[j+k*sizemp] = -even*mjk;
            k++;
        }
    }  /*    end of loop through columns */

    for (i=0; i<size; i++) {  /* finally, repermute columns */
        for(k=0; k<size; k++)
            hv[per[k]]= matrix[i+k*sizemp];
        for(k=0; k<size; k++)
            matrix[i+k*sizemp]=hv[k];
    }
}


#if defined(TESTBED)

real a[] = { 1.0, 2.0, 3.0, 4.0 };

main()
{
    real det;

    matinv(a,2,2,&det);
    printf("det = %f\n a^-1 = %f %f %f %f\n",
        det, a[0], a[1], a[2], a[3]);

}
        
#endif
