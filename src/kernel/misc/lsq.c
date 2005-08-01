/*
 *  Linear Least Squares Fitting:  lsq_zero, lsq_accum, lsq_solve, lsq_cfill
 *	using normalized equations
 *
 *	29-sep-90	Created			Peter Teuben
 *	 8-dec-90	Isolated matinv.c	PJT
 *	15-nov-91	added lsq_cfill()	PJT
 *	18-feb-92       lsq_accum: treat w as weight, not 1/w*w	PJT
 *	25-feb-92	happy gcc2.0				PJT
 *	13-jun-94       some error() calls for obvious mistakes PJT
 *	22-jan-95	ansi prototypes				pjt
 *      16-feb-97       extern proto instead of nested		pjt
 */

#include <stdinc.h>

extern void matinv(real *, int, int,real *);

/* 
 * LSQ_ZERO:   reset LSQ matrix and r.h.s. vector to zero
 */

void lsq_zero(int n, real *mat, real *vec)
{
    int i;

    if (n<1) error("lsq_zero: n=%d",n);
    for (i=0; i<n*n; i++)
        mat[i] = 0.0;
    for (i=0; i<n; i++)
        vec[i] = 0.0;
}

/* 
 *  LSQ_ACCUM:   accumulate data in array a[] into the matrix mat[]
 *               The observed quantity is stored in the last
 *               slot [n] of a[], and is accumulated into vec[]
 *               Matrix mat[] is symmetric by definition.
 *               w is the weight factor.
 */

void lsq_accum(int n, real *mat, real *vec, real *a, real w)
	       /* DIM: mat[n*n]  vec[n]  a[n+1]  */
{
    int i, j;
  
    if (n<1) error("lsq_accum: n=%d",n);
    /* w = 1.0 / (w*w);		No, let's really call 'w' a weight */
    for (i=0; i<n; i++) {
        vec[i] += a[i]*a[n]*w;
        for (j=0; j<n; j++)
            mat[j*n+i] += a[i]*a[j]*w;
    }
}

/*
 *  LSQ_SOLVE:	solve 'sol' for  mat*sol=vec assuming lsq_accum filled
 *		the 'mat' and 'vec' properly.
 */

void lsq_solve(int n, real *mat, real *vec, real *sol)
{
    real det;
    int i, j;

    if (n<1) error("lsq_solve: n=%d",n);
    matinv(mat,n,n,&det);
    if (det == 0.0) {
        dprintf(1,"lsq_solve: singular matrix of order %d",n);
        return;
    }
    for (i=0; i<n; i++) {
        sol[i] = 0.0;
        for (j=0; j<n; j++)
            sol[i] +=  mat[i+j*n] * vec[j];
    }
}

/*
 * LSQ_CFILL:  Fill a particul column in the matrix
 */

void lsq_cfill(int n, real * mat, int c, real *vec)
{
    int i, off = c*n;

    if (n<1) error("lsq_cfill: n=%d",n);
    for (i=0; i<n; i++)
       mat[off+i] = vec[i];
}

#if defined(TESTBED)

#include <getparam.h>

string defv[] = {
        "VERSION=0",
	NULL,
};

real mat[2*2];
real vec[2];
real sol[2];

nemo_main()
{
    mat[0] = 2;
    mat[1] = 1;
    mat[2] = 2;
    mat[3] = 2;
    vec[0] = 4;
    vec[1] = 1;
    lsq_solve(2,mat,vec,sol);
    printf("sol = %g %g\n",sol[0],sol[1]);
}

#endif
