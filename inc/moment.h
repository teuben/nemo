/*
 * moment.h	data stucture to aid in moment & minmax calculations
 */


typedef struct moment {
    int n;		    /* number of data accumulated */
    int mom;		    /* highest moment (use -1 if minmax is all) */
    real *sum;		    /* mom+1 length array with moment^{0..mom} sums */
    real datamin, datamax;  /* min & max of data */
} Moment, *MomentPtr; 

void ini_moment   (Moment *, int);		/* allocates */
void accum_moment (Moment *, real, real);	/* accumulates */
void decr_moment  (Moment *, real, real);	/* decrements (dangerous) */
void reset_moment (Moment *);       	/* resets */

real show_moment  (Moment *, int);     /* general case to peek at (special) values */

int  n_moment     (Moment *);	/* number of moments */
real sum_moment   (Moment *);	/* computes sum0 */
real mean_moment  (Moment *);	/* computes mean (mom=-1) */
real sigma_moment (Moment *);	/* computes dispersion around mean (mom=-2) */
real skewness_moment (Moment *);	/* computes special 3rd moment (mom=-3) */
real kurtosis_moment (Moment *);	/* computes special 4th moment (mom=-4) */

real min_moment (Moment *);
real max_moment (Moment *);
