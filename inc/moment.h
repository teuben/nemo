/*
 * moment.h	data stucture to aid in moment & minmax calculations
 */


typedef struct moment {
    int mom;		    /* highest moment (use -1 if minmax is all) */
    int n;		    /* number of data accumulated so far */
    int ndat;               /* max number of data for moving moments */
    int idat;               /* index to last written data for moving moments */
    real *dat;              /* data[ndata] if moving moments used */
    real *wgt;              /* weight[ndata] if moving moments used */
    real *sum;		    /* mom+1 length array with moment^{0..mom} sums */
    bool *msk;              /* optional mask (@todo) */
    real datamin, datamax;  /* min & max of data */
    real sumn, sump;        /* separate sum of negative and positive numbers */
} Moment, *MomentPtr; 

void ini_moment   (Moment *, int, int);		/* allocates */
void accum_moment (Moment *, real, real);	/* accumulates */
void decr_moment  (Moment *, real, real);	/* decrements (dangerous) */
void reset_moment (Moment *);       	        /* resets */
void free_moment  (Moment *);                   /* frees allocs from ini_ */

real show_moment  (Moment *, int);     /* general case to peek at (special) values */

int  n_moment     (Moment *);	/* number of moments */
real sum_moment   (Moment *);	/* computes sum0 */
real mean_moment  (Moment *);	/* computes mean (mom=-1) */
real median_moment(Moment *);   /* only works if ndat > 0 */
real sigma_moment (Moment *);	/* computes weighted dispersion around mean (mom=-2) */
real rms_moment   (Moment *);	/* computes rms */
real mad_moment   (Moment *);   /* MAD  = median absolute deviation */
real mard_moment  (Moment *);	/* MARD = mean absolute relative difference */
real sratio_moment(Moment *);   /* sratio = (sump+sumn)/(sump-sumn) */

real skewness_moment (Moment *);	/* computes special 3rd moment (mom=-3) */
real kurtosis_moment (Moment *);	/* computes special 4th moment (mom=-4) */
real h3_moment (Moment *);              /* see kinemetry */
real h4_moment (Moment *);              /* see kinemetry */

real min_moment (Moment *);
real max_moment (Moment *);

void compute_robust_moment(Moment *);   /* call this first to grab the next 4 values */
int  n_robust_moment(Moment *);
real min_robust_moment(Moment *);
real max_robust_moment(Moment *);
real mean_robust_moment(Moment *);
real sigma_robust_moment(Moment *);
real median_robust_moment(Moment *);
void robust_range(Moment *, real *range);

