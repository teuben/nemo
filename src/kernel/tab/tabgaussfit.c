/*
 * TABGAUSSFIT: compare some gaussfitting routines, using various methods:
 *              nemo (as in tablsqfit), gsl, and xpnt 
 *
 *      23-nov-05  V1.0
 *      15-feb-2025  V1.1 - tried resurrecting
 *
 */

#include <stdinc.h>  
#include <getparam.h>
#include <lsq.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>

string defv[] = {
    "in=???\n           input (table) file name",
    "xcol=1,2\n         column(s) for x,y",
    "ycol=3\n           column(s) for z",
    "dycol=4\n          column for sigma-z, if used",
    "fit=gauss2d\n      fitmode (xpnt, nemo, gsl)",
    "nsigma=-1\n        delete points more than nsigma away?",
    "nmax=10000\n       Default max allocation",
    "tab=f\n            short one-line output?",
    "VERSION=1.1\n      15-feb-2025 PJT",
    NULL
};

string usage="testing 2d gaussian fitting";


/**************** SOME GLOBAL VARIABLES ************************/

#if !defined(HUGE)
#define HUGE 1e20
#endif

#ifndef MAXCOL
#define MAXCOL 16
#endif

typedef struct column {
    int maxdat;     /* allocated length of data */          /* not used */
    int ndat;       /* actual length of data */             /* not used */
    real *dat;      /* pointer to data */
    int colnr;      /* column number this data came from */ /* not used */
} a_column;

int nxcol, nycol, xcolnr[MAXCOL], ycolnr[MAXCOL], dycolnr; 
a_column            xcol[MAXCOL],   ycol[MAXCOL],   dycol;

real xrange[MAXCOL*2];      /* ??? */

string method;              /* fit method (line, ellipse, ....) */
stream instr;               /* input  file */


int    nmax;                /* allocated space */
int    npt;                 /* actual number of points from table */
real   nsigma;              /* fractional sigma removal */
real  a,b;                  /* fit parameters in: y=ax+b  */
bool Qtab;                  /* do table output ? */

struct data {
  size_t n;
  double *x;
  double *y;
  double *z;
  double *s;
};

void setparams(void);
void setrange(real *rval, string exp);
void read_data(void);
void print_state(size_t iter, gsl_multifit_fdfsolver *s);
void do_xpnt(void);
void do_gsl(void);
void do_nemo(void);
void do_peak(void);


/* XPNT's fortran code, slightly cleaned up for modern g77 */
extern void mqgaus2d_(float *X,float *Y,float *Z,float *S,int *N,float *XC,float *YC,
		      float *PEAK,float *FWHM,float *SDEV,float *CHISQ,float *BADFIT);


int gauss2d_f(const gsl_vector *x, void *params, gsl_vector *f)
{
  size_t n = ((struct data *)params)->n;
  double *x1 = ((struct data *)params)->x;
  double *x2 = ((struct data *)params)->y;
  double *y = ((struct data *)params)->z;
  double *s = ((struct data *)params)->s;

  double A  = gsl_vector_get (x, 0);
  double x0 = gsl_vector_get (x, 1);
  double y0 = gsl_vector_get (x, 2);
  double b  = gsl_vector_get (x, 3);

  size_t i;

  for (i=0; i<n; i++) {
    double Yi = A * exp(    -((x1[i]-x0)*(x1[i]-x0)  +  (x2[i]-y0)*(x2[i]-y0))/(2*b*b) );
    gsl_vector_set(f,i,(Yi-y[i])/s[i]);
  }
  return GSL_SUCCESS;
}

int gauss2d_df(const gsl_vector *x, void *params, gsl_matrix *J)
{
  size_t n = ((struct data *)params)->n;
  double *x1 = ((struct data *)params)->x;
  double *x2 = ((struct data *)params)->y;
  double *y = ((struct data *)params)->z;
  double *s = ((struct data *) params)->s;
  
  double A  = gsl_vector_get (x, 0);
  double x0 = gsl_vector_get (x, 1);
  double y0 = gsl_vector_get (x, 2);
  double b  = gsl_vector_get (x, 3);

  size_t i;

  for (i=0; i<n; i++) {
    double x11=(x1[i]-x0)*(x1[i]-x0);
    double x22=(x2[i]-y0)*(x2[i]-y0);
    double bb2=2*b*b;
    double e  = exp(-(x11+x22)/bb2);
    double sigma=s[i];

    gsl_matrix_set (J, i, 0, e/sigma);
    gsl_matrix_set (J, i, 1, A*e*x11/(bb2*sigma));
    gsl_matrix_set (J, i, 2, A*e*x22/(bb2*sigma));
    gsl_matrix_set (J, i, 3, A*e*(x11+x11)*2/(bb2*sigma*b));
  }
  return GSL_SUCCESS;
}

int gauss2d_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J)
{
  gauss2d_f(x,params,f);
  gauss2d_df(x,params,J);
  return GSL_SUCCESS;
}




/****************************** START OF PROGRAM **********************/

void nemo_main()
{

    setparams();
    read_data();

    if (scanopt(method,"xpnt")) {
        do_xpnt();
    } else if (scanopt(method,"gsl")) {
    	do_gsl();
    } else if (scanopt(method,"nemo")) {
    	do_nemo();
    } else
      error("fit=%s invalid", getparam("fit"));
}

void setparams()
{
    string inname = getparam("in");
    nmax = nemo_file_lines(inname,getiparam("nmax"));
    if (nmax<0) error("Error opening %s",inname);
    if (nmax==0) error("No data?");
    instr = stropen (inname,"r");

    nxcol = nemoinpi(getparam("xcol"),xcolnr,MAXCOL);
    if (nxcol<0) error("Illegal xcol= nxcol=%d",nxcol);
    nycol = nemoinpi(getparam("ycol"),ycolnr,MAXCOL);
    if (nycol<0) error("Illegal ycol= nycol=%d",nycol);
    if (hasvalue("dycol"))
        dycolnr = getiparam("dycol");
    else
        dycolnr = 0;

    xrange[0] = -HUGE;
    xrange[1] = HUGE;
    
    method = getparam("fit");
    nsigma = getdparam("nsigma");
    Qtab = getbparam("tab");
}

void setrange(real *rval, string rexp)
{
    char *cptr;

    cptr = strchr(rexp, ':');
    if (cptr) {
        rval[0] = natof(rexp);
        rval[1] = natof(cptr+1);
    } else {
        rval[0] = 0.0;
        rval[1] = natof(rexp);
    	warning("Range taken from 0 - %g",rval[1]);
    }
}

void read_data()
{
    real *coldat[2*MAXCOL+1];
    int colnr[2*MAXCOL+1], ncols = 0, i, j;

    for (i=0; i<nxcol; i++) {
        coldat[ncols] = xcol[i].dat = (real *) allocate(nmax * sizeof(real));
        colnr[ncols] = xcolnr[i];        
        ncols++;
    }
    for (i=0; i<nycol; i++) {
        coldat[ncols] = ycol[i].dat = (real *) allocate(nmax * sizeof(real));
        colnr[ncols] = ycolnr[i];        
        ncols++;
    }
    if (dycolnr>0) {
        coldat[ncols] = dycol.dat = (real *) allocate(nmax * sizeof(real));
        colnr[ncols] = dycolnr;
        ncols++;
    }
    
    npt = get_atable(instr,ncols,colnr,coldat,nmax);
    if (npt < 0) {
        npt = -npt;
       	warning("Could only read %d data",npt);
    }


    /* special case for nxcol=1  ... what to do for nxcol > 1 ??? */
    /* should also handle nycol > 1  but does not yet             */

    if (nxcol == 1 && nycol == 1) {
        for(i=0, j=0; i<npt; i++) {
          if(xrange[0] <= xcol[0].dat[i] && xcol[0].dat[i] <= xrange[1]) {    /* sub-select on X */
              xcol[0].dat[j] = xcol[0].dat[i];
              ycol[0].dat[j] = ycol[0].dat[i];
              if (dycolnr>0) dycol.dat[j] = dycol.dat[i];
              j++;
           }
        }
        dprintf(1,"Copied over %d/%d data within xrange's\n",j,npt);
	npt = j;
    }
       
    if (npt==0) error("No data");
}




void
print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f % 15.8f "
          "|f(x)| = %g\n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2), 
          gsl_vector_get (s->x, 3), 
          gsl_blas_dnrm2 (s->f));
}



void do_xpnt()
{
  float *x, *y, *z, *s, peak, fwhm, xc, yc, chisq, sdev[4], badfit;
  int n;
  int i, j, gorder=3, neg=0;

  if (nycol<1) error("Need 1 value for ycol=");
  if (nxcol<2) error("Need 2 values for xcol=");
  if (npt<4) error("Need at least 4 data points for gauss2d fit");

  x = (float *) allocate(npt*sizeof(float));
  y = (float *) allocate(npt*sizeof(float));
  z = (float *) allocate(npt*sizeof(float));
  s = (float *) allocate(npt*sizeof(float));

  peak = -1.0;
  for (i=0; i<npt; i++) {
    x[i] = xcol[0].dat[i];
    y[i] = xcol[1].dat[i];
    z[i] = ycol[0].dat[i];
    if (z[i] > peak) peak=z[i];
#if 0
    s[i] = dycol.dat[i];
#else
    s[i] = 0.1;
#endif
    dprintf(1,"%g %g %g %g\n",x[i],y[i],z[i],s[i]);
  }
  fwhm = 1.0;
  xc = yc = 1.0;
  if (0) {
    fwhm *= -1.0;
    peak *= -1.0;
  }
  dprintf(0,"initial guess: peak=%g  fwhm=%g\n",peak,fwhm);
         /* X,Y,Z,S ,N  ,XC  ,YC  ,PEAK  ,FWHM ,SDEV  ,CHISQ  ,BADFIT */
  // need to resurrect this
  mqgaus2d_(x,y,z,s,&n, &xc, &yc, &peak, &fwhm, sdev, &chisq, &badfit); 

  printf("mggauss:\n");
  printf("  center: %g %g   (%g %g)\n",xc,yc,sdev[0],sdev[1]);
  printf("  peak:   %g  (%g)\n",peak,sdev[2]);
  printf("  fwhm:   %g  (%g)\n",fwhm,sdev[3]);
  printf("  chisq:  %g\n",chisq);
  printf("  badfit: %g\n",badfit);

}

void do_gsl()
{
  size_t i, iter=0;
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  const size_t p = 4;
  const size_t n = npt;
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  struct data d;
  gsl_multifit_function_fdf f;
  double x_init[4] = {3.1, 1.1, 0.9, 1.1};
  gsl_vector_view x = gsl_vector_view_array (x_init, p);

  d.n = n;
  d.x = xcol[0].dat;
  d.y = xcol[1].dat;
  d.z = ycol[0].dat;
  d.s = dycol.dat;

  f.f = &gauss2d_f;
  f.df = &gauss2d_df;
  f.fdf = &gauss2d_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);

  print_state (iter, s);
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate (s);
    printf ("status = %s\n", gsl_strerror (status));
    print_state (iter, s);
    if (status)
      break;

    status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
  } while (status == GSL_CONTINUE && iter < 500);

  //gsl_multifit_covar (s->J, 0.0, covar);     // no member?


#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  printf ("A      = %.5f +/- %.5f\n", FIT(0), ERR(0));
  printf ("x0     = %.5f +/- %.5f\n", FIT(1), ERR(1));
  printf ("y0     = %.5f +/- %.5f\n", FIT(2), ERR(2));
  printf ("b      = %.5f +/- %.5f\n", FIT(3), ERR(3));

  { 
    double chi = gsl_blas_dnrm2(s->f);
    printf("chisq/dof = %g\n",  pow(chi, 2.0)/ (n - p));
  }

  printf ("status = %s\n", gsl_strerror (status));
  
  gsl_multifit_fdfsolver_free (s);
}


/*
 *   GAUSS2D:   y = A * exp( -[(x-x0)^2 + (y-y0)^2]/2b^2 )
 *        needs min. 4 points for a linear fit
 *            
 */

void do_nemo()
{
  real mat[(MAXCOL+1)*(MAXCOL+1)], vec[MAXCOL+1], sol[MAXCOL+1], a[MAXCOL+2], sum;
  int i, j, gorder=3, neg=0;
  real A, b, x0, y0, x,y,z;

  if (nycol<1) error("Need 1 value for ycol=");
  if (nxcol<2) error("Need 2 values for xcol=");
  if (npt<4) error("Need at least 4 data points for gauss2d fit");

  lsq_zero(gorder+1, mat, vec);
  for (i=0; i<npt; i++) {
    a[0] = 1.0;                     /* ln A - (x0^2+y0^2)/2b^2 */
    a[1] = xcol[0].dat[i];          /* x0/b^2  */
    a[2] = xcol[1].dat[i];          /* y0/b^2  */
    a[3] = sqr(a[1]) + sqr(a[2]);   /* -1/2b^2 */
    if (ycol[0].dat[i] <= 0) {
      neg++;
      continue;
    }
    a[4] = log(ycol[0].dat[i]);
    lsq_accum(gorder+1,mat,vec,a,1.0);
  }
  if (neg > 0) {
    warning("Ignored %d negative data",neg);
    if (npt-neg < 4) error("Too many points rejected for a gauss2d fit");
  }
  lsq_solve(gorder+1,mat,vec,sol);
  printf("gauss2d fit:\n");
  for (j=0; j<=gorder; j++) printf("%g ",sol[j]);
  printf("\n\n");
  printf("  y = A * exp( -[(x-x0)^2 + (y-y0)^2]/2b^2 ):\n\n");
  if (sol[3] > 0) {
    warning("Bad gauss2d fit: 1/b^2 = %g\n",sol[3]);
    return;
  }

  x0 = -sol[1]/(2*sol[3]);
  y0 = -sol[2]/(2*sol[3]);
  b = sqrt(-1/(2*sol[3]));
  A = exp(sol[0] - sol[3]*(x0*x0+y0*y0));
  printf("     A  = %g\n",A);
  printf("     x0 = %g\n",x0);
  printf("     y0 = %g\n",y0);
  printf("     b  = %g\n",b);
  for (i=0; i<npt; i++) {
    x = xcol[0].dat[i];
    y = xcol[1].dat[i];
    z = sol[0] + sol[1]*x + sol[2]*y + sol[3]*(x*x+y*y);
    if (ycol[0].dat[i] <= 0) continue;
    dprintf(1,"%g %g %g %g => %g\n",x,y,log(ycol[0].dat[i]),z,log(ycol[0].dat[i])-z);
  }
}




/* fit a peak, 
 * for now this is the my_poly() code, though forced with order=2 
 * first find the peak, then take the two points on either side
 * to find an exact solution
 */

void do_peak()
{
    real mat[(MAXCOL+1)*(MAXCOL+1)], vec[MAXCOL+1], sol[MAXCOL+1], a[MAXCOL+2];
    real *x, *y, ymax;
    int i, j, k, range=1, order;
    
    order = 2;

    if (nycol<1) error("Need 1 value for ycol=");

    x = xcol[0].dat;
    y = ycol[0].dat;
    for (i=1, j=0; i<npt; i++)			/* find the peak, at j */
        if (y[j] < y[i]) j = i;
    if (j==0 || j==npt-1) {			/* handle edge cases */
        warning("Peak at the edge");
        printf("%g %g\n",x[j],y[j]);
        return;
    }
    if (range==2) {
    	if (j==1 || j==npt-2) {
    	    warning("Peak too close to edge");
    	    return;
    	}
    }

    lsq_zero(order+1, mat, vec);
    for (i=j-range; i<=j+range; i++) {
        a[0] = 1.0;
        for (k=0; k<order; k++) {
            a[k+1] = a[k] * (x[i]-x[j]);
        }
        a[order+1] = y[i];
        lsq_accum(order+1,mat,vec,a,1.0);
    }
    lsq_solve(order+1,mat,vec,sol);
    dprintf(1,"Poly2 fit near j=%d (%g,%g)\n",j+1,x[j],y[j]);
    printf("Peak:x,y= %g %g\n",
            x[j] - sol[1]/(2*sol[2]),
            sol[0]+sqr(sol[1])/(2*sol[2]));
}

