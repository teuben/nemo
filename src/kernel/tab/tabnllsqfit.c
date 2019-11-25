/*
 * TABNLLSQFIT: a general (non-linear) fitting program for tabular data 
 *
 *      12-jul-02  1.0  derived from tablsqfit, but using nllsqfit() now
 *      17-jul-02  1.1  allow NR's mrqmin() to be used via an emulator (nr_nllsqfit)
 *                      added initial code for dynamic object loading functions 
 *      11-sep-02  1.1e error/warning if fit is bad, but can still write out residuals
 *      31-oct02   1.2  add "fourier m=1" mode as 'arm' for Rahul - problems fitting  as 
 *                      AMPHASE though
 *      21-nov-02  1.3  implemented nsigma for "line" and "arm"
 *      22-dec-02  1.3a fixing code to make it work for -DSINGLEPREC
 *      12-feb-03  1.3b add "fourier m=2' mode as arm3 for Rahul by Rahul
 *      14-feb-03  1.6  back to PJT style (version # got screwed up)
 *      18-feb-03  1.6a changed named of output variables in arm and arm3 
 *      21-mar-03  1.7  optional bootstrapping to check on errors
 *       4-apr-03  1.8  added dypow= keyword, and fixed bug in handling dycol=
 *      10-mar-04  1.8b added Lorentzian fitting, fixed setting lab= for loadable functions
 *      17-apr-04  1.9a added printing out the chi-squared or RMS or whatever it can do
 *       3-may-05  1.9b add x/x0+y/y0=1 variant for a linear fit 
 *      25-aug-05  1.10 poly2 for Rahul
 *      21-nov-05  2.0  gauss2d, fix error in gauss1d
 *      24-apr-08  2.1  psf (phase structure function)
 *       7-may-10  2.2  grow
 *      14-may-11  2.3  different options for bootstrap method
 *       9-dec-12  3.0  xrange= now allows separate sections   a:b,c:d,....
 *      10-oct-13  4.0  method=
 *      26-may-16  4.1  the fit=grow recoded
 *
 *  line       a+bx
 *  plane      p0+p1*x1+p2*x2+p3*x3+.....     up to 'order'   (a 2D plane in 3D has order=2)
 *  poly       p0+p1*x+p2*x^2+p3*x^3+.....    up to 'order'   (paraboloid has order=2)
 *  exp        p0+p1*exp(-(x-p2)/p3)  
 *  grow       p1*(1-exp(x/p3))               lyaponov?
 *  arm        p0+p1*cos(x)+p2*sin(x)         special version for rahul 
 *  arm3       p0+p1*cos(x)+p2*sin(x)+p3*cos(3*x)+p4*sin(3*x) 
 *  loren      (p1/PI) / ( (x-p0)^2 + p1^2 )
 *  gauss1d    p0+p1*exp(-(x-p2)^2/(2*p3^2))
 *  gauss2d    p0+p1*exp(-[(x-x0)^2 + (y-y0)^2]/2b^2)
 *  psf        p0*x**p1*sin(y)**p2+p3
 *  
 */ 

#include <stdinc.h>  
#include <getparam.h>
#include <loadobj.h>
#include <filefn.h>
#include <moment.h>

string defv[] = {
    "in=???\n           input (table) file name",
    "xcol=1\n           column(s) for x, the independant variable(s)",
    "ycol=2\n           column(s) for y, the dependant variable(s)",
    "dycol=\n           optional column for sigma-y (weight = (1/dy**2)**<dypow>)",
    "dypow=1\n          optional extra power to control the weights",
    "xrange=\n          in case restricted range is used (1D only)",
    "fit=line\n         fitmode (line, plane, poly, gauss, exp, loren)",
    "order=2\n		Order of plane/poly fit",
    "out=\n             output file for some fit modes",
    "nsigma=-1\n        delete points more than nsigma away?",
    "par=\n             initial estimates of parameters (p0,p1,p2,...)",
    "free=\n            free(1) or fixed(0) parameters? [1,1,1,....]",
    "load=\n            If used, uses this dynamic object function (full path)",
    "x=\n               X-values to test the function for, then exit",
    "nmax=10000\n       Default max allocation",
    "tol=\n             Tolerance for convergence of nllsqfit",
    "lab=\n             Mixing parameter for nllsqfit",
    "itmax=50\n         Maximum number of allowed nllsqfit iterations",
    "format=%g\n        Output format for fitted values and their errors",
    "bootstrap=0\n      Bootstrapping to estimate errors",
    "seed=0\n           Random seed initializer",
    "method=gipsy\n     method:   Gipsy(nllsqfit), Numrec(mrqfit), MINPACK(mpfit)"
    "VERSION=4.1a\n     7-jul-2019 PJT",
    NULL
};

string usage="a non-linear least square fitting program for tabular data";

string cvsid="$Id$";

/**************** SOME GLOBAL VARIABLES ************************/

#if !defined(HUGE)
#define HUGE 1e20
#endif

#define MAXCOL 10
#define MAXPAR 32
#define MAXSIG 10

typedef struct column {
    int maxdat;     /* allocated length of data */          /* not used */
    int ndat;       /* actual length of data */             /* not used */
    real *dat;      /* pointer to data */
    int colnr;      /* column number this data came from */ /* not used */
} a_column;

#define MAXR 16

typedef struct range {
  int nr;
  real *rmin;
  real *rmax;
} a_range;

int nxcol, nycol, xcolnr[MAXCOL], ycolnr[MAXCOL], dycolnr;
real dypow;
a_column            xcol[MAXCOL],   ycol[MAXCOL], dycol,  bcol;
a_range    xrange;

/* real xrange[MAXCOL*2];      /* ??? */

string fit_method;          /* fit method (gipsy, numrec, minpack) */
string fit_object;          /* fit method (line, poly, ....) */

stream instr, outstr;       /* input / output file */


int    nmax;                /* allocated space */
int    npt;                 /* actual number of points from table */
real   nsigma[MAXSIG];      /* fractional sigma removals */
int    msigma;              /* max number of iterations in nsigma removals */
real   tol = -1;            /* tolerance for M-L convergence, if used */ 
real   lab = -1;            /* mixing parameter for M-L convergence, if used  */
int    itmax;               /* max. number of iterations */

int order;                  /* order of fits that have an order (poly's/hyperplane's) */

int mask[MAXPAR];           /* 1=free  0=fixed parameter */
real par[MAXPAR];           /* initial parameters */
int npar; 

char fmt[256];
string format;

bool Qtab;                  /* do table output ? */

int  nboot;

typedef real (*my_proc1)(real *, real *, int);
typedef void (*my_proc2)(real *, real *, real *, int);
typedef int  (*my_proc3)(real *, int, real *, real *, real *, int, real *, real *, int *, 
			 int, real, int, real, my_proc1, my_proc2);


my_proc1 fitfunc;
my_proc2 fitderv;


extern int nr_nllsqfit(real *, int, real *, real *, real *, int, real *, real *, int *, 
		       int, real, int, real, my_proc1, my_proc2);
extern int mp_nllsqfit(real *, int, real *, real *, real *, int, real *, real *, int *, 
	 	       int, real, int, real, my_proc1, my_proc2);
extern int    nllsqfit(real *, int, real *, real *, real *, int, real *, real *, int *, 
		       int, real, int, real, my_proc1, my_proc2);

extern double  xrandom(double a, double b);

extern real gammq(real a, real x);    /* from nr */


int setrange(a_range *r, string rexp);
int  inrange(a_range *r, real rval);


my_proc3 my_nllsqfit;    /* set via numrec= to be the Gipsy or NumRec routine */



static real func_gauss1d(real *x, real *p, int np)
{
  real a,b,arg;
  a = p[2]-x[0];
  b = p[3];
  arg = a*a/(2*b*b);
  return p[0] + p[1] * exp(-arg);
}

static void derv_gauss1d(real *x, real *p, real *e, int np)
{
  real a,b,arg;
  a = p[2]-x[0];
  b = p[3];
  arg = a*a/(2*b*b);
  e[0] = 1.0;
  e[1] = exp(-arg);
  e[2] = -p[1]*e[1] *  a  /  (b*b);
  e[3] =  p[1]*e[1] * a*a / (b*b*b);
}

static real func_gauss2d(real *x, real *p, int np)
{
  real a,b,c,arg;
  a = p[2]-x[0];
  b = p[3]-x[1];
  c = p[4];
  arg = (a*a+b*b)/(2*c*c);
  return p[0] + p[1] * exp(-arg);
}

static void derv_gauss2d(real *x, real *p, real *e, int np)
{
  real a,b,c,arg;
  a = p[2]-x[0];
  b = p[3]-x[1];
  c = p[4];
  arg = (a*a+b*b)/(2*c*c);
  e[0] = 1.0;
  e[1] = exp(-arg);
  e[2] = -p[1]*e[1] * a / (c*c);
  e[3] = -p[1]*e[1] * b / (c*c);
  e[4] =  p[1]*e[1] * (a*a+b*b) / (c*c*c);
}


static real func_exp(real *x, real *p, int np)
{
  real a,b,arg;
  a = x[0]-p[2];
  b = p[3];
  arg = a/b;
  return p[0] + p[1] * exp(-arg);
}

static void derv_exp(real *x, real *p, real *e, int np)
{
  real a,b,arg;
  a = x[0]-p[2];
  b = p[3];
  arg = a/b;
  e[0] = 1.0;
  e[1] = exp(-arg);
  e[2] = p[1]*e[1] / b;
  e[3] = e[2] * arg;
}

static real func_grow(real *x, real *p, int np)
{
  real arg;
  arg = x[0]/p[1];
  return p[0] * (1-exp(-arg));
}

static void derv_grow(real *x, real *p, real *e, int np)
{
  real arg,arg1;
  arg = x[0]/p[1];
  arg1 = exp(-arg);
  e[0] = 1-arg1;
  e[1] = -p[0]*arg1*x[0]/(p[1]*p[1]);
}

#ifdef TESTMP

int mp_func_line(int npt, int npar, double *par, double *dev, double **der, void *pri)
{ /* 1D line, npar=2  should be */
  int i;
  double f, *x, *y;

  x = pri->x;
  y = pri->y;

  for (i=0; i<npt; i++) {
    f = p[0] + p[1]*x[i];     /* Linear fit function */
    dy[i] = y[i] - f;
  }


  if (derivs) {
    int j;
    
  }
}

#endif

static real func_line(real *x, real *p, int np)
{
  return p[0] + p[1]*x[0];
}

static void derv_line(real *x, real *p, real *e, int np)
{
  e[0] = 1.0;
  e[1] = x[0];
}

static real func_plane(real *x, real *p, int np)
{
  int i;
  real r = p[0];

  for (i=0; i<order; i++)  r += p[i+1]*x[i];
    
  return r;
}

static void derv_plane(real *x, real *p, real *e, int np)
{
  int i;

  e[0] = 1.0;
  for (i=0; i<order; i++)  e[i+1] = x[i];

}

static real func_poly(real *x, real *p, int np)
{
  real r = p[order];
  int i;

  for (i=order; i>0; i--)
    r = x[0]*r + p[i-1];
  return r;
}

static void derv_poly(real *x, real *p, real *e, int np)
{
  real r;
  int i;

  r = e[0] = 1.0;
  for (i=0; i<order; i++) {
    r *= x[0];
    e[i+1] = r;
  }
}

static real func_poly2(real *x, real *p, int np)
{
  real r = x[0] - p[0];

  if (np != 4) error("func_poly2: np=%d ?",np);

  /* hardcoded for order=2 */
  r = p[1] * (1 + p[2]*r + p[3]*r*r);
  return r;
}

static void derv_poly2(real *x, real *p, real *e, int np)
{
  real r = x[0] - p[0];

  if (np != 4) error("derv_poly2: np=%d ?",np);

  /* hardcoded for order=2 */
  e[0] = 0.0;
  e[1] = (1 + p[2]*r + p[3]*r*r);
  e[2] = p[1] * r;
  e[3] = p[1] * r * r;
}


/* testing for Rahul - oct 2002 */

#define DPR 57.29577951308232

/* don't use PHASEAMP, somehow it's not working .. */
/* #define PHASEAMP */

static real func_arm(real *x, real *p, int np)
{
#ifdef PHASEAMP
    return p[0] + p[1]*cos( (x[0]-p[2])/DPR);
#else
  real y = x[0]/DPR;

  return p[0] + p[1]*cos(y) + p[2]*sin(y);
#endif
}

static void derv_arm(real *x, real *p, real *e, int np)
{
#ifdef PHASEAMP
    e[0] = 1.0;
    e[1] = cos( (x[0]-p[2])/DPR);
    e[2] = p[1]*sin ((x[0]-p[2])/DPR) / DPR;
#else
  real y = x[0]/DPR;

  e[0] = 1.0;
  e[1] = cos(y);
  e[2] = sin(y);
#endif
}

static real func_loren(real *x, real *p, int np)
{
  return (p[1]/PI)/((x[0]-p[0])*(x[0]-p[0])+p[1]*p[1]);
}

static void derv_loren(real *x, real *p, real *e, int np)
{
  real tmp = (x[0]-p[0])*(x[0]-p[0]) + p[1]*p[1];

  tmp = 1.0/(tmp*tmp*PI);
  e[0] =  2*(x[0]-p[0])*p[1]*tmp;
  e[1] =   ((x[0]-p[0])*(x[0]-p[0])-p[1]*p[1])*tmp;
}

static real func_arm3(real *x, real *p, int np)
{
  real y = x[0]/DPR;

  return p[0] + p[1]*cos(y) + p[2]*sin(y) + p[3]*cos(3*y) + p[4]*sin(3*y);
}

static void derv_arm3(real *x, real *p, real *e, int np)
{
  real y = x[0]/DPR;

  e[0] = 1.0;
  e[1] = cos(y);
  e[2] = sin(y);
  e[3] = cos(3*y);
  e[4] = sin(3*y);
}

static real func_psf(real *x, real *p, int np)
{
  real x1 = pow(x[0],p[1]);
  real x2 = pow(sin(x[1]),p[2]);

  return p[0] * x1 * x2   +   p[3];
}

static void derv_psf(real *x, real *p, real *e, int np)
{
  real x1 = pow(x[0],p[1]);
  real x2 = pow(sin(x[1]),p[2]);

  e[0] =      x1*x2;
  e[1] = p[0]*x1*x2*log(x[0]);
  e[2] = p[0]*x1*x2*log(sin(x[1]));
  e[3] =                        1.0;
}




/****************************** START OF PROGRAM **********************/

nemo_main()
{

    setparams();
    read_data();

    if (hasvalue("load")) {
      load_function(getparam("load"),fit_object);
      if (hasvalue("x"))
	do_function_test(getparam("x"));
      else
	do_function(fit_object);
    } else if (scanopt(fit_object,"line")) {
        do_line();
    } else if (scanopt(fit_object,"plane")) {
    	do_plane();
    } else if (scanopt(fit_object,"poly")) {
    	do_poly();
    } else if (scanopt(fit_object,"poly2")) {
    	do_poly2();
    } else if (scanopt(fit_object,"gauss1d")) {
    	do_gauss1d();
    } else if (scanopt(fit_object,"gauss2d")) {
    	do_gauss2d();
    } else if (scanopt(fit_object,"exp")) {
    	do_exp();
    } else if (scanopt(fit_object,"grow")) {
    	do_grow();
    } else if (scanopt(fit_object,"arm")) {
    	do_arm();
    } else if (scanopt(fit_object,"loren")) {
    	do_loren();
    } else if (scanopt(fit_object,"arm3")) {
    	do_arm3();
    } else if (scanopt(fit_object,"psf")) {
    	do_psf();
    } else
        error("fit=%s invalid; try [line,plane,poly,poly2,gauss1d,gauss2d,exp,arm,loren,arm3]",
	      getparam("fit"));
}

setparams()
{
    int i;
    string inname = getparam("in");
  
    nmax = nemo_file_lines(inname,getiparam("nmax"));
    if (nmax<0) error("Error opening %s",inname);
    if (nmax==0) error("No data?");
    instr = stropen (inname,"r");

    if (hasvalue("out"))
        outstr=stropen(getparam("out"),"w");
    else
        outstr=NULL;

    nxcol = nemoinpi(getparam("xcol"),xcolnr,MAXCOL);
    if (nxcol<0) error("Illegal xcol= nxcol=%d",nxcol);
    nycol = nemoinpi(getparam("ycol"),ycolnr,MAXCOL);
    if (nycol<0) error("Illegal ycol= nycol=%d",nycol);
    if (hasvalue("dycol"))
        dycolnr = getiparam("dycol");
    else
        dycolnr = 0;
    dypow = getrparam("dypow");
    dypow *= -2.0;

    if (hasvalue("xrange"))
      setrange(&xrange,getparam("xrange"));
    else
      xrange.nr = 0;

    if (hasvalue("tol")) tol = getdparam("tol");
    if (hasvalue("lab")) lab = getdparam("lab");
    itmax = getiparam("itmax");
    
    fit_object = getparam("fit");
    msigma = nemoinpr(getparam("nsigma"),nsigma,MAXSIG);
    order = getiparam("order");
    if (order<0) error("order=%d of %s cannot be negative",order,fit_object);

    if (hasvalue("par")) {         /* get the initial estimates of parameters */
      npar = nemoinpr(getparam("par"),par,MAXPAR);
      if (npar < 0) error("bad par=");
    } else
      npar = 0;
    for (i=npar; i<MAXPAR; i++)
      par[i] = 0.0;

    if (hasvalue("free")) {        /* determine which parameters are free */
      int nfree;
      nfree = nemoinpi(getparam("free"),mask,MAXPAR);
      if (nfree < 0) error("bad free=");
      for (i=nfree; i<MAXPAR; i++)
	mask[i] = 1;
    } else {
      for (i=0; i<MAXPAR; i++)
	mask[i] = 1;
    }
    fit_method = getparam("method");
    switch (*fit_method) {
    case 'g':
    case 'G':
      my_nllsqfit = nllsqfit;
      break;
    case 'n':
    case 'N':
      my_nllsqfit = nr_nllsqfit;
      break;
    case 'm':
    case 'M':
      my_nllsqfit = mp_nllsqfit;
      break;
    default:
      error("method=%s not supported, try Gipsy, Numrec, MINPACK",fit_method);
    }
    format = getparam("format");
    nboot = getiparam("bootstrap");
    init_xrandom(getparam("seed"));
}

/*
 * parse    some kind of range=min1:max1,min2:max2,....
 */

int setrange(a_range *r, string rexp)
{
  char *cptr;
  string *bs;
  int i,nr;

  bs = burststring(rexp,", ");
  nr = xstrlen(bs, sizeof(string)) - 1;
  r->nr = nr;
  if (nr==0) return 0;
  r->rmin = (real *) allocate(nr*sizeof(real));
  r->rmax = (real *) allocate(nr*sizeof(real));

  for (i=0; i<nr; i++) {
    dprintf(0,"bs=%s\n",bs[i]);
    cptr = strchr(bs[i], ':');
    if (cptr) {
      r->rmin[i] = natof(bs[i]);
      r->rmax[i] = natof(cptr+1);
      dprintf(1,"xrange(%d)= %g %g\n",i+1,r->rmin[i],r->rmax[i]);
    } else 
      error("xrange needs a comma separates series of xmin:xmax");
  }
  freestrings(bs);
  return nr;
}


int inrange(a_range *r, real rval)
{
  int i,  nr = r->nr;
  if (nr==0) return 1;
  for (i=0; i<nr; i++) {
    dprintf(1,"%g <? %g <? %g\n",r->rmin[i],rval,r->rmax[i]);
    if (r->rmin[i] <= rval && rval <= r->rmax[i]) return 1;
  }
  dprintf(1,"%g not in range\n",rval);
  return 0;
}


real data_rms(int n, real *d, real *dy, int m)
{
  int i;
  real sum = 0.0;

  if (dy)
    for (i=0; i<n; i++) {
      dprintf(1,"DEBUG(i,d,w) %d %g %g\n",i,d[i],dy[i]);
      sum += sqr(d[i])*dy[i];    /* dy was converted to 1/sigma^2 */
    }
  else {
    for (i=0; i<n; i++)
      sum += sqr(d[i]);
    printf("rms= %g\n",sqrt(sum/n));
  }

  printf("rms2/chi2= %g\n",sum);
#if 0
  return sqrt(sum);
#else
  if (n-m <= 0) warning("no DOF to compute rms/chi");
  return gammq(0.5*(n-m), 0.5*sum);
#endif
    
}


read_data()
{
    real *coldat[2*MAXCOL+1];
    int colnr[2*MAXCOL+1], ncols = 0, i, j, ncount;

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
    if (nboot>0) {
      bcol.dat = (real *) allocate(nmax * sizeof(real));
    }
    
    npt = get_atable(instr,ncols,colnr,coldat,nmax);
    if (npt < 0) {
        npt = -npt;
       	warning("Could only read %d data",npt);
    }
    if (dycolnr>0) {         /* convert dy such that it's a weight now */
      ncount = 0;
      for(i=0; i<npt; i++) {
	if (dycol.dat[i] < 0) 
	  ncount++;
	else
	  dycol.dat[i] = pow(dycol.dat[i], dypow);
      }
    }


    /* special case for nxcol=1  ... what to do for nxcol > 1 ??? */
    /* should also handle nycol > 1  but does not yet             */

    if (nxcol == 1 && nycol == 1) {
        for(i=0, j=0; i<npt; i++) {
	  if(inrange(&xrange,xcol[0].dat[i])) {      /* sub-select on X */
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

    for (i=0; i<npt; i++)
      dprintf(2,"DATA(x[0],y,wt): %g %g %g\n",xcol[0].dat[i],ycol[0].dat[i], dycolnr>0 ? dycol.dat[i] : 1.0);
}


load_function(string fname,string method)
{
  char func_name[80], derv_name[80];
  string path;

  mysymbols(getargv0());
  path = pathfind(".",fname);
  if (path == NULL) error("Cannot open %s",fname);

  if (method) {
    sprintf(func_name,"func_%s",method);
    sprintf(derv_name,"derv_%s",method);
  } else {
    sprintf(func_name,"func_loadobj");
    sprintf(derv_name,"derv_loadobj");
  }
  dprintf(0,"load_function: %s with %s\n",path,func_name);
  loadobj(path);
  
  fitfunc = (my_proc1) findfn(func_name);
  fitderv = (my_proc2) findfn(derv_name);
  if (fitfunc==NULL) error("Could not find func_loadobj in %s",fname);
  if (fitderv==NULL) error("Could not find derv_loadobj in %s",fname);
}

do_function(string method)
{
  real *x, *y, *dy, *d;
  int i,j,k,nrt, mpar[MAXPAR];
  real fpar[MAXPAR], epar[MAXPAR];
  int lpar = npar;        /* MUST be set */

  if (npar == 0) error("You must specify initial conditions for all parameters");
  sprintf(fmt,"p%%d= %s %s\n",format,format);
    
  if (nxcol < 1) error("nxcol=%d",nxcol);
  if (nycol < 1) error("nycol=%d",nycol);
  if (tol < 0) tol = 0.0;
  if (lab < 0) lab = 0.01;

  x = xcol[0].dat;
  y = ycol[0].dat;
  dy = (dycolnr>0 ? dycol.dat : NULL);
  d = (real *) allocate(npt * sizeof(real));
  
  for (i=0; i<lpar; i++) {
    mpar[i] = mask[i];
    fpar[i] = par[i];
  }

  nrt = (*my_nllsqfit)(x,1,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,itmax,lab, fitfunc,fitderv);
  printf("nrt=%d\n",nrt);
  printf("Fitting LOADED function \"%s\":  \n",method);
  for (k=0; k<lpar; k++)
    printf(fmt, k,fpar[k],epar[k]);
  if (nrt==-2)
    warning("No free parameters");
  else if (nrt<0)
    error("Bad fit, nrt=%d",nrt);

  if (outstr)
    for (i=0; i<npt; i++)
      fprintf(outstr,"%g %g %g\n",x[i],y[i],d[i]);

  printf("rms/chi = %g\n",data_rms(npt,d,dy,2));
}


do_function_test(string xvals)
{
  real *x, y, y0, dp, dyda[MAXPAR], lpar[MAXPAR];
  int i,j,k,l,n = getiparam("nmax");

  x = (real *) allocate(n * sizeof(real));
  n = nemoinpr(xvals,x,n);
  if (n < 0) error("Parsing x, try larger nmax=%d",getiparam("nmax"));
  for (i=0; i<n; i++) {
    y0 = y = (*fitfunc)(&x[i],par,npar);
    (*fitderv)(&x[i],par,dyda,npar);
    printf("%g %g  ",x[i],y);
    for (j=0; j<npar; j++)
      printf(" %g",dyda[j]);
    printf("\n");
  }
  if (n>1 && x[0]==x[1]) {
    printf("# Now computing numerical derivatives to check your functions (itermax=%d)\n",n);
    printf("# par iter Y dP  (Y-Y0)/dP   [(Y-Y0)/dP  - dY/dP]\n");
    for(j=0; j<npar; j++) {
      for(k=0;k<npar;k++) lpar[k] = par[k];
      y0 = (*fitfunc)(&x[0],lpar,npar);
      dp = 0.1;
      for(l=0;l<n; l++, dp /= 2.0) {
	lpar[j] = par[j] + dp; 
	y = (*fitfunc)(&x[0],lpar,npar);
	printf("%d %d %g %g %g [%g]\n",j,l,y,dp,(y-y0)/dp,(y-y0)/dp-dyda[j]);
      }
    }
  }
}



int remove_data(real *x, int nx, real *y, real *dy, real *d, int npt, real nsigma) 
{
  int i, j;
  real sigma, s;

  if (nsigma <= 0.0) return npt;           /* pass data through unchanged */

  for(i=0, sigma=0.0; i<npt; i++)
    sigma += sqr(d[i]);
  sigma /= (real) npt;
  sigma = nsigma*sqrt(sigma);
  for (i=0, j=0; i<npt; i++) {             /* loop over data, shift data in arrays */
    s = ABS(d[i]);
    if (s > sigma) continue;
    x[j] = x[i];
    y[j] = y[i];
    if (dy) dy[j] = dy[i];
    j++;
  }
  dprintf(0,"%d/%d points outside %g*sigma (%g)\n",npt-j,npt,nsigma,sigma);
  return j;
}

/*
 *  random_permute:  make a new random permutation
 *
 *  pick two random slots, and exchange them. repeat this N times
 *  
 */

#define random_permute random_permute1

random_permute1(int n, int *idx) 
{
  int i, j, k, tmp;
  double xn = n;

  for (i=0; i<n; i++) {
    j = (int) xrandom(0.0,xn);
    k = (int) xrandom(0.0,xn);
    if (j < 0 || j >= n) error("Error random_permute: %d %d",i,j);
    tmp = idx[j];
    idx[j] = idx[k];
    idx[k] = tmp;
  }
}

/*
 * iter:
 *   0     pick between 1...n-1 to permute  (n-1)
 *   1     pick between 2...n-1 to permute  (n-2)
 *
 *  n-3    pick between n-2..n-1 to permute (2)

# exp   nboot    mean    sigma
100     1000     3.08964 0.0035
1000    1000         962     30
                     987
1000    1000
         100         976    100
                     

*/

random_permute2(int n, int *idx) 
{
  int i, j, k, tmp;
  double xn = n;

  for (i=0; i<n-2; i++) {
    j = (int) xrandom(i+1.0,n+0.0);
    if (j <= i || j >= n) error("Error random_permute: %d %d",i,j);
    dprintf(1,"permuting (%d,%d)\n",i,j);
    tmp = idx[j];
    idx[j] = idx[i];
    idx[i] = tmp;
  }
  for (i=0; i<n; i++)
    dprintf(1,"%d ",idx[i]);
  dprintf(1,"\n");
}

/*
 * pick between 1..n, but with replacement, the official bootstrap way
 *
 */

random_permute3(int n, int *idx)
{
  int i;
  double xn = n;

  for (i=0; i<n; i++) {
    idx[i] = (int) xrandom(0.0,xn);  
    if (idx[i]<0 || idx[i]>=n) error("random_permute3 range error");
  }
}


#define bootstrap  bootstrap3

/* 
 * bootstrap1: take a number of new samples of the errors and distribute them 
 *             on the first fit. then refit and see what the distribution of
 *             the errors is.
 *             AKA resampling residuals, this incoorporates knowning a model
 */

void bootstrap1(int nboot, 
	       int npt, int ndim, real *x, real *y, real *dy, real *d, 
	       int npar, real *fpar, real *epar, int *mpar)
{
  real *y1, *d1, *bpar;
  int *perm, i, j, nrt;
  Moment *m;

  if (nboot < 1) return;

  perm = (int *) allocate(npt*sizeof(int));
  y1 = (real *) allocate(npt*sizeof(real));
  d1 = (real *) allocate(npt*sizeof(real));
  bpar = (real *) allocate(npar*sizeof(real));
  m = (Moment *) allocate(npar*sizeof(Moment));

  for (i=0; i<npt; i++)
    perm[i] = i;
  for (i=0; i<npar; i++) {
    bpar[i] = fpar[i];
    ini_moment(&m[i],2,0);
  }
  
  for (j=0; j<nboot; j++) {
    random_permute(npt,perm);
    for (i=0; i<npt; i++) {
      y1[i] = (*fitfunc)(&x[i],bpar,npar) + d[perm[i]];
    }
    nrt = (*my_nllsqfit)(x,ndim,y1,dy,d1,npt,fpar,epar,mpar,npar,tol,itmax,lab, fitfunc,fitderv);
    dprintf(1,"%g %g %g %g\n", fpar[0],fpar[1],epar[0],epar[1]);
    for (i=0; i<npar; i++) {
      accum_moment(&m[i],fpar[i],1.0);
    }
  }
  printf("bootstrap1= ");
  for (i=0; i<npar; i++)
    printf("%g %g ",mean_moment(&m[i]),sigma_moment(&m[i]));
  printf("\n");

  free(y1);
  free(d1);
  free(bpar);
  free(perm);
  free(m);
}
/* 
 * bootstrap3: take a new permution of the (x,y) sample, with 
 *             replacement and refit. 
 *             do this nboot time and see what the
 *             distribution of fitted values is
 *             this is a non-parametric way
 */

void bootstrap3(int nboot, 
	       int npt, int ndim, real *x, real *y, real *dy, real *d, 
	       int npar, real *fpar, real *epar, int *mpar)
{
  real *x1, *y1, *dy1, *d1, *bpar;
  int *perm, i, j, nrt;
  Moment *m;

  if (nboot < 1) return;

  perm = (int *) allocate(npt*sizeof(int));
  x1  = (real *) allocate(npt*sizeof(real));
  y1  = (real *) allocate(npt*sizeof(real));
  if (dy)
    dy1 = (real *) allocate(npt*sizeof(real));
  else
    dy1 = NULL;
  d1  = (real *) allocate(npt*sizeof(real));
  bpar = (real *) allocate(npar*sizeof(real));
  m  = (Moment *) allocate(npar*sizeof(Moment));

  for (i=0; i<npar; i++) {
    bpar[i] = fpar[i];
    ini_moment(&m[i],2,0);
  }
  
  for (j=0; j<nboot; j++) {
    random_permute3(npt,perm);
    for (i=0; i<npt; i++) {
      x1[i] = x[perm[i]];
      y1[i] = y[perm[i]];
      if (dy) dy1[i] = dy[perm[i]];
    }
    nrt = (*my_nllsqfit)(x1,ndim,y1,dy1,d1,npt,fpar,epar,mpar,npar,tol,itmax,lab, fitfunc,fitderv);
    dprintf(1,"%g %g %g %g\n", fpar[0],fpar[1],epar[0],epar[1]);
    for (i=0; i<npar; i++) {
      accum_moment(&m[i],fpar[i],1.0);
    }
  }
  printf("bootstrap3= ");
  for (i=0; i<npar; i++)
    printf("%g %g ",mean_moment(&m[i]),sigma_moment(&m[i]));
  printf("\n");

  free(y1);
  free(d1);
  free(bpar);
  free(perm);
  free(m);
}

/*
 * LINE:     y = a + b * x
 *    See also: poly(order=1)
 */

do_line()
{
  real *x, *y, *dy, *d, *y1, *d1;
  int i,j, nrt, npt1, iter, mpar[2], *perm;
  real fpar[2], epar[2], sigma, s, bpar[2];
  int lpar = 2;
    
  if (nxcol < 1) error("nxcol=%d",nxcol);
  if (nycol < 1) error("nycol=%d",nycol);
  if (tol < 0) tol = 0.0;
  if (lab < 0) lab = 0.0;
  sprintf(fmt,"Fitting a+bx:  \na= %s %s \nb= %s %s\nx0= %s %s\ny0= %s %s\n", 
	  format,format,format,format,format,format,format,format);

  x = xcol[0].dat;
  y = ycol[0].dat;
  dy = (dycolnr>0 ? dycol.dat : NULL);
  d = (real *) allocate(npt * sizeof(real));
  
  for (i=0; i<lpar; i++) {
    mpar[i] = mask[i];
    fpar[i] = par[i];
  }

  fitfunc = func_line;
  fitderv = derv_line;

  for (iter=0; iter<=msigma; iter++) {
    nrt = (*my_nllsqfit)(x,1,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,itmax,lab, fitfunc,fitderv);
    if (nrt==-2)
      warning("No free parameters");
    else if (nrt<0)
      error("Bad fit, nrt=%d",nrt);
    printf("nrt=%d\n",nrt);
    printf(fmt, fpar[0],epar[0],fpar[1],epar[1],
	   -fpar[0]/fpar[1], sqrt( sqr(epar[0]/fpar[1]) + sqr(epar[1]*fpar[0]/(fpar[1]*fpar[1])) ),
	   fpar[0],epar[0]);

    npt1 = remove_data(x,1,y,dy,d,npt,nsigma[iter]);
    if (npt1 == npt) break;
    npt = npt1;
  }
  bootstrap(nboot, npt,1,x,y,dy,d, lpar,fpar,epar,mpar);

  if (outstr)
    for (i=0; i<npt; i++)
      fprintf(outstr,"%g %g %g\n",x[i],y[i],d[i]);
  printf("rms/chi = %g\n",data_rms(npt,d,dy,2));
  free(d);
}


/*
 * PLANE:       y = b_0 + b_1 * x_1 + b_2 * x_2 + ... + b_n * x_n
 *
 *      used:   n = dimensionality of space in which hyper plane is fit
 */
 
/*  BUG:  plane,order=1 does not reproduce line   !!  */

do_plane()
{
  real *x1, *x2, *x, *y, *dy, *d;
  real **xp;
  int i,j,k,nrt, npt1,iter,mpar[MAXPAR];
  real fpar[MAXPAR], epar[MAXPAR];
  int lpar = order+1;

  if (nxcol != order) error("Need exactly %d value(s) for xcol=",order);
  if (nycol<1) error("Need 1 value for ycol=");
  if (tol < 0) tol = 0.0;
  if (lab < 0) lab = 0.0;
  sprintf(fmt,"p%%d= %s %s\n",format,format);

  y = ycol[0].dat;
  dy = (dycolnr>0 ? dycol.dat : NULL);
  d = (real *) allocate(npt * sizeof(real));
  x = (real *) allocate(order * npt * sizeof(real));
  for (i=0, j=0; i<npt; i++) {
    for (k=0; k<order; k++)
      x[j++] = xcol[k].dat[i];
  }

  for (i=0; i<lpar; i++) {
    mpar[i] = mask[i];
    fpar[i] = par[i];
  }

  fitfunc = func_plane;
  fitderv = derv_plane;

  for (iter=0; iter<=msigma; iter++) {
    /* should the 2 be order ?? */
    nrt = (*my_nllsqfit)(x,2,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,itmax,lab, fitfunc,fitderv);
    printf("nrt=%d\n",nrt);
    printf("Fitting p0+p1*x1+p2*x2+.....pN*xN: (N=%d)\n",order);
    for (k=0; k<lpar; k++)
      printf(fmt, k,fpar[k],epar[k]);
    if (nrt==-2)
      warning("No free parameters");
    else if (nrt<0)
      error("Bad fit, nrt=%d",nrt);
#if 0    
    npt1 = remove_data(x,2,y,dy,d,npt,nsigma[iter]);
    if (npt1 == npt) iter=msigma+1;       /* signal early bailout */
    npt = npt1;
#else
    break;
#endif
  }
  bootstrap(nboot, npt,2,x,y,dy,d, lpar,fpar,epar,mpar);

  if (outstr)
    for (i=0; i<npt; i++)
      /* bad : x1 and x2 not initialized */
      fprintf(outstr,"%g %g %g %g\n",x1[i],x2[i],y[i],d[i]);

#if 0
  outdparam("a",fpar[0]);
  outdparam("b",fpar[0]);
  outdparam("c",fpar[0]);
  outdparam("siga",epar[0]);
  outdparam("sigb",epar[0]);
  outdparam("sigc",epar[0]);
#endif
}

/*
 * GAUSS1d:       y = a + b * exp( - (x-c)^2/(2*d^2) )
 *
 */
 
do_gauss1d()
{
  real *x1, *x2, *x, *y, *dy, *d;
  int i,j, nrt, npt1, iter, mpar[4];
  real sum, dmin, dmax, xmin, xmax, fpar[4], epar[4];
  int lpar = 4;

  if (nxcol < 1) error("nxcol=%d",nxcol);
  if (nycol<1) error("Need 1 value for ycol=");
  if (tol < 0) tol = 0.0;
  if (lab < 0) lab = 0.01;
  sprintf(fmt,"Fitting a+b*exp(-(x-c)^2/(2*d^2)):  \na= %s %s \nb= %s %s \nc= %s %s\nd= %s %s\n",
	  format,format,format,format,format,format,format,format);

  x = xcol[0].dat;
  y = ycol[0].dat;
  dy = (dycolnr>0 ? dycol.dat : NULL);
  d = (real *) allocate(npt * sizeof(real));

  if (npar==0) {
    warning("No initial estimates for gauss1d, attempting to get some");
    /* tricky, this assumes X is sorted */
    par[0] = par[1] = y[0];
    par[2] = x[0];
    par[3] = x[1]-x[0];
    if (par[3] < 0) warning("Xcol not sorted, estimates may be lousy");
    sum = 0.0;
    for (i=1; i<npt; i++) {
      sum += y[i]*(x[i]-x[i-1]);          /* sum of emission */
      if (y[i] < par[0]) {                /* store min */
	par[0] = y[i];
	xmin = x[i];
      }
      if (y[i] > par[1]) {                /* store max + loc */
	par[1] = y[i];
	xmax   = x[i];
      }
    }
    dmin = 0.5*(y[0]+y[npt-1]) - par[0];
    dmax = 0.5*(y[0]+y[npt-1]) - par[1];
    dmin = ABS(dmin);
    dmax = ABS(dmax);
    printf("par01,dmin/max = %g %g %g %g\n",par[0],par[1],dmin,dmax);
    if (dmax > dmin) {                    /* positive peak */
      par[1] -= par[0];
      par[2] = xmax;
    } else {                              /* negative peak */
      dmin = par[1];
      dmax = par[0]-par[1];
      par[0] = dmin;
      par[1] = dmax;
      par[2] = xmin;
    }

    par[3] = sum / (par[1]*sqrt(TWO_PI)); /* sigma */
    printf("par=%g,%g,%g,%g\n",par[0],par[1],par[2],par[3]);
  }
  
  for (i=0; i<lpar; i++) {
    mpar[i] = mask[i];
    fpar[i] = par[i];
  }

  fitfunc = func_gauss1d;
  fitderv = derv_gauss1d;

  for (iter=0; iter<=msigma; iter++) {
    nrt = (*my_nllsqfit)(x,1,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,itmax,lab, fitfunc,fitderv);
    printf("nrt=%d\n",nrt);
    printf(fmt,fpar[0],epar[0],fpar[1],epar[1],fpar[2],epar[2],fpar[3],epar[3]);
    if (nrt==-2)
      warning("No free parameters");
    else if (nrt<0)
      error("Bad fit, nrt=%d",nrt);

    npt1 = remove_data(x,1,y,dy,d,npt,nsigma[iter]);
    if (npt1 == npt) iter=msigma+1;       /* signal early bailout */
    npt = npt1;
  }
  bootstrap(nboot, npt,1,x,y,dy,d, lpar,fpar,epar,mpar);    

  if (outstr)
    for (i=0; i<npt; i++)
      fprintf(outstr,"%g %g %g\n",x[i],y[i],d[i]);
  printf("rms/chi = %g\n",data_rms(npt,d,dy,4));
}

/*
 * GAUSS2d:       y = a + b * exp( -[(x-c)^2+(y-d)^2]/(2*e^2) )
 *
 */
 
do_gauss2d()
{
  real *x1, *x2, *x, *y, *dy, *d;
  int i,j,k, nrt, npt1, iter, mpar[5];
  real fpar[5], epar[5];
  int lpar = 5;

  if (nxcol != 2) error("bad nxcol=%d",nxcol);
  if (nycol<1) error("Need 1 value for ycol=");
  if (tol < 0) tol = 0.0;
  if (lab < 0) lab = 0.01;
  sprintf(fmt,"Fitting a+b*exp(-[(x-c)^2+(y-d)^2]/(2*e^2)):  \n"
	  "a= %s %s \nb= %s %s \nc= %s %s\nd= %s %s\ne= %s %s\n",
	  format,format,format,format,format,format,format,format,format,format);

  y = ycol[0].dat;
  dy = (dycolnr>0 ? dycol.dat : NULL);
  d = (real *) allocate(npt * sizeof(real));
  x = (real *) allocate(2 * npt * sizeof(real));
  for (i=0, j=0; i<npt; i++) {
    for (k=0; k<2; k++)
      x[j++] = xcol[k].dat[i];
  }

  for (i=0; i<lpar; i++) {
    mpar[i] = mask[i];
    fpar[i] = par[i];
  }
  
  fitfunc = func_gauss2d;
  fitderv = derv_gauss2d;

  for (iter=0; iter<=msigma; iter++) {
    nrt = (*my_nllsqfit)(x,2,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,itmax,lab, fitfunc,fitderv);
    printf("nrt=%d\n",nrt);
    printf(fmt,fpar[0],epar[0],fpar[1],epar[1],fpar[2],epar[2],fpar[3],epar[3],fpar[4],epar[4]);
    if (nrt==-2)
      warning("No free parameters");
    else if (nrt<0)
      error("Bad fit, nrt=%d",nrt);

    npt1 = remove_data(x,2,y,dy,d,npt,nsigma[iter]);
    if (npt1 == npt) iter=msigma+1;       /* signal early bailout */
    npt = npt1;
  }
  bootstrap(nboot, npt,2,x,y,dy,d, lpar,fpar,epar,mpar);    

  if (outstr)
    for (i=0; i<npt; i++)
      fprintf(outstr,"%g %g %g\n",x[i],y[i],d[i]);
  printf("rms/chi = %g\n",data_rms(npt,d,dy,5));
}


/*
 * EXP:       y = a + b * exp(-(x-c)/d)
 *
 */
 
do_exp()
{
  real *x1, *x2, *x, *y, *dy, *d;
  int i,j, nrt, npt1, iter, mpar[4];
  real fpar[4], epar[4];
  int lpar = 4;

  warning("fit=exp does not seem to work when all parameters free");

  if (nxcol < 1) error("nxcol=%d",nxcol);
  if (nycol<1) error("Need 1 value for ycol=");
  if (tol < 0) tol = 0.0;
  if (lab < 0) lab = 0.01;
  sprintf(fmt,"Fitting a+b*exp(-(x-c)/d):  \na= %s %s \nb= %s %s \nc= %s %s\nd= %s %s\n",
	  format,format,format,format,format,format,format,format);
  x = xcol[0].dat;
  y = ycol[0].dat;
  dy = (dycolnr>0 ? dycol.dat : NULL);
  d = (real *) allocate(npt * sizeof(real));

  for (i=0; i<lpar; i++) {
    mpar[i] = mask[i];
    fpar[i] = par[i];
  }
  
  fitfunc = func_exp;
  fitderv = derv_exp;

  for (iter=0; iter<=msigma; iter++) {
    nrt = (*my_nllsqfit)(x,1,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,itmax,lab, fitfunc,fitderv);
    printf("nrt=%d\n",nrt);
    printf(fmt,fpar[0],epar[0],fpar[1],epar[1],fpar[2],epar[2],fpar[3],epar[3]);
    if (nrt==-2)
      warning("No free parameters");
    else if (nrt<0)
      error("Bad fit, nrt=%d",nrt);

    npt1 = remove_data(x,1,y,dy,d,npt,nsigma[iter]);
    if (npt1 == npt) iter=msigma+1;       /* signal early bailout */
    npt = npt1;
  }
  bootstrap(nboot, npt,1,x,y,dy,d, lpar,fpar,epar,mpar);

  if (outstr)
    for (i=0; i<npt; i++)
      fprintf(outstr,"%g %g %g\n",x[i],y[i],d[i]);

}


/*
 * GROW:       y = a * (1-exp(-x/b))
 *
 */
 
do_grow()
{
  real *x1, *x2, *x, *y, *dy, *d;
  int i,j, nrt, npt1, iter, mpar[2];
  real fpar[2], epar[2];
  int lpar = 2;

  warning("fit=grow does not seem to work when all parameters free");

  if (nxcol < 1) error("nxcol=%d",nxcol);
  if (nycol<1) error("Need 1 value for ycol=");
  if (tol < 0) tol = 0.0;
  if (lab < 0) lab = 0.01;
  sprintf(fmt,"Fitting a*(1-exp(-x/b)):  \na= %s %s \nb= %s %s \n",
	  format,format,format,format);
  x = xcol[0].dat;
  y = ycol[0].dat;
  dy = (dycolnr>0 ? dycol.dat : NULL);
  d = (real *) allocate(npt * sizeof(real));

  for (i=0; i<lpar; i++) {
    mpar[i] = mask[i];
    fpar[i] = par[i];
  }
  
  fitfunc = func_grow;
  fitderv = derv_grow;

  for (iter=0; iter<=msigma; iter++) {
    nrt = (*my_nllsqfit)(x,1,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,itmax,lab, fitfunc,fitderv);
    printf("nrt=%d\n",nrt);
    printf(fmt,fpar[0],epar[0],fpar[1],epar[1],fpar[2],epar[2],fpar[3],epar[3]);
    if (nrt==-2)
      warning("No free parameters");
    else if (nrt<0)
      error("Bad fit, nrt=%d",nrt);

    npt1 = remove_data(x,1,y,dy,d,npt,nsigma[iter]);
    if (npt1 == npt) iter=msigma+1;       /* signal early bailout */
    npt = npt1;
  }
  //bootstrap(nboot, npt,1,x,y,dy,d, lpar,fpar,epar,mpar);

  if (outstr)
    for (i=0; i<npt; i++)
      fprintf(outstr,"%g %g %g\n",x[i],y[i],d[i]);

}

/*
 * POLYNOMIAL:  y = b_0 + b_1 * x^1 + b_2 * x^2 + ... + b_n * x^n
 *
 *      used:   n = order of polynomial
 */

do_poly()
{
  real *x, *y, *dy, *d;
  int i,j, nrt, npt1, iter, mpar[MAXPAR];
  real fpar[MAXPAR], epar[MAXPAR];
  int lpar = order+1;
    
  if (nxcol < 1) error("nxcol=%d",nxcol);
  if (nycol < 1) error("nycol=%d",nycol);
  if (tol < 0) tol = 0.0;
  if (lab < 0) lab = 0.0;
  sprintf(fmt,"p%%d= %s %s\n",format,format);

  x = xcol[0].dat;
  y = ycol[0].dat;
  dy = (dycolnr>0 ? dycol.dat : NULL);
  d = (real *) allocate(npt * sizeof(real));
  
  for (i=0; i<lpar; i++) {
    mpar[i] = mask[i];
    fpar[i] = par[i];
  }

  fitfunc = func_poly;
  fitderv = derv_poly;

  for (iter=0; iter<=msigma; iter++) {
    nrt = (*my_nllsqfit)(x,1,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,itmax,lab, fitfunc,fitderv);
    printf("nrt=%d\n",nrt);
    printf("Fitting p0+p1*x+p2*x^2+.....pN*x^N: (N=%d)\n",order);
    for (i=0; i<lpar; i++)
      printf(fmt,i,fpar[i],epar[i]);
    if (nrt==-2)
      warning("No free parameters");
    else if (nrt<0)
      error("Bad fit, nrt=%d",nrt);
    npt1 = remove_data(x,1,y,dy,d,npt,nsigma[iter]);
    if (npt1 == npt) iter=msigma+1;       /* signal early bailout */
    npt = npt1;
  }
  bootstrap(nboot, npt,1,x,y,dy,d, lpar,fpar,epar,mpar);

  if (outstr)
    for (i=0; i<npt; i++)
      fprintf(outstr,"%g %g %g\n",x[i],y[i],d[i]);
}


/*
 * POLYNOMIAL2:  y = a ( 1 + b*(x-x0) + c*(x-x0)^2)
 *
 *   hardcoded polynomial
 */

do_poly2()
{
  real *x, *y, *dy, *d;
  int i,j, nrt, npt1, iter, mpar[MAXPAR];
  real fpar[MAXPAR], epar[MAXPAR];
  int lpar = 4;
  
  order = 2;
  warning("testing a new poly2 mode, order=%d",order);
    
  if (nxcol < 1) error("nxcol=%d",nxcol);
  if (nycol < 1) error("nycol=%d",nycol);
  if (tol < 0) tol = 0.0;
  if (lab < 0) lab = 0.0;
  sprintf(fmt,"p%%d= %s %s\n",format,format);

  x = xcol[0].dat;
  y = ycol[0].dat;
  dy = (dycolnr>0 ? dycol.dat : NULL);
  d = (real *) allocate(npt * sizeof(real));
  
  for (i=0; i<lpar; i++) {
    mpar[i] = mask[i];
    fpar[i] = par[i];
  }

  fitfunc = func_poly2;
  fitderv = derv_poly2;

  for (iter=0; iter<=msigma; iter++) {
    nrt = (*my_nllsqfit)(x,1,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,itmax,lab, fitfunc,fitderv);
    printf("nrt=%d\n",nrt);
    printf("Fitting p1(1+p2*(x-p0)+p3*(x-p0)^2): (fixed order=%d)\n",order);
    for (i=0; i<lpar; i++)
      printf(fmt,i,fpar[i],epar[i]);
    if (nrt==-2)
      warning("No free parameters");
    else if (nrt<0)
      error("Bad fit, nrt=%d",nrt);
    npt1 = remove_data(x,1,y,dy,d,npt,nsigma[iter]);
    if (npt1 == npt) iter=msigma+1;       /* signal early bailout */
    npt = npt1;
  }
  bootstrap(nboot, npt,1,x,y,dy,d, lpar,fpar,epar,mpar);

  if (outstr)
    for (i=0; i<npt; i++)
      fprintf(outstr,"%g %g %g\n",x[i],y[i],d[i]);
}



/*
 * ARM:     y = a + b * cos [(x-c)/DPR]
 *      simple m=1 fourier analysis
 *  BUG:  for some reason the AMPHASE model does just not work. Write it out
 *        and fit it linearly as a separate cos and sin term, it works fine.
 */

do_arm()
{
  real *x, *y, *dy, *d;
  int i,j, nrt, npt1, iter, mpar[3];
  real fpar[3], epar[3],amp,pha,amperr,phaerr;
  int lpar = 3;
    
  if (nxcol < 1) error("nxcol=%d",nxcol);
  if (nycol < 1) error("nycol=%d",nycol);
  if (tol < 0) tol = 0.0;
  if (lab < 0) lab = 0.0;

  x = xcol[0].dat;
  y = ycol[0].dat;
  dy = (dycolnr>0 ? dycol.dat : NULL);
  d = (real *) allocate(npt * sizeof(real));
  
  for (i=0; i<lpar; i++) {
    mpar[i] = mask[i];
    fpar[i] = par[i];
  }

  fitfunc = func_arm;
  fitderv = derv_arm;

  
  for (iter=0; iter<=msigma; iter++) {
    nrt = (*my_nllsqfit)(x,1,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,itmax,lab, fitfunc,fitderv);
    /*if (nrt==-2)
      warning("No free parameters");
    else if (nrt<0)
      error("Bad fit, nrt=%d",nrt);
    
      printf("nrt=%d\n",nrt); */
#ifdef PHASEAMP
      printf("Fitting a+b.cos(x-c)/DPR:  \na= %g %g \nb= %g %g\nc= %g %g\n", 
      fpar[0],epar[0],fpar[1],epar[1],fpar[2],epar[2]);
#else
    printf("Fitting a0 + c1.cos(x/DPR) + s1.sin(y/DPR): ");
    printf("[x now in degrees]  \na0= %g %g \nc1= %g %g\ns1= %g %g\n", 
	 fpar[0],epar[0],fpar[1],epar[1],fpar[2],epar[2]);
    printf("Converting to amp/phase: (Fitting as a0 + a1.cos((x-p1)/DPR)\n");
    amp = sqrt(sqr(fpar[1])+sqr(fpar[2]));
    pha = atan2(fpar[2],fpar[1])*DPR;
    amperr = sqrt(sqr(fpar[1]*epar[1]) + sqr(fpar[2]*epar[2]))/amp;
    phaerr = sqrt(sqr(fpar[1]*epar[2]) + sqr(fpar[2]*epar[1]))/(amp*amp)*DPR;
    printf("a1= %g %g\n",amp,amperr); 
    printf("p1= %g %g\n",pha,phaerr);
#endif

    npt1 = remove_data(x,1,y,dy,d,npt,nsigma[iter]);
    if (npt1 == npt) iter=msigma+1;       /* signal early bailout */
    npt = npt1;
  }
  bootstrap(nboot, npt,1,x,y,dy,d, lpar,fpar,epar,mpar);


  if (outstr)
    for (i=0; i<npt; i++)
      fprintf(outstr,"%g %g %g\n",x[i],y[i],d[i]);

}

do_loren()
{
  real *x, *y, *dy, *d, *y1, *d1;
  int i,j, nrt, npt1, iter, mpar[2], *perm;
  real fpar[2], epar[2], sigma, s, bpar[2];
  int lpar = 2;
    
  if (nxcol < 1) error("nxcol=%d",nxcol);
  if (nycol < 1) error("nycol=%d",nycol);
  if (tol < 0) tol = 0.0;
  if (lab < 0) lab = 0.01;
  sprintf(fmt,"Fitting (a/PI) / ((x-b)^2 + a^2):  \na= %s %s \nb= %s %s\n", format,format,format,format);

  x = xcol[0].dat;
  y = ycol[0].dat;
  dy = (dycolnr>0 ? dycol.dat : NULL);
  d = (real *) allocate(npt * sizeof(real));
  
  for (i=0; i<lpar; i++) {
    mpar[i] = mask[i];
    fpar[i] = par[i];
  }

  fitfunc = func_loren;
  fitderv = derv_loren;

  for (iter=0; iter<=msigma; iter++) {
    nrt = (*my_nllsqfit)(x,1,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,itmax,lab, fitfunc,fitderv);
    if (nrt==-2)
      warning("No free parameters");
    else if (nrt<0)
      error("Bad fit, nrt=%d",nrt);
    printf("nrt=%d\n",nrt);
    printf(fmt, fpar[0],epar[0],fpar[1],epar[1]);

    npt1 = remove_data(x,1,y,dy,d,npt,nsigma[iter]);
    if (npt1 == npt) break;
    npt = npt1;
  }
  bootstrap(nboot, npt,1,x,y,dy,d, lpar,fpar,epar,mpar);

  if (outstr)
    for (i=0; i<npt; i++)
      fprintf(outstr,"%g %g %g\n",x[i],y[i],d[i]);
  printf("rms/chi = %g\n",data_rms(npt,d,dy,2));
  free(d);

}



do_arm3()
{
  real *x, *y, *dy, *d;
  int i,j, nrt, npt1, iter, mpar[5];
  real fpar[5], epar[5],amp,pha,amperr,phaerr,amp3,pha3,amperr3,phaerr3;
  int lpar = 5;
    
  if (nxcol < 1) error("nxcol=%d",nxcol);
  if (nycol < 1) error("nycol=%d",nycol);
  if (tol < 0) tol = 0.0;
  if (lab < 0) lab = 0.0;

  x = xcol[0].dat;
  y = ycol[0].dat;
  dy = (dycolnr>0 ? dycol.dat : NULL);
  d = (real *) allocate(npt * sizeof(real));
  
  for (i=0; i<lpar; i++) {
    mpar[i] = mask[i];
    fpar[i] = par[i];
  }

  fitfunc = func_arm3;
  fitderv = derv_arm3;

  for (iter=0; iter<=msigma; iter++) {
    nrt = (*my_nllsqfit)(x,1,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,itmax,lab, fitfunc,fitderv);
    if (nrt==-2)
      warning("No free parameters"); 
    else if (nrt<0)
      error("Bad fit, nrt=%d",nrt); 
    
    printf("nrt=%d\n",nrt);  
    printf("Fitting a0 + c1.cos(x/DPR) + s1.sin(y/DPR) + c3.sin(3x/DPR) + s3.cos(3x/DPR): ");
    printf("[x now in degrees]  \na0= %g %g \nc1= %g %g\ns1= %g %g \nc3= %g %g \ns3= %g %g \n", 
	 fpar[0],epar[0],fpar[1],epar[1],fpar[2],epar[2],fpar[3],epar[3],fpar[4],epar[4]);
    printf("Converting to amp/phase: (Fitting as a0 + a1.cos((x-p1)/DPR) + a3.cos(3(x-p3)/DPR)\n");
    amp = sqrt(sqr(fpar[1])+sqr(fpar[2]));
    pha = atan2(fpar[2],fpar[1])*DPR;
    amperr = sqrt(sqr(fpar[1]*epar[1]) + sqr(fpar[2]*epar[2]))/amp;
    phaerr = sqrt(sqr(fpar[1]*epar[2]) + sqr(fpar[2]*epar[1]))/(amp*amp)*DPR;
    printf("a1= %g %g\n",amp,amperr);
    printf("p1= %g %g\n",pha,phaerr);

    amp3 = sqrt(sqr(fpar[3])+sqr(fpar[4]));
    pha3 = atan2(fpar[4],fpar[3])*DPR/3;
    amperr3 = sqrt(sqr(fpar[3]*epar[3]) + sqr(fpar[4]*epar[4]))/amp3;
    phaerr3 = sqrt(sqr(fpar[3]*epar[4]) + sqr(fpar[4]*epar[3]))/(amp3*amp3)*DPR/3;
    printf("a3= %g %g\n",amp3,amperr3);
    printf("p3= %g %g\n",pha3,phaerr3);


    npt1 = remove_data(x,1,y,dy,d,npt,nsigma[iter]);
    if (npt1 == npt) iter=msigma+1;       /* signal early bailout */
    npt = npt1;
  }
  bootstrap(nboot, npt,1,x,y,dy,d, lpar,fpar,epar,mpar);


  if (outstr)
    for (i=0; i<npt; i++)
      fprintf(outstr,"%g %g %g\n",x[i],y[i],d[i]);

}


/*
 * PSF:       z = a * x^b * sin(y)^c + d
 *
 */
 
do_psf()
{
  real *x1, *x2, *x, *y, *dy, *d;
  int i,j,k, nrt, npt1, iter, mpar[4];
  real fpar[4], epar[4];
  int lpar = 4;

  if (nxcol != 2) error("bad nxcol=%d",nxcol);
  if (nycol<1) error("Need 1 value for ycol=");
  if (tol < 0) tol = 0.0;
  if (lab < 0) lab = 0.01;
  sprintf(fmt,"Fitting a * x^b * sin(y)^c + d::  \n"
	  "a= %s %s \nb= %s %s \nc= %s %s\nd= %s %s\n",
	  format,format,format,format,format,format,format,format);

  y = ycol[0].dat;
  dy = (dycolnr>0 ? dycol.dat : NULL);
  d = (real *) allocate(npt * sizeof(real));
  x = (real *) allocate(2 * npt * sizeof(real));
  for (i=0, j=0; i<npt; i++) {
    for (k=0; k<2; k++)
      x[j++] = xcol[k].dat[i];
  }

  for (i=0; i<lpar; i++) {
    mpar[i] = mask[i];
    fpar[i] = par[i];
    epar[i] = 0.0;
  }
  
  fitfunc = func_psf;
  fitderv = derv_psf;

  for (iter=0; iter<=msigma; iter++) {
    nrt = (*my_nllsqfit)(x,2,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,itmax,lab, fitfunc,fitderv);
    printf("nrt=%d\n",nrt);
    printf(fmt,fpar[0],epar[0],fpar[1],epar[1],fpar[2],epar[2],fpar[3],epar[3]);
    if (nrt==-2)
      warning("No free parameters");
    else if (nrt<0)
      error("Bad fit, nrt=%d",nrt);

    npt1 = remove_data(x,2,y,dy,d,npt,nsigma[iter]);
    if (npt1 == npt) iter=msigma+1;       /* signal early bailout */
    npt = npt1;
  }
  bootstrap(nboot, npt,2,x,y,dy,d, lpar,fpar,epar,mpar);    

  if (outstr)
    for (i=0; i<npt; i++)
      fprintf(outstr,"%g %g %g\n",x[i],y[i],d[i]);
  printf("rms/chi = %g\n",data_rms(npt,d,dy,4));
}
