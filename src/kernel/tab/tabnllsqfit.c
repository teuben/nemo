/*
 * TABNLLSQFIT: a general (non-linear) fitting program for tabular data 
 *
 *      12-jul-02  derived from tablsqfit, but using nllsqfit() now
 *
 *
 *  line       a+bx
 *  plane      p0+p1*x1+p2*x2+p3*x3+.....     up to 'order'   (2d plane has order=2)
 *  poly       p0+p1*x+p2*x^2+p3*x^3+.....    up to 'order'   (paraboloid has order=2)
 *  gauss      p0+p1*exp(-(x-p2)^2/(2*p3^2))
 *  exp        p0+p1*exp(-(x-p2)/p3)   
 */

#include <stdinc.h>  
#include <getparam.h>

string defv[] = {
    "in=???\n           input (table) file name",
    "xcol=1\n           column(s) for x, the independant variable(s)",
    "ycol=2\n           column(s) for y, the dependant variable(s)",
    "dycol=\n           optional column for sigma-y (weight = 1/dy**2)",
    "xrange=\n          in case restricted range is used (1D only)",
    "fit=line\n         fitmode (line, plane, poly, gauss, exp)",
    "order=2\n		Order of plane/poly fit",
    "out=\n             output file for some fit modes",
    "nsigma=-1\n        ** delete points more than nsigma away?",
    "par=\n             initial estimates of parameters (p0,p1,p2,...)",
    "free=\n            free(1) or fixed(0) parameters? [1,1,1,1,....]",
    "nmax=10000\n       Default max allocation",
    "VERSION=1.0b\n     13-jul-02 PJT",
    NULL
};

string usage="a non-linear least square fitting program";

/**************** SOME GLOBAL VARIABLES ************************/

#if !defined(HUGE)
#define HUGE 1e20
#endif

#define MAXCOL 10
#define MAXPAR 10

typedef struct column {
    int maxdat;     /* allocated length of data */          /* not used */
    int ndat;       /* actual length of data */             /* not used */
    real *dat;      /* pointer to data */
    int colnr;      /* column number this data came from */ /* not used */
} a_column;

int nxcol, nycol, xcolnr[MAXCOL], ycolnr[MAXCOL], dycolnr; 
a_column            xcol[MAXCOL],   ycol[MAXCOL],   dycol;

real xrange[MAXCOL*2];      /* ??? */

string method;              /* fit method (line, poly, ....) */
stream instr, outstr;       /* input / output file */


int    nmax;                /* allocated space */
int    npt;                 /* actual number of points from table */
real   nsigma;              /* fractional sigma removal */

int order;                  /* order of poly's/hyperplane's */

int mask[MAXPAR];           /* 1=free  0=fixed parameter */
real par[MAXPAR];           /* initial parameters */
int npar; 

bool Qtab;                  /* do table output ? */

rproc fitfunc;
proc  fitderv;



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
  e[2] = -p[1]*e[1] * a / (b*b);
  e[3] = p[1] * e[1] * a * a / (b*b*b);
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
  e[3] = p[1] * e[1] * arg / b;
}

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



/****************************** START OF PROGRAM **********************/

nemo_main()
{

    setparams();
    read_data();

    if (scanopt(method,"line")) {
        do_line();
    } else if (scanopt(method,"plane")) {
    	do_plane();
    } else if (scanopt(method,"poly")) {
    	do_poly();
    } else if (scanopt(method,"gauss")) {
    	do_gauss();
    } else if (scanopt(method,"exp")) {
    	do_exp();
    } else
        error("fit=%s invalid; try [line,plane,poly,gauss]",
	      getparam("fit"));

    if (outstr) strclose(outstr);
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

    if (hasvalue("xrange"))
        setrange(xrange,getparam("xrange"));
    else {
        xrange[0] = -HUGE;
        xrange[1] = HUGE;
    } 
    
    method = getparam("fit");
    nsigma = getdparam("nsigma");
    order = getiparam("order");
    if (order<0) error("order=%d of %s cannot be negative",order,method);

    if (hasvalue("free")) {
      int nfree;
      nfree = nemoinpi(getparam("free"),mask,MAXPAR);
      if (nfree < 0) error("bad free=");
      for (i=nfree; i<MAXPAR; i++)
	mask[i] = 1;
    } else {
      for (i=0; i<MAXPAR; i++)
	mask[i] = 1;
    }
    if (hasvalue("par")) {
      npar = nemoinpr(getparam("par"),par,MAXPAR);
      if (npar < 0) error("bad par=");
      for (i=npar; i<MAXPAR; i++)
	par[i] = 0.0;
    } else
      npar = 0;
}

setrange(real *rval, string rexp)
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

read_data()
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



/*
 * LINE:     y = a + b * x
 *    See also: poly(order=1)
 */

do_line()
{
  real *x, *y, *dy, *d;
  int i,j, nrt, mpar[2];
  real fpar[2], epar[2];
  int its = 50;
  real tol = 0.0, lab = 0.0;
  int lpar = 2;
    
  if (nxcol < 1) error("nxcol=%d",nxcol);
  if (nycol < 1) error("nycol=%d",nycol);
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

  nrt = nllsqfit(x,1,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,its,lab, fitfunc,fitderv);
  printf("nrt=%d\n",nrt);
  printf("a+bx:  \na= %g %g \nb= %g %g\n", fpar[0],epar[0],fpar[1],epar[1]);
  if (outstr)
    for (i=0; i<npt; i++)
      fprintf(outstr,"%g %g %g\n",x[i],y[i],d[i]);

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
  int i,j,k,nrt, mpar[MAXPAR];
  real fpar[MAXPAR], epar[MAXPAR];
  int its = 50;
  int lpar = order+1;
  real tol = 0.0, lab = 0.0;

  if (nxcol<order) error("Need %d value(s) for xcol=",order);
  if (nycol<1) error("Need 1 value for ycol=");

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

  nrt = nllsqfit(x,2,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,its,lab, fitfunc,fitderv);
  printf("nrt=%d\n",nrt);
  printf("Fitting p0+p1*x1+p2*x2+.....pN*xN: (N=%d)\n",order);
  for (k=0; k<lpar; k++)
    printf("p%d= %g %g\n",k,fpar[k],epar[k]);
  if (outstr)
    for (i=0; i<npt; i++)
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
 * GAUSS:       y = a + b * exp( - (x-c)^2/(2*d^2) )
 *
 */
 
do_gauss()
{
  real *x1, *x2, *x, *y, *dy, *d;
  int i,j, nrt, mpar[4];
  real fpar[3], epar[4];
  int its = 50;
  real tol = 0.0, lab = 0.01;
  int lpar = 4;

  if (nxcol < 1) error("nxcol=%d",nxcol);
  if (nycol<1) error("Need 1 value for ycol=");

  x = xcol[0].dat;
  y = ycol[0].dat;
  dy = (dycolnr>0 ? dycol.dat : NULL);
  d = (real *) allocate(npt * sizeof(real));

  for (i=0; i<lpar; i++) {
    mpar[i] = mask[i];
    fpar[i] = par[i];
  }
  
  fitfunc = func_gauss1d;
  fitderv = derv_gauss1d;

  nrt = nllsqfit(x,1,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,its,lab, fitfunc,fitderv);
  printf("nrt=%d\n",nrt);
  printf("Fitting a+b*exp(-(x-c)^2/(2*d^2)):  \na= %g %g \nb= %g %g \nc= %g %g\nd= %g  %g\n",
	 fpar[0],epar[0],fpar[1],epar[1],fpar[2],epar[2],fpar[3],epar[3]);
  if (nrt < 0) error("Bad gauss fit nrt=%d",nrt);
  if (outstr)
    for (i=0; i<npt; i++)
      fprintf(outstr,"%g %g %g\n",x[i],y[i],d[i]);

}


/*
 * EXP:       y = a + b * exp(-(x-c)/d)
 *
 */
 
do_exp()
{
  real *x1, *x2, *x, *y, *dy, *d;
  int i,j, nrt, mpar[4];
  real fpar[3], epar[4];
  int its = 50;
  real tol = 0.0, lab = 0.01;
  int lpar = 4;

  warning("fit=exp does not work yet");

  if (nxcol < 1) error("nxcol=%d",nxcol);
  if (nycol<1) error("Need 1 value for ycol=");

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

  nrt = nllsqfit(x,1,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,its,lab, fitfunc,fitderv);
  printf("nrt=%d\n",nrt);
  printf("Fitting a+b*exp(-(x-c)/d):  \na= %g %g \nb= %g %g \nc= %g %g\nd= %g  %g\n",
	 fpar[0],epar[0],fpar[1],epar[1],fpar[2],epar[2],fpar[3],epar[3]);
  if (nrt < 0) error("Bad exp fit nrt=%d",nrt);
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
  int i,j, nrt, mpar[MAXPAR];
  real fpar[MAXPAR], epar[MAXPAR];
  int its = 50;
  real tol = 0.0, lab = 0.0;
  int lpar = order+1;
    
  if (nxcol<order) error("Need %d value(s) for xcol=",order);
  if (nycol < 1) error("nycol=%d",nycol);

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

  nrt = nllsqfit(x,1,y,dy,d,npt,  fpar,epar,mpar,lpar,  tol,its,lab, fitfunc,fitderv);
  printf("nrt=%d\n",nrt);
  printf("Fitting p0+p1*x+p2*x^2+.....pN*x^N: (N=%d)\n",order);
  for (i=0; i<lpar; i++)
    printf("p%d= %g %g\n",i,fpar[i],epar[i]);
  if (outstr)
    for (i=0; i<npt; i++)
      fprintf(outstr,"%g %g %g\n",x[i],y[i],d[i]);
}

