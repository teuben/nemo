/*WIP: adding mode=plane */
/*TODO:adding mode=poly*/
/*BUG: large tables will fail in old table input method */
/*
 * TABNLLSQFIT: a general (non-linear) fitting program for tabular data 
 *
 *      12-jul-02  derived from tablsqfit, but using nllsqfit() now
 *
 */

#include <stdinc.h>  
#include <getparam.h>

string defv[] = {
    "in=???\n           input (table) file name",
    "xcol=1\n           column(s) for x, the independant variable(s)",
    "ycol=2\n           column(s) for y, the dependant variable(s)",
    "dycol=\n           column for sigma-y, if used",
    "xrange=\n          in case restricted range is used (1D only)",
    "fit=line\n         fitmode (line, ellipse, imageshift, plane, poly, peak, area)",
    "order=0\n		Order of plane/poly fit",
    "out=\n             optional output file for some fit modes",
    "nsigma=-1\n        delete points more than nsigma away?",
    "par=\n             initial estimates of parameters, if truely non-linear fit",
    "free=\n            free(1) or fixed(0) parameters?",
    "npar=\n            override number of parameters",
    "nmax=10000\n       Default max allocation",
    "tab=f\n            short one-line output?",
    "VERSION=0.9\n      12-jul-02 PJT",
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

string method;              /* fit method (line, ellipse, ....) */
stream instr, outstr;       /* input / output file */


int    nmax;                /* allocated space */
int    npt;                 /* actual number of points from table */
real   nsigma;              /* fractional sigma removal */

real  a,b;                  /* fit parameters in: y=ax+b  */
int order;

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
  return p[0] + p[1]*x[0] + p[2]*x[1];
}

static void derv_plane(real *x, real *p, real *e, int np)
{
  e[0] = 1.0;
  e[1] = x[0];
  e[2] = x[1];
}



/****************************** START OF PROGRAM **********************/

nemo_main()
{

    setparams();
    read_data();

    if (scanopt(method,"line")) {
        do_line();
    } else if (scanopt(method,"ellipse")) {
        do_ellipse();
    } else if (scanopt(method,"imageshift")) {
        do_imageshift();
    } else if (scanopt(method,"plane")) {
    	do_plane();
    } else if (scanopt(method,"poly")) {
    	do_poly();
    } else if (scanopt(method,"area")) {
        do_area();
    } else if (scanopt(method,"peak")) {
    	do_peak();
    } else
        error("fit=%s invalid; try [line,ellipse,imageshift,plane,poly,peak,area]",
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
    Qtab = getbparam("tab");

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
    if (hasvalue("npar"))
      npar = getiparam("npar");
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


write_data(outstr)              /* only for straight line fit */
stream outstr;
{
#if 0    
    int i;
    real d;

    for (i=0; i<npt; i++) {
        d = y[i] - a*x[i] - b;
        fprintf (outstr,"%f %f %f\n",x[i],y[i],d);
    }
#endif    
}

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

  nrt = nllsqfit(x,1,y,dy,d,npt,  fpar,epar,mpar,2,  tol,its,lab, fitfunc,fitderv);
  printf("a+bx:  a=%g (+/-%g) b=%g (+/-%g) \n",fpar[0],epar[0],fpar[1],epar[1]);
  for (i=0; i<npt; i++)
    dprintf(1,"%g %g %g\n",x[i],y[i],d[i]);

}

/*
 * PLANE:       y = b_0 + b_1 * x_1 + b_2 * x_2 + ... + b_n * x_n
 *
 *      used:   n = dimensionality of space in which hyper plane is fit
 */
 
do_plane()
{
  real *x1, *x2, *x, *y, *dy, *d;
  int i,j, nrt, mpar[3];
  real fpar[3], epar[3];
  int its = 50;
  real tol = 0.0, lab = 0.0;

  if (nycol<1) error("Need 1 value for ycol=");
  if (nxcol<order) error("Need %d value(s) for xcol=",order);

  x1 = xcol[0].dat;
  x2 = xcol[1].dat;
  y = ycol[0].dat;
  dy = (dycolnr>0 ? dycol.dat : NULL);
  d = (real *) allocate(npt * sizeof(real));
  x = (real *) allocate(2 * npt * sizeof(real));
  for (i=0, j=0; i<npt; i++) {
    x[j++] = x1[i];
    x[j++] = x2[i];
  }
  
  mpar[0] = mpar[1] = mpar[2] = 1;    /* 1 is which parameters to fit, 0= fixed */

  fitfunc = func_plane;
  fitderv = derv_plane;

  nrt = nllsqfit(x,2,y,dy,d,npt,  fpar,epar,mpar,3,  tol,its,lab, fitfunc,fitderv);
  printf("a+bx+cy:  a=%g (+/-%g) b=%g (+/-%g) c=%g (+/- %g)\n",
	 fpar[0],epar[0],fpar[1],epar[1],fpar[2],epar[2]);

#if 0
  outdparam("a",fpar[0]);
  outdparam("b",fpar[0]);
  outdparam("c",fpar[0]);
  outdparam("siga",epar[0]);
  outdparam("sigb",epar[0]);
  outdparam("sigc",epar[0]);
#endif
}



char name[10] = "ABCDE";

do_ellipse()
{
    real *xdat, *ydat;
    real mat[5*5], vec[5], sol[5], a[6], cnt, xmean, ymean, x, y;
    real aa,bb,cc,dd,ee,pa,pp,cospp,sinpp,cospa,sinpa,ecc,al,r,ab,
         s1,s2,s3,y1,y2,y3,x0,y0, sum0, sum1, sum2, dx, dy, rr,
	 radmean, radsig, dr;
    real aaa, bbb, xp, yp, delta, sigfac;
    real siga, sigb, sigc, sigd, sige, fac1, fac2, fac3, dr1da, dr1db, dr1dc, sigr;
    int i, j;

    if (nxcol < 1) error("nxcol=%d",nxcol);
    if (nycol < 1) error("nycol=%d",nycol);
    xdat = xcol[0].dat;
    ydat = ycol[0].dat;

    if (npt < 5) {
        warning("Got %d; need minimum 5 points for an ellipse",npt);
        return 0;
    }
    xmean = ymean = 0.0;
    for (i=0; i<npt; i++) {             /* get mean of (x,y) as first guess */
      xmean += xdat[i];                 /* of center of ellips */
      ymean += ydat[i];                 /* to prevent all too much roundoff */
      dprintf(2,"Data: %f %f\n",xdat[i], ydat[i]);
    } 
    xmean /= npt;
    ymean /= npt;
    dprintf(1,"Estimate for center of ellips using %d points: %f %f\n",
					npt,xmean,ymean);
    if (hasvalue("estimate")) {
        if (nemoinpr(getparam("estimate"),vec,2) != 2) 
            error("estimate=%s needs two values for center of ellipse",
                   getparam("estimate"));
        xmean = vec[0];
        ymean = vec[1];
        dprintf(1,"Reset center of ellips: %f %f\n", xmean,ymean);
    }

    lsq_zero(5,mat,vec);

    for (i=0; i<npt; i++) {       /* gather all the stuff in matrix */
        x = xdat[i]-xmean;          /* treat (x,y) w.r.t. the mean center */
        y = ydat[i]-ymean;          /* of all points to prevent rounding err */
        a[0] = sqr(x);
        a[1] = 2*x*y;
        a[2] = sqr(y);
        a[3] = x;
        a[4] = y;
        a[5] = 1.0;
        lsq_accum(5,mat,vec,a,1.0); /* spread a[] into mat[] and vec[] */
    }

    for (i=0; i<5; i++)		/* print input matrix */
      dprintf(1,"( %9.3e %9.3e %9.3e %9.3e %9.3e  ) * ( %c ) = ( %9.3e )\n",
        mat[i*5+0], mat[i*5+1], mat[i*5+2], mat[i*5+3], mat[i*5+4],
        name[i], vec[i]);

    lsq_solve(5,mat,vec,sol);

    for (i=0; i<5; i++)		/* print inverse matrix & solution */
      dprintf(1,"( %9.3e %9.3e %9.3e %9.3e %9.3e  ) ; %c = %9.3e\n",
        mat[i*5+0], mat[i*5+1], mat[i*5+2], mat[i*5+3], mat[i*5+4],
        name[i], sol[i]);

    sigfac = 0.0;
    for (i=0; i<npt; i++) {
        x = xdat[i]-xmean;          /* treat (x,y) w.r.t. the mean center */
        y = ydat[i]-ymean;          /* of all points to prevent rounding err */
	sigfac += sqr(sol[0]*x*x+2*sol[1]*x*y+sol[2]*y*y+sol[3]*x+sol[4]*y-1);
    }
    sigfac /= (npt-5);
    sigfac = sqrt(sigfac);
    dprintf(1,"Sigma factor=%g\n",sigfac);
    siga = sigfac * sqrt(mat[0]);     /* trace elements are errors */
    sigb = sigfac * sqrt(mat[6]);
    sigc = sigfac * sqrt(mat[12]);
    sigd = sigfac * sqrt(mat[18]);
    sigd = sigfac * sqrt(mat[24]);

    /* Now that we have the coefficient a,b,c,d,e in:
     *      a.x^2 + 2bxy + cy^2 + dx + ey = 1
     * They have to be converted to human readable ones
     */

    aa = sol[0]; bb = sol[1]; cc = sol[2]; dd = sol[3]; ee = sol[4];
    dprintf(1,"Solutions  a..e: %g %g %g %g %g\n",aa,bb,cc,dd,ee);
    dprintf(1,"Errors  in a..e: %g %g %g %g %g\n",siga,sigb,sigc,sigd,sige);
    dprintf(1,"FracErr in a..e: %g %g %g %g %g\n",siga/aa,sigb/bb,sigc/cc,sigd/dd,sige/ee);

    delta = aa*cc-bb*bb;
    if (delta < 0)
        warning("You seem to have an hyperbola, not an ellipse");
    pp = atan(2.0*bb/(aa-cc));
    pa = 0.5*pp;        /* P.A. of undetermined axis */
    cospp = cos(pp);
    cospa = cos(pa);
    sinpp = sin(pp);
    sinpa = sin(pa);
    fac1 = sinpp*(aa+cc) + 2*bb;
    fac2 = sinpp*(aa+cc) - 2*bb;
    fac3 = sqr(aa-cc)+4*sqr(bb);
    r = fac1/fac2;
    r = sqrt(ABS(r));           /* r is now axis ratio b/a or a/b for ell/hyp */
    s1 = bb*ee-cc*dd;
    s2 = bb*dd-aa*ee;
    s3 = aa*cc-bb*bb;
    x0 = 0.5*s1/s3;
    y0 = 0.5*s2/s3;

    lsq_zero(3,mat,vec);
    vec[0] = aa*sqr(x0)+1; 
    vec[1] = bb*sqr(x0);
    vec[2] = cc*sqr(x0);
    lsq_cfill(3,mat,0,vec);
    vec[0] = 2*aa*x0*y0; 
    vec[1] = 2*bb*x0*y0+1;
    vec[2] = 2*cc*x0*y0; 
    lsq_cfill(3,mat,1,vec);
    vec[0] = aa*sqr(y0);
    vec[1] = bb*sqr(y0);
    vec[2] = cc*sqr(y0)+1;
    lsq_cfill(3,mat,2,vec);
    vec[0] = aa;
    vec[1] = bb;
    vec[2] = cc;
    lsq_solve(3,mat,vec,sol);
    dprintf(1,"Real A,B,C = %f %f %f \n",sol[0],sol[1],sol[2]);
    if (ABS(cospp) > ABS(sinpp))    /* two ways to find ab: which is a^2 now */
        ab = (2/(sol[0]+sol[2] + (sol[0]-sol[2])/cospp));
    else                            /* use the biggest divisor */
        ab = (2/(sol[0]+sol[2] + 2*sol[1]/sinpp));

    /* compute error in axis ratio : for this one we only need 3 partial derivitives */
#if 1
    /* first Mousumi version */
    warning("old math");
    dr1da = 4*bb*(sinpp-(2*bb*(aa+cc)*cospp)/fac3)/sqr(fac2);
    dr1dc = 4*bb*(sinpp+(2*bb*(aa+cc)*cospp)/fac3)/sqr(fac2);
#else
    /* new version, but it's really the same math, just refactorized:-) */
    warning("new math");
    dr1da = 4*bb*(sinpp-((aa+cc)*(aa-cc)*sinpp)/fac3)/sqr(fac2);
    dr1dc = 4*bb*(sinpp+((aa+cc)*(aa-cc)*sinpp)/fac3)/sqr(fac2);
#endif
    dr1db = (8*bb*cospp*(sqr(aa)-sqr(cc))+4*sinpp*(aa+cc)*fac3)/(fac3*sqr(fac2)); 
    sigr  = sqrt((sqr(siga*dr1da)+sqr(sigb*dr1db)+sqr(sigc*dr1dc))/(2*r));
    dprintf(1,"sigr=%g\n",sigr);
    dprintf(2,"aa,bb,cc,sinpp,cospp=%g %g %g %g %g\n",
  	       aa,bb,cc,sinpp,cospp);

    pa *= 180/PI;           /* PA is now in degrees */
    dprintf(1,"Before re-arranging: ab, r, pa = %f %f %f\n",ab,r,pa);
    if (delta > 0) {            /* ELLIPSE */
        if (r>1.0) { /* ellipse with r>1 means we've got a & b mixed up */
            r = 1/r;            /* so make r < 1 */
            ab = sqrt(ab)/r;    /* and ab is now the true semi major axis */
            pa -= 90.0;         /* and change PA by 90 degrees */
	    sigr *= sqr(r);     /* correct errors */
       }
       ecc = sqrt(1-r*r);       /* eccentricity of the ellipse ecc=sqrt(a^2-b^2)/a */
    } else {                    /* HYPERBOLA */
        if (ab < 0) {             /* if ab<0 means we've got a and b mixed up */
            r = 1/r;              /* so make r < 1 */
            ab = sqrt(ABS(ab))/r; /* so take abs value as semi major real axis */
            pa -= 90.0;           /* and change PA by 90 degrees */
        } else {
            ab = sqrt(ab); 
        }
        ecc = sqrt(1+r*r);      /* eccentricity of the hyperbola ecc=sqrt(a^2+b^2)/a */
    }

    x0 += xmean;        /* and finally correct to the real center */
    y0 += ymean;    

    if (nsigma<0) {         /* if no rejection of 'bad' points */
        if (delta>0) {
	  if(Qtab) {
	    printf("%f %f %f %f %f %f\n",ab,ecc,r,x0,y0,pa+90);
	  } else {
            printf("Ellips fit: \n");
            printf(" semi major axis: %g\n",ab);
            printf(" eccentricity:    %g axis ratio: %g  error: %g\n",ecc,r,sigr);
            printf(" x-center:        %g\n",x0);
            printf(" y-center:        %g\n",y0);
            printf(" P.A.:            %g\n",pa+90);
	  }
        } else {
            printf("Hyperbole fit: \n");
            printf(" semi real axis: %g\n",ab);
            printf(" eccentricity:   %g\n",ecc);
            printf(" x-center:       %g\n",x0);
            printf(" y-center:       %g\n",y0);
            printf(" P.A.:           %g\n",pa+90);
        }
    } else {
        dprintf(1,"Ellfit[%d]: a,b/a,x,y,pa=%g %g %g %g %g\n",npt,ab,r,x0,y0,pa+90);
			
        dprintf(0,"Now looking to delete points > %f sigma\n",nsigma);
        sum0 = sum1 = sum2 = 0.0;
        pa *= PI/180;           /* pa now in radians */
        for (i=0; i<npt; i++) {
            x = xdat[i] - x0;
            y = ydat[i] - y0;
            pp = atan2(y,x)-pa;
            rr = sqrt((sqr(x)+sqr(y))*(sqr(cos(pp))+sqr(sin(pp)/r)));
            sum0 += 1.0;
            sum1 += rr;
            sum2 += sqr(rr);
        }
        radmean = sum1/sum0;
        radsig = sqrt(sum2/sum0 - sqr(radmean));
	dprintf(0,"Deprojected radii:  mean: %f  sigma: %f\n",
			radmean, radsig);
	
        for (i=0; i<npt; i++) {
            x = xdat[i] - x0;
            y = ydat[i] - y0;
            pp = atan2(y,x)-pa;
            rr = sqrt((sqr(x)+sqr(y))*(sqr(cos(pp))+sqr(sin(pp)/r)));
            dr = (rr - radmean)/radsig;
            if (ABS(dr) < nsigma)
                printf("%g %g %g\n",xdat[i],ydat[i],dr);
            else
                dprintf(0,"Data Point %d deviates %f sigma from mean\n",
			i,dr);
        }
#if 0
        for (i=0; i<npt; i++) {
            pp = atan2(ydat[i]-y0,xdat[i]-x0)*180.0/PI;
            x = xdat[i] - xmean;
            y = ydat[i] - ymean;
            rr = aa*x*x + 2*bb*x*y + cc*y*y + dd*x + ee*y;
            printf("%d %g %g\n",i+1,(pp<0?pp+360:pp),rr);
        }
#endif    
    }

    if (outstr) {
       if (delta<0) error("Can't compute output hyperbola table yet");
       bbb = r*ab;        /* b = semi minor axis  */
       aaa = 1-r*r;            /* eps^2 */
       for (i=0; i<=360; i++) {
           cospp = cos(PI*i/180.0);
           rr = bbb/sqrt(1-aaa*cospp*cospp);
           xp = x0 + rr*cos(PI*(i+pa)/180.0);
           yp = y0 + rr*sin(PI*(i+pa)/180.0);
           fprintf(outstr,"%d %g %g\n", i+1, xp, yp);
       }
    }
}

do_imageshift()
{
    real *x, *y, *u, *v;
    /* this code is on 3b1 ... which we've lot by now */
}


my_poly(bool Qpoly)
{
  real mat[(MAXCOL+1)*(MAXCOL+1)], vec[MAXCOL+1], sol[MAXCOL+1], a[MAXCOL+2], sum;
  int i, j;

  if (nycol<1) error("Need 1 value for ycol=");
  if (nxcol<order && !Qpoly) error("Need %d value(s) for xcol=",order);

  lsq_zero(order+1, mat, vec);
  for (i=0; i<npt; i++) {
    a[0] = 1.0;
    for (j=0; j<order; j++) {
      if (Qpoly)
	a[j+1] = a[j] * xcol[0].dat[i];     /* polynomial */
      else
	a[j+1] = xcol[j].dat[i];            /* plane */
    }
    a[order+1] = ycol[0].dat[i];
    lsq_accum(order+1,mat,vec,a,1.0);
  }
  if (order==0) printf("TEST = %g %g\n",mat[0], vec[0]);
  lsq_solve(order+1,mat,vec,sol);
  printf("%s fit of order %d:\n", Qpoly ? "Polynomial" : "Planar" , order);
  for (j=0; j<=order; j++) printf("%g ",sol[j]);
  printf("\n");

  if (outstr && Qpoly) {           /* output fitted values, if need be */
    for(i=0; i<npt; i++) {
      sum=sol[order];
      for (j=order-1; j>=0; j--)
	sum = (xcol[0].dat[i] * sum + sol[j]);
      fprintf(outstr,"%g %g %g\n",xcol[0].dat[i], ycol[0].dat[i], sum);
    }
  }
}

/*
 * POLYNOMIAL:  y = b_0 + b_1 * x^1 + b_2 * x^2 + ... + b_n * x^n
 *
 *      used:   n = order of polynomial
 */
 
do_poly()
{
    my_poly(TRUE);
}



/* fit a peak, 
 * for now this is the my_poly() code, though forced with order=2 
 * first find the peak, then take the two points on either side
 * to find an exact solution
 */

do_peak()
{
    real mat[(MAXCOL+1)*(MAXCOL+1)], vec[MAXCOL+1], sol[MAXCOL+1], a[MAXCOL+2];
    real *x, *y, ymax;
    int i, j, k, range=1;

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

do_area()
{
    real *x, *y, ymax, xmean, ymean;
    int i, j, k;

    order = 2;

    if (nycol<1) error("Need 1 value for ycol=");

    x = xcol[0].dat;
    y = ycol[0].dat;
    for (i=0; i<npt; i++) {
      xmean += x[i];
      ymean += y[i];
    }
    xmean /= npt;
    ymean /= npt;

    error("ah, not done coding here yet");

    /* allocate temp array, 
       compute angles 
       sort an index array by angles
       compute half the sum of x[i+1]*y[i]-x[i]*y[i+1]
       that's the area
    */
       

}


myfit(x,y,npt,dy,mwt,b,a,sigb,siga,chi2,q)    /* NumRec emulator - no chi^2 */
real *x, *y, *dy;
int npt, mwt;
real *b, *a, *siga, *sigb, *chi2, *q;
{
    real mat[4], vec[2], sol[2], aa[3];
    int i;

    dprintf(0,"local fit\n");
    lsq_zero(2,mat,vec);
    for (i=0; i<npt; i++) {
        aa[0] = 1.0;       /* mat */
        aa[1] = x[i];
        aa[2] = y[i];       /* rsh */
        lsq_accum(2, mat,  vec, aa, 1.0);    /* 1.0 -> dy */
    }
    dprintf(0,"fit: mat = %f %f %f %f\n", mat[0], mat[1], mat[2], mat[3]);
    dprintf(0,"fit: vec = %f %f\n", vec[0], vec[1]);
    lsq_solve(2, mat, vec, sol);
    *a = sol[1];
    *b = sol[0];
}

/* NR version of fit(), using gamma functions */

extern real gammq(real a, real x);

void fit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
real x[],y[],sig[],*a,*b,*siga,*sigb,*chi2,*q;
int ndata,mwt;
{
        int i;
        real wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;



        *b=0.0;
        if (mwt) {
                ss=0.0;
                for (i=0;i<ndata;i++) {
                        wt=1.0/sqr(sig[i]);
                        ss += wt;
                        sx += x[i]*wt;
                        sy += y[i]*wt;
                }
        } else {
                for (i=0;i<ndata;i++) {
                        sx += x[i];
                        sy += y[i];
                }
                ss=ndata;
        }
        sxoss=sx/ss;
        if (mwt) {
                for (i=0;i<ndata;i++) {
                        t=(x[i]-sxoss)/sig[i];
                        st2 += t*t;
                        *b += t*y[i]/sig[i];
                }
        } else {
                for (i=0;i<ndata;i++) {
                        t=x[i]-sxoss;
                        st2 += t*t;
                        *b += t*y[i];
                }
        }
        *b /= st2;
        *a=(sy-sx*(*b))/ss;
        *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
        *sigb=sqrt(1.0/st2);
        *chi2=0.0;
        if (mwt == 0) {
                for (i=0;i<ndata;i++)
                        *chi2 += sqr(y[i]-(*a)-(*b)*x[i]);
                *q=1.0;
                sigdat=sqrt((*chi2)/(ndata-2));
                *siga *= sigdat;
                *sigb *= sigdat;
        } else {
                for (i=0;i<ndata;i++)
                        *chi2 += sqr((y[i]-(*a)-(*b)*x[i])/sig[i]);
                *q=gammq(0.5*(ndata-2),0.5*(*chi2));	
        }
}


void fit1(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
real x[],y[],sig[],*a,*b,*siga,*sigb,*chi2,*q;
int ndata,mwt;
{
  real mat[4], vec[2], sol[2], aa[3];
  int i;

  warning("testing fit using lsq_ routines");

  lsq_zero(2,mat,vec);
  for (i=0; i<ndata; i++) {
    aa[0] = 1.0;
    aa[1] = x[i];
    aa[2] = y[i];
    lsq_accum(2,mat,vec,aa,1.0);
  }
  printf("mat   : %g %g => %g \n",mat[0],mat[1],vec[0]);
  printf("      : %g %g => %g \n",mat[2],mat[3],vec[1]);
  lsq_solve(2,mat,vec,a);
  printf("a+bx  :  a=%g  b=%g\n", aa[0],aa[1]);
  printf("mat^-1: %g %g \n       %g %g\n", mat[0],mat[1],mat[2],mat[3]);


  stop(0);       
}

