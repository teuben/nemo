/*******************************************************************************
 *    ROTCURTAB prints a table of a rotation curve - big hack
 *
 *     History: 29-jan-08 : quickly hacked off rotcurshape                   pjt
 *
 ******************************************************************************/

#include <stdinc.h>
#include <getparam.h>
#include <image.h>
#include <filefn.h>
#include <extstring.h>
#include <loadobj.h>

#define GPARAM  5           /* number of geometric parameters for a disk */
#define PARAMS  50          /* (max) number of parameters */
#define MAXMOD  5           /* number of rotation curve models (should match rotcur#) */
#define MAXPAR  5           /* number of parameters per model (can be made larger) */

#define MAXRAD  1000        /* maximum number of radii */


string defv[] = {
    "radii=\n        Radii where to sample the rotation curve",
    "rotcur1=\n      Rotation curve <NAME>, and parameters",
    "rotcur2=\n      Rotation curve <NAME>, and parameters",
    "rotcur3=\n      Rotation curve <NAME>, and parameters",
    "rotcur4=\n      Rotation curve <NAME>, and parameters",
    "rotcur5=\n      Rotation curve <NAME>, and parameters",
    "load=\n         dynamically loadobject file with rotcur_<NAME>",
    "VERSION=1.0\n   29-jan-08 PJT",
    NULL,
};

string usage="prints a table of a rotation curve";

string cvsid="$Id$";


typedef real (*rcproc)(real, int, real *, real *);


int    npar[MAXMOD];            /* active number of parameters per model */
int    ipar[MAXMOD];
real   mpar[MAXMOD][MAXPAR];    /* parameters per model */
rcproc rcfn[MAXMOD];

int    nmod = 0;                /* number of models actually used */
int    nparams = 0;             /* total number of parameters (5 + # for models) */



int nrad = 0;
real rad[MAXRAD];


extern string *burststring(string,string);
extern void freestrings(string *);


real rotcur(real r);




/* a bunch of rotation curves and parameter derivatives */


real rotcur_flat(real r, int n, real *p, real *d)
{
  d[0] = 1.0;
  return p[0];
}

real rotcur_linear(real r, int n, real *p, real *d)
{
  d[0] = r;
  return p[0] * r;
}

real rotcur_poly(real r, int np, real *p, real *d)
{
  real v, dp, x = r/p[1];
  int i;

  i = np-1;
  v = 0;
  dp = 0;
  d[1] = p[0] * x;   /* fake placeholder for recursion coming up */
  while (i > 1) {    /* p[0] and p[1] are special, p[2] last one in loop */
    v = v*x + p[i];
    dp = dp*x + i*p[i];
    d[np+1-i] = d[np-i] * x;
    i--;
  }
  v  = x*(1+x*v);
  dp = x*(1+x*dp);

  d[0] = v;
  d[1] = -p[0]*dp/p[1];

  return p[0] * v;
}

/*
 *  v = p0 * x / (1+x)
 */

real rotcur_core1(real r, int n, real *p, real *d)     /* power, with c=1 fixed */
{
  real x = r / p[1];
  d[0] = x/(1+x);
  d[1] = -p[0]*d[0]/(p[1]*(1+x));
  return p[0] * d[0];
}

/*
 *  v = p0 * x / sqrt(1+x^2)
 */

real rotcur_core2(real r, int n, real *p, real *d)     /* power, with c=2 fixed */
{
  real x = r/p[1];
  real y1 = 1+x*x;
  real y2 = sqrt(y1);
  real v = x / y2;

  d[0] = v;
  d[1] = -p[0]*v/(y1*p[1]);
  return p[0] * v;
}

real rotcur_core(real r, int np, real *p, real *d)
{
  real x = r/p[1];
  real c = p[2];
  real q1 = pow(x,c);
  real q = 1+q1;
  real lnx = log(x);
  real lnq = log(q);
  real y = pow(q,1/c);

  d[0] = x / y;
  d[1] = -p[0]*d[0]/(p[1]*q);
  d[2] = (-((q1*lnx)/(c*q)) + lnq/(c*c))/y;     /* CForm[D[(1+x^c)^(-1/c),c]]  */
  d[2] *= p[0] * x;
  return p[0] * d[0];
}

real rotcur_plummer(real r, int np, real *p, real *d)
{
  real x = r/p[1];
  real y = pow(1+x*x,-0.75);
  d[0] = y;
  d[1] = -x*p[0]/p[1]*(1-x*x/2)/(1+x*x)/y;
  return p[0] * x * y;
}

#if 0
real rotcur_tanh(real r, int np, real *p, real *d)
{
  v = tanh(x);
  dvdx = sqr(sech(x));
}
#endif

/*
 * softened iso-thermal sphere: (a.k.a. pseudo-isothermal)
 *    rho  = rho0/(1+x^2)                    x = r/r0
 *    vrot = vrot0*(1-atan(x)/x)^(1/2)   ,   vrot0 = sqrt(4.pi.G.rho0*r0^2)
 *    

In[1]:= Integrate[x^2/(1+x^2),x]

Out[1]= x - ArcTan[x]


In[7]:=D[Sqrt[1-ArcTan[x]/x],x]

              1         ArcTan[x]
        -(----------) + ---------
                  2         2
          x (1 + x )       x
Out[7]= -------------------------
                     ArcTan[x]
          2 Sqrt[1 - ---------]
                         x
*/


real rotcur_iso(real r, int np, real *p, real *d)
{
  real x = r/p[1];
  real v = sqrt(1-atan(x)/x);
  d[0] = v;
  d[1] = -p[0]/p[1]*(1/(1+x*x) - v);
  return p[0] * v;
}

/*
 * used in van Moorsel & Wells, AJ 90, 1038 (1985)
 *
 *    V/Vmax = 1 - e^{-ln{100) R/Rmax}
 *  or as we write:
 *    V = Vmax ( 1 - e^{-R/Rmax} )
 *
 */

real rotcur_exp(real r, int np, real *p, real *d)
{
  real x = r/p[1];
  real y = exp(-x);
  d[0] = 1-y;
  d[1] = -p[0]*y/p[1]*x;
  return p[0] * d[0];
}

/*
 *  NFW profile:  pars = V_200,R_200,c
 *  V_c^2(x)=V_{200}^2 \frac{\ln(1+cx)-cx(1+cx)^{-1}}
 *           {x[\ln(1+c)-c(1+c)^{-1}]}
 *
 *  In[7]:=D[ Sqrt[(Log[1+c*x]-c*x/(1+c*x))/x] ,x]
 */

real rotcur_nfw(real r, int np, real *p, real *d)
{
  real x = r/p[1];
  real c = p[2];
  real cx = x*c;
  real lncx = log(1+cx);
  real a = -cx/(1+cx) + lncx;
  real v2 = a/(log(1+c)-c/(1+c))/x;
  real v=sqrt(v2);
  d[0] = v;
  d[1] = -p[0]/p[1]*x  * (sqr(c/(1+cx)) - a/x)/2/sqrt(a);
  d[2] = 0.0;  /* Note, don't allow derivatives w.r.t. c -- otherwise degenerate */
  return p[0] * d[0];
}

/*
 * Moore et al (1999, MNRAS 310, 1147:
 *
 * rho \propto  r^{-3/2}
 * i.e.
 * v   \propto  r^{1/4}      ->   see 'power' with P3=0.25
 */

real rotcur_moore(real r, int np, real *p, real *d)
{
  return 0.0;
}


/*
 * Brandt, J.C. (1960, ApJ 131, 293)
 * Brandt, J.C. & Scheer, L.S.  1965 AJ 70, 471
 *
 * v = v_0   x /   (1/3 + 3/2*x^n)^(3/2n)
 * 
 * 
 */

real rotcur_brandt(real r, int np, real *p, real *d)
{
  real x = r/p[1];
  real n = p[2];
  real dn = pow(1.0/3.0+1.5*pow(x,n), 1.5/n);
  real v = x/dn;
  d[0] = v;
  d[1] = 0.0;   /* to do !! */
  d[2] = 0.0;   /* to do !! */
  return p[0] * v;
}

/* 
 * simple Power Law rotation curve
 *
 *  v = v_0 x^a
 */

real rotcur_power(real r, int np, real *p, real *d)
{
  real x = r/p[1];
  real a = p[2];
  real v = pow(x,a);
  d[0] = v;
  d[1] = -a*p[0]*v/p[1];
  d[2] = p[0]*v*log(x);
  return p[0] * d[0];
}

/*
 * some kind of toy disk for max disk degeneracy simulations
 *
 */

real rotcur_disk1(real r, int np, real *p, real *d)
{
  /* this #if 0 is needed for the 2.96 compiler on mdk81 :-) */
#if 0
  error("disk1 not implemented yet");
#endif
}


/******************************************************************************/
nemo_main()
{
    int  mask[PARAMS];   /* mask to define the free(1) or fixed(0) parameters */
    int  nring;  /* number of rings defined by users */
    int  side;   /* denotes which side of galaxy to be used */
    int  wpow;   /* denotes weigthing funtion to be used */
    int  cor[2];         /* plot error ellipses ? */
    real p[PARAMS],e[PARAMS]; /* arrays containing resp. pars and the errors */
    real ri,ro,r;     /* vars denoting inner, outer and mean radius */
    real x0,y0,vsys;  /* vars for init. estim. of xpos, ypos and vsys */
    real thf;     /* var  denoting free angle around minor axis */
    real nsigma;
    real old_factor, factor;    /* factor > 1, by which errors need be multiplied */
    int i,j,k;

    rotcurparse();
    if (nmod==0) error("No rotcur models specified");

    for (i=0; i<nrad; i++) 
      printf("%g  %g\n",rad[i],rotcur(rad[i]));
}

rotcurparse()
{
  string *sp, fname, path;
  char keyname[30];
  int i, j, nsp;
  bool Qext = hasvalue("load");
  char func_name[80];

  nrad = nemoinpr(getparam("radii"),rad,MAXRAD);

  if (Qext) {                          /* load= was used, load potentials from this file */
    fname = getparam("load");
    mysymbols(getargv0());
    path = pathfind(".",fname);
    if (path == NULL) error("Cannot open %s",fname);
    loadobj(path);
  }

  nmod = 0;                                 /* number of rotcur# keywords used */
  for (i=0; i<MAXMOD; i++) {                 /* process all rotcur#= keywords */
    sprintf(keyname,"rotcur%d",i+1);
    if (hasvalue(keyname)) {
      sp = burststring(getparam(keyname),", ");
      nsp = xstrlen(sp,sizeof(string))-2;  /* one for name, one for terminating 0 */

      if (streq(sp[0],"linear")) {             /* first check for predefined ones */
	if (nsp != 1) error("linear needs 1 parameter");
	npar[nmod] = 1;
	mpar[nmod][0] = natof(sp[1]);
	rcfn[nmod] = rotcur_linear;
      } else if (streq(sp[0],"flat")) {
	if (nsp != 1) error("flat needs 1 parameter");
	npar[nmod] = 1;
	mpar[nmod][0] = natof(sp[1]);
	rcfn[nmod] = rotcur_flat;
      } else if (streq(sp[0],"plummer")) {
	if (nsp != 2) error("plummer needs 2 parameters");
	npar[nmod] = 2;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	rcfn[nmod] = rotcur_plummer;
      } else if (streq(sp[0],"core1")) {
	if (nsp != 2) error("core1 needs 2 parameters");
	npar[nmod] = 2;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	rcfn[nmod] = rotcur_core1;
      } else if (streq(sp[0],"core2")) {
	if (nsp != 2) error("core2 needs 2 parameters");
	npar[nmod] = 2;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	rcfn[nmod] = rotcur_core2;
      } else if (streq(sp[0],"iso")) {
	if (nsp != 2) error("iso needs 2 parameters");
	npar[nmod] = 2;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	rcfn[nmod] = rotcur_iso;
      } else if (streq(sp[0],"power")) {
	if (nsp != 3) error("power needs 3 parameters");
	npar[nmod] = 3;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	mpar[nmod][2] = natof(sp[3]);
	rcfn[nmod] = rotcur_power;
      } else if (streq(sp[0],"brandt")) {
	if (nsp != 3) error("brandt needs 3 parameters");
	npar[nmod] = 3;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	mpar[nmod][2] = natof(sp[3]);
	rcfn[nmod] = rotcur_brandt;
      } else if (streq(sp[0],"exp")) {
	if (nsp != 2) error("exp needs 2 parameter");
	npar[nmod] = 2;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	rcfn[nmod] = rotcur_exp;
      } else if (streq(sp[0],"nfw")) {
	if (nsp != 3) error("nfw needs 3 parameters, although 3rd cannot be varied");
	npar[nmod] = 3;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	mpar[nmod][2] = natof(sp[3]);
	rcfn[nmod] = rotcur_nfw;
      } else if (streq(sp[0],"core")) {
	if (nsp != 3) error("core needs 3 parameters");
	npar[nmod] = 3;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	mpar[nmod][2] = natof(sp[3]);
	rcfn[nmod] = rotcur_core;
      } else if (streq(sp[0],"poly")) {
	if (nsp % 2) error("poly really needs an even number of numbers");
	npar[nmod] = nsp/2;
	for (j=0; j<npar[nmod]; j++) {
	  mpar[nmod][j] = natof(sp[j+1]);
	}
	rcfn[nmod] = rotcur_poly;
      } else {                /* else, if load= was used, try and find it there */
	if (Qext) {
	  sprintf(func_name,"rotcur_%s",sp[0]);
	  rcfn[nmod] = (rcproc) findfn(func_name);
	  if (rcfn[nmod]==NULL) error("Could not find %s in %s",func_name,fname);

	  if (nsp % 2) error("%s really needs an even number of numbers",sp[0]);
	  npar[nmod] = nsp/2;
	  for (j=0; j<npar[nmod]; j++)
	    mpar[nmod][j] = natof(sp[j+1]);
	} else
	  error("Cannot find %s, perhaps need load=",sp[0]);
      }
      dprintf(0,"ROTCUR%d: name=%s parameters=",i+1,sp[0]);
      for (j=0; j<npar[nmod]; j++) {
	dprintf(0,"%g ",mpar[nmod][j]);
      }
      dprintf(0,"\n");
      nmod++;
      freestrings(sp);
    } /* hasvalue */
  } /* i */
}



real rotcur(real r)
{
  int i;
  real vc, rcderv[MAXPAR];
  
  vc = 0;
  for (i=0; i<nmod; i++)
    vc += sqr((*rcfn[i])(r,npar[i],mpar[i],rcderv));
  return sqrt(vc);
}
