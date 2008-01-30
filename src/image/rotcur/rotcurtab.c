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
#include <rotcurshape.h>

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
    "VERSION=1.1\n   30-jan-08 PJT",
    NULL,
};

string usage="prints a table of a rotation curve";

string cvsid="$Id$";



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
