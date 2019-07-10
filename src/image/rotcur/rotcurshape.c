/*******************************************************************************
 *    ROTCURSHAPE derives the kinematical parameters from the observed
 *    velocity field by fitting tilted-rings to the velocity field.  It
 *    does a least-squares-fitting to the function:
 *
 *                v(x,y) = VSYS + VROT(r) * cos(theta) * sin(INC)
 *
 *                           - (x-XPOS) * sin(PA) + (y-YPOS) * cos(PA)
 *     with:    cos(theta) = -----------------------------------------
 *                                            r
 * 
 *     and a functional parameterized form for VROT(r)
 *
 *     History: 19/jul/02 : cloned off rotcur                                  pjt
 *              10-sep-02 : implemented resid= for images, more pts for disk   pjt
 *              19-sep-02 : added a few more rotcur's (exp, nfw)               pjt
 *              13-dec-02 : 1.0h: added power law rotation curve for Josh Simon  pjt
 *              30-jan-03 : 1.1: added error bars in V (unfinished)              pjt
 *              12-feb-03 : 1.1a: fixed residual field computation bug           pjt
 *              19-mar-03 : 1.2  clean up weights map dens= if beam= not given   pjt
 *               2-sep-03 : 1.2a  fix rotcurmode=t default masks                 pjt
 *              22-dec-03 : 1.2b  fix error estimates for the shape pars         pjt
 *              25-may-04 : 1.2c  fix sigma computation (sqrt(N))                pjt
 *                             e  plus many improvement when nsigma given
 *              11-jun-04 : 1.3   write fit directly , with or without masked    pjt
 *              13-jan-05 :    b  rachel mode: inc=0 -> cosi=sini=1.0 to get circles    pjt
 *              26-may-05 :    d  fixed bad NFW bug (used v^2, not v)                   pjt
 *              30-jan-08   1.4   shapes in rotcurs.c in library now                    pjt
 *
 *
 ******************************************************************************/

#include <stdinc.h>
#include <getparam.h>
#include <image.h>
#include <filefn.h>
#include <extstring.h>
#include <loadobj.h>
#include <rotcurshape.h>

/*     Set this appropriately if you want to use NumRec's mrqmin() based engine */
/*     right now it appears as if the NR routine does not work as well          */
/*     See also my tabnllsqfit experiments where you can test both side by side */
#if 0
#define nllsqfit nr_nllsqfit
#endif

#define GPARAM  5           /* number of geometric parameters for a disk */
#define PARAMS  50          /* (max) number of parameters */
#define MAXMOD  5           /* number of rotation curve models (should match rotcur#) */
#define MAXPAR  5           /* number of parameters per model (can be made larger) */

#define RING         10     /* maximum number of rings (17 arrays) */
#define MAXPTS   100000     /* maximum number of pixels per ring (9 arrays) - only for tabular info*/

#define DEF_TOL   0.001     /* tolerance for fit */
#define DEF_LAB   0.001     /* mixing parameter */
#define DEF_T     50        /* maximum number of iterations for fit */

#define F 0.0174532925  /* deg to rad conversion     a.k.a. pi / 180 */
#define G 0.4246609001  /* FWHM to sigma conversion, a.k.a. 1 / 2sqrt(2ln2))  */

string defv[] = {
    "in=???\n        Input image velocity field",
    "radii=\n        Radii of rings (arcsec)",
    "pa=\n           Position angle (degrees)",
    "inc=\n          Inclination (degrees)",
    "vsys=\n         Systemic velocity",
    "center=\n       Rotation center (grids w.r.t. 0,0) [center of map]",
    "frang=0\n       Free angle around minor axis (degrees)",
    "side=\n         Side to fit: receding, approaching or [both]",
    "weight=u\n      Weighting function: {[uniform],cosine,cos-squared}",
    "fixed=\n        Geometric parameters to be kept fixed {vsys,xpos,ypos,pa,inc}",
    "ellips=\n       ** Parameters for which to plot error ellips",
    "beam=\n         ** Beam (arcsec) for beam correction [no correction]",
    "dens=\n         Image containing containing density map to be used as weight",
    "tab=\n          If specified, this output table is used in append mode",
    "resid=\n        Output of residual field",
    "fit=f\n         Output the fit? or the residuals",
#if 0
    "masked=t\n      Masked fitted field or Full fitted field IFF fit=true",
#endif
    "tol=0.001\n     Tolerance for convergence of nllsqfit",
    "lab=0.001\n     Mixing parameter for nllsqfit",
    "itmax=50\n      Maximum number of allowed nllsqfit iterations",
    "units=deg,1\n   Units of input {deg, arcmin, arcsec, rad, #},{#} for length and velocity",
    "blank=0.0\n     Value of the blank (pixel) value to be ignored",
    "nsigma=-1\n     Iterate once by rejecting points more than nsigma resid",
    "imagemode=t\n   Input image mode? (false means ascii table)",
    "rotcurmode=f\n  Full velocity field, or rotcur (r,v) fit only",
    "load=\n         dynamically loadobject file with rotcur_<NAME>",
    "rotcur1=\n      Rotation curve <NAME>, parameters and set of free(1)/fixed(0) values",
    "rotcur2=\n      Rotation curve <NAME>, parameters and set of free(1)/fixed(0) values",
    "rotcur3=\n      Rotation curve <NAME>, parameters and set of free(1)/fixed(0) values",
    "rotcur4=\n      Rotation curve <NAME>, parameters and set of free(1)/fixed(0) values",
    "rotcur5=\n      Rotation curve <NAME>, parameters and set of free(1)/fixed(0) values",
    "VERSION=1.4c\n  10-jul-2019 PJT",
    NULL,
};

string usage="nonlinear fit of kinematical parameters to the velocity field of a coplanar disk";

string cvsid="$Id$";



imageptr denptr, velptr, resptr;      /* pointers to the various Images -- some recycling */
int   lmin, mmin, lmax, mmax;              /* boundaries of map */
real  grid[2];    /* grid separations in x and y (arcsec.) */
real  beam[2];    /* size of beam in arcseconds            */
real  dx,dy;      /* grid separation in x and y */
real  undf;       /* undefined value in map */
real  pamp;       /* position angle of map */

typedef struct {
  int    npar;
  int    ipar;
  real   mpar[MAXPAR];
  int    mmsk[MAXPAR];
  rcproc rcfn;
} rotcurz[MAXMOD];

int    npar[MAXMOD];            /* active number of parameters per model */
int    ipar[MAXMOD];
real   mpar[MAXMOD][MAXPAR];    /* parameters per model */
int    mmsk[MAXMOD][MAXPAR];    /* mask per model (1=free 0=fixed) */
rcproc rcfn[MAXMOD];

int    nmod = 0;                /* number of models actually used */
int    nparams = 0;             /* total number of parameters (5 + # for models) */


bool Qimage;                             /* input mode (false means tables are used) */
bool Qrotcur;                            /* rotcur (rv) vs. velocity field (xyv) table mode */
bool Qbeam;                              /* attempt to do beam correction */
bool Qfit;                               /* write fit instead of residual */
bool Qmasked;                            /* in fitted output image, use observation mask ? */
bool Qrachel = FALSE;                    /* special mode where inc=0 and sini=cosi=1 */
real  *xpos_vel, *ypos_vel, *vrad_vel, *verr_vel;   /* pointer to tabular information */
real  *vsig_vel;                         
int  n_vel = 0;                          /* length of tabular arrays */

real tol,lab;     /* parameters that go into nllsqfit */
int  itmax;

/* Compute engine, derivative and beam smearing correction functions:
 * _c1 = cos(theta)	
 * _s1 = sin(theta)
 * _c2 = cos(2*theta)
 * _s2 = sin(2*theta)  etc.
 */

real vobs_c1 (real *c, real *p, int m);
void vobsd_c1(real *c, real *p, real *d, int m);
void vcor_c1 (real *c, real *p, real *vd, real *dn);
real vobs_s1 (real *c, real *p, int m);
void vobsd_s1(real *c, real *p, real *d, int m);
void vcor_s1 (real *c, real *p, real *vd, real *dn);


int rotinp(real *rad, real pan[], real inc[], real vro[], int *nring, int ring, real *vsys, 
	   real *x0, real *y0, real *thf, int *wpow, int mask[], int *side, int cor[], 
	   real *nsigma, stream lunpri);
int rotfit(real ri, real ro, real p[], real e[], int mask[], int wpow, int side, real thf, 
	   real elp4[], int cor[], int *npt, real nsigma, stream lunres);
int perform_out(int h, real p[6], int n, real q);
void rotplt(real rad[], real vsy[], real evs[], real pan[], real epa[], 
	   real inc[], real ein[], real xce[], real exc[], real yce[], real eyc[], 
	   real p[], real e[],
	   int mask[], int ifit, real elp[][4], stream lunpri, int cor[], int npt[], real factor);
int set_blank(int n, real *res, real *w, int *iblank, int nfr, real nsigma, real *q);
void stat2(real a[], int n, int nf, real *mean, real *sig);
void stat2b(real a[], int b[], int n, int nf, real *mean, real *sig);
void getdat(real x[], real y[], real w[], int idx[], real res[], int *n, int nmax, 
	   real p[], real ri, real ro, real thf, int wpow, real *q, int side, bool *full, int nfr, int mode);
real bmcorr(real xx[2], real p[], int l, int m);
int perform_init(real *p, real *c);


typedef real (*my_proc1)(real *, real *, int);
typedef void (*my_proc2)(real *, real *, real *, int);
typedef void (*my_proc3)(real *, real *, real *, real *);
 
extern int nllsqfit(real *xdat, int xdim, real *ydat, real *wdat, real *ddat, int ndat, real *fpar, 
		    real *epar, int *mpar, int npar, real tol, int its, real lab, my_proc1 f, my_proc2 df);


my_proc1 vobs;
my_proc2 vobsd;
my_proc3 vcor;			/* pointers to the correct functions */

extern string *burststring(string,string);
extern void freestrings(string *);
extern bool scanopt(string option, string key);

extern int match(string, string, int *);


/******************************************************************************/
void nemo_main()
{
    int  ier;            /* error return code */
    int  ifit=0;         /* counter for number of succesful fits */
    int  n;		 /* number of points in a ring */
    int  irng;           /* loop-counter */
    stream lunpri;       /* file for table output */
    stream lunres;       /* file for residual output */
    int  mask[PARAMS];   /* mask to define the free(1) or fixed(0) parameters */
    int  nring;  /* number of rings defined by users */
    int  side;   /* denotes which side of galaxy to be used */
    int  wpow;   /* denotes weigthing funtion to be used */
    int  cor[2];         /* plot error ellipses ? */
    real elp[RING][4],elp4[4];  /* array containing ellipse parameters */
    real rad[RING+1];         /* array contains radii of rings */
    real vro[RING],evr[RING]; /* arrays containing resp. vrot and its error */
    real vsy[RING],evs[RING]; /* arrays containing resp. vsys and its error */
    real pan[RING],epa[RING]; /* arrays containing resp. p.a. and its error */
    real inc[RING],ein[RING]; /* arrays containing resp. inc. and its error */
    real xce[RING],exc[RING]; /* arrays containing resp. xpos and its error */
    real yce[RING],eyc[RING]; /* arrays containing resp. ypos and its error */
    int  npt[RING];	      /* array containing number of points in ring */
    real p[PARAMS],e[PARAMS]; /* arrays containing resp. pars and the errors */
    real ri,ro,r;     /* vars denoting inner, outer and mean radius */
    real x0,y0,vsys;  /* vars for init. estim. of xpos, ypos and vsys */
    real thf;     /* var  denoting free angle around minor axis */
    real nsigma;
    real old_factor, factor;    /* factor > 1, by which errors need be multiplied */
    int i,j,k;

    rotcurparse();
    if (nmod==0) error("No rotcur models specified");

    for (i=0, nparams=GPARAM; i<nmod; i++)     /* count total number of parameters */
      nparams += npar[i];

    if (hasvalue("tab"))
        lunpri=stropen(getparam("tab"),"a");  /* pointer to table stream output */
    else
        lunpri=NULL;                    /* no table output */
    Qimage = getbparam("imagemode");
    Qrotcur = getbparam("rotcurmode");
    if (Qrotcur) Qimage=FALSE;
    Qfit = getbparam("fit");
#if 0
    Qmasked = getbparam("masked");
#endif

    if (hasvalue("resid")) {
      if (Qimage) {
	lunres=stropen(getparam("resid"),"w");  /* pointer to table stream output */
        resptr = NULL;                          /* signal that we need to allocate */
      } else
        lunres=stropen(getparam("resid"),"a");  /* pointer to table stream output */
    } else
        lunres=NULL;                    /* no residual table or image output */

    rotinp(rad,pan,inc,vro,&nring,RING,     /* get input parameters */
           &vsys,&x0,&y0,&thf,
           &wpow,mask,&side,cor,&nsigma,lunpri);

    old_factor = sqrt((beam[0]+grid[0])*(beam[1]+grid[1])/(grid[0]*grid[1]));
    if (beam[0] > 0 && beam[1] > 0)
      factor = sqrt(FOUR_PI*beam[0]*beam[1]/(grid[0]*grid[1])); /* Sicking 1997 !!! */
    else
      factor = 1.0;
    dprintf(1,"Sicking (1997)'s error multiplication factor=%g  (old_factor=%g)\n",
	    factor,old_factor);

    for (irng=0; irng<nring-1; irng++) {  /* loop for each ring */
         ri=rad[irng];          /* inner radius of ring */
         ro=rad[irng+1];        /* outer radius of ring */
         r=0.5*(ri+ro);         /* mean radius of ring */

         p[0] = vsys;
         p[1] = x0;
         p[2] = y0;
         p[3] = pan[irng];
         p[4] = inc[irng];
	 for (i=0, k=0; i<nmod; i++) {
	   for (j=0; j<npar[i]; j++, k++) {
	     p[GPARAM+k] = mpar[i][j];
	     mask[GPARAM+k] = mmsk[i][j];
	   }
	 }

         ier = rotfit(ri,ro,p,e,mask,wpow,side,thf,elp4,cor,&n,-1.0,lunres);
	 if (ier > 0 && nsigma > 0)
	   ier = rotfit(ri,ro,p,e,mask,wpow,side,thf,elp4,cor,&n,nsigma,lunres);
	 else if (nsigma > 0 && ier < 0)
	   warning("bad fit; ier=%d; 2nd improved residual image not written",ier);
         if (ier>0) {           /* only if fit OK, store fit */
            rad[ifit]=r;            /*  radius of ring */
            vsy[ifit]=p[0];         /*  systemic velocity */
            evs[ifit]=e[0]*factor;  /*  error in systemic velocity */
            xce[ifit]=p[1];         /*  x-position */
            exc[ifit]=e[1]*factor;  /*  error in x-position */
            yce[ifit]=p[2];         /*  y-position */
            eyc[ifit]=e[2]*factor;  /*  error in y-position */
            pan[ifit]=p[3];         /*  position angle */
            epa[ifit]=e[3]*factor;  /*  error in position angle */
            inc[ifit]=p[4];         /*  inclination */
            ein[ifit]=e[4]*factor;  /*  error in inclination */
            vro[ifit]=p[5];         /*  circular velocity */
            evr[ifit]=e[5]*factor;  /*  error in circular velocity */
            elp[ifit][0]=elp4[0];   /*  save ellipse parameters */
            elp[ifit][1]=elp4[1];   /* NOT corrected by 'factor' */
            elp[ifit][2]=elp4[2];
            elp[ifit][3]=elp4[3];
	    npt[ifit] = n;
            ifit++;
         }
    } /* end of loop through rings */

    rotplt(rad,vsy,evs,pan,epa,         /* output the results */
           inc,ein,xce,exc,yce,eyc,p,e,
           mask,ifit,elp,lunpri,cor,npt,factor);
}

rotcurparse()
{
  string *sp, fname, path;
  char keyname[30];
  int i, j, nsp;
  bool Qext = hasvalue("load");
  char func_name[80];

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
      if (nsp % 2) warning("%s= should use an even number of numbers",keyname);

      if (streq(sp[0],"linear")) {             /* first check for predefined ones */
	if (nsp != 2) error("linear needs 1 parameter");
	npar[nmod] = 1;
	mpar[nmod][0] = natof(sp[1]);
	mmsk[nmod][0] = natoi(sp[2]);
	rcfn[nmod] = rotcur_linear;
      } else if (streq(sp[0],"flat")) {
	if (nsp != 2) error("flat needs 1 parameter");
	npar[nmod] = 1;
	mpar[nmod][0] = natof(sp[1]);
	mmsk[nmod][0] = natoi(sp[2]);
	rcfn[nmod] = rotcur_flat;
      } else if (streq(sp[0],"plummer")) {
	if (nsp != 4) error("plummer needs 2 parameters");
	npar[nmod] = 2;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	mmsk[nmod][0] = natoi(sp[3]);
	mmsk[nmod][1] = natoi(sp[4]);
	rcfn[nmod] = rotcur_plummer;
      } else if (streq(sp[0],"core1")) {
	if (nsp != 4) error("core1 needs 2 parameters");
	npar[nmod] = 2;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	mmsk[nmod][0] = natoi(sp[3]);
	mmsk[nmod][1] = natoi(sp[4]);
	rcfn[nmod] = rotcur_core1;
      } else if (streq(sp[0],"core2")) {
	if (nsp != 4) error("core2 needs 2 parameters");
	npar[nmod] = 2;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	mmsk[nmod][0] = natoi(sp[3]);
	mmsk[nmod][1] = natoi(sp[4]);
	rcfn[nmod] = rotcur_core2;
      } else if (streq(sp[0],"iso")) {
	if (nsp != 4) error("iso needs 2 parameters");
	npar[nmod] = 2;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	mmsk[nmod][0] = natoi(sp[3]);
	mmsk[nmod][1] = natoi(sp[4]);
	rcfn[nmod] = rotcur_iso;
      } else if (streq(sp[0],"power")) {
	if (nsp != 6) error("power needs 3 parameters");
	npar[nmod] = 3;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	mpar[nmod][2] = natof(sp[3]);
	mmsk[nmod][0] = natoi(sp[4]);
	mmsk[nmod][1] = natoi(sp[5]);
	mmsk[nmod][2] = natoi(sp[6]);
	rcfn[nmod] = rotcur_power;
      } else if (streq(sp[0],"brandt")) {
	if (nsp != 6) error("brandt needs 3 parameters");
	npar[nmod] = 3;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	mpar[nmod][2] = natof(sp[3]);
	mmsk[nmod][0] = natoi(sp[4]);
	mmsk[nmod][1] = natoi(sp[5]);
	mmsk[nmod][2] = natoi(sp[6]);
	rcfn[nmod] = rotcur_brandt;
      } else if (streq(sp[0],"exp")) {
	if (nsp != 4) error("exp needs 2 parameter");
	npar[nmod] = 2;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	mmsk[nmod][0] = natoi(sp[3]);
	mmsk[nmod][1] = natoi(sp[4]);
	rcfn[nmod] = rotcur_exp;
      } else if (streq(sp[0],"nfw")) {
	if (nsp != 6) error("nfw needs 3 parameters, although 3rd cannot be varied");
	npar[nmod] = 3;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	mpar[nmod][2] = natof(sp[3]);
	mmsk[nmod][0] = natoi(sp[4]);
	mmsk[nmod][1] = natoi(sp[5]);
	mmsk[nmod][2] = natoi(sp[6]);
	rcfn[nmod] = rotcur_nfw;
      } else if (streq(sp[0],"core") || streq(sp[0],"rotcurm")) {
	if (nsp != 6) error("core/rotcurm needs 3 parameters");
	npar[nmod] = 3;
	mpar[nmod][0] = natof(sp[1]);
	mpar[nmod][1] = natof(sp[2]);
	mpar[nmod][2] = natof(sp[3]);
	mmsk[nmod][0] = natoi(sp[4]);
	mmsk[nmod][1] = natoi(sp[5]);
	mmsk[nmod][2] = natoi(sp[6]);
	rcfn[nmod] = rotcur_core;
      } else if (streq(sp[0],"poly")) {
	if (nsp % 2) error("poly really needs an even number of numbers");
	npar[nmod] = nsp/2;
	for (j=0; j<npar[nmod]; j++) {
	  mpar[nmod][j] = natof(sp[j+1]);
	  mmsk[nmod][j] = natoi(sp[j+1+npar[nmod]]);
	}
	if (mmsk[nmod][0] && mmsk[nmod][1]) 
	  error("Polynomial needs at least one of P1 and P2 fixed");
	rcfn[nmod] = rotcur_poly;
      } else {                /* else, if load= was used, try and find it there */
	if (Qext) {
	  sprintf(func_name,"rotcur_%s",sp[0]);
	  rcfn[nmod] = (rcproc) findfn(func_name);
	  if (rcfn[nmod]==NULL) error("Could not find %s in %s",func_name,fname);

	  if (nsp % 2) error("%s really needs an even number of numbers",sp[0]);
	  npar[nmod] = nsp/2;
	  for (j=0; j<npar[nmod]; j++) {
	    mpar[nmod][j] = natof(sp[j+1]);
	    mmsk[nmod][j] = natoi(sp[j+1+npar[nmod]]);
	  }
	} else
	  error("Cannot find %s, perhaps need load=",sp[0]);
      }
      dprintf(0,"ROTCUR%d: name=%s parameters=",i+1,sp[0]);
      for (j=0; j<npar[nmod]; j++) {
	dprintf(0,"%g (%s) ",mpar[nmod][j], mmsk[nmod][j] ? "free" : "fixed");
      }
      dprintf(0,"\n");
      nmod++;
      freestrings(sp);
    } /* hasvalue */
  } /* i */
}

/*
 *    ROTINP: This function inputs parameters from the commandline
 *
 *    RAD      real array       radii of rings
 *    PAN      real array       position angles of rings
 *    INC      real array       inclinations of rings
 *    VRO      real array       circular velocities of rings
 *    NRING    integer          number of rings
 *    RING     integer          maximum number of rings
 *    VSYS     real             systemic velocity
 *    X0       real             x-position of centre
 *    Y0       real             y-position of centre
 *    THF      real             free angle around minor axis
 *    WPOW     integer          weighting mode
 *    MASK     integer array    keep parameters fixed(0) or free(1)
 *    SIDE     integer          receding, approaching or both sides
 *    COR      integer array    plot error ellipses ?
 *    LUNPRI   integer          LUN for print output
 */

int rotinp(rad,pan,inc,vro,nring,ring,vsys,x0,y0,thf,wpow,mask,side,cor,nsigma,lunpri)
int  *wpow;           /* weighting function to be applied */
int  ring, *nring;    /* max. number of rings and wanted number of rings */
int  mask[];          /* mask to define free and fixed parameters */
int  *side;           /* variable denotes which part of galaxy to fit */
int  cor[];           /* plot error ellipses ? */
real *rad;            /* user defined radii of rings */
real *vsys,vro[],pan[],inc[],*x0,*y0;  /* initial estimates (pars. 1 <--> PARAMS ) */
real *thf;            /* user defined free angle around minor axis */
real *nsigma;
stream  lunpri;       /* LUN for print output */
{
    char *input;
    string *inputs;
    int iret, i, j, k, n, nfixed, fixed, ninputs;
    real center[2], toarcsec, tokms;
    stream velstr, denstr;
    bool Qdens;
    real *coldat[3];
    int colnr[3];

    dprintf(0,"%s: NEMO VERSION %s\n", 
                        getparam("argv0"), getparam("VERSION"));
    if (lunpri) fprintf(lunpri,"%s: VERSION %s [NEMO]\n\n",
                        getparam("argv0"), getparam("VERSION"));

    for (i=0; i<PARAMS; i++) mask[i] = 1;       /* default: set all free */

    input = getparam("in");
    if (lunpri) fprintf(lunpri," file                : %s\n",input);
    if (lunpri) fprintf(lunpri," velocity field file : %s (%s)\n",input,
			Qimage ? "image" : "ascii table");

    velstr = stropen(input,"r");                /* open velocity field */  
    Qdens = hasvalue("dens");
    if (Qimage) {
      read_image(velstr,&velptr);                 /* get data */
    } else {
      Qdens = getbparam("dens");              /* redo this, since now it's t/f */
      n_vel = nemo_file_lines(input,100000);
      xpos_vel = (real *) allocate(n_vel * sizeof(real));
      ypos_vel = (real *) allocate(n_vel * sizeof(real));
      vrad_vel = (real *) allocate(n_vel * sizeof(real));
      vsig_vel = (real *) allocate(n_vel * sizeof(real));
      if (Qdens) verr_vel = (real *) allocate(n_vel * sizeof(real));
      if (Qrotcur) {                         /* 1dim rotation curve */
	colnr[0] = 1;    coldat[0] = ypos_vel;
	colnr[1] = 2;    coldat[1] = vrad_vel;
	if (Qdens) {
	  colnr[2] = 3;  coldat[2] = verr_vel;
	  n_vel = get_atable(velstr,3,colnr,coldat,n_vel);
	  dprintf(0,"[Found %d points in rotation curve (R,V,DV)]\n",n_vel);
	} else {
	  n_vel = get_atable(velstr,2,colnr,coldat,n_vel);
	  verr_vel = NULL;
	  dprintf(0,"[Found %d points in rotation curve (R,V)]\n",n_vel);
	}
	for (i=0; i<n_vel; i++) {
	  xpos_vel[i] = 0.0;
	  vsig_vel[i] = 1.0;  /* not used yet */
	}
      } else {                            /* 2dim velocity field */
	colnr[0] = 1;    coldat[0] = xpos_vel;
	colnr[1] = 2;    coldat[1] = ypos_vel;
	colnr[2] = 3;    coldat[2] = vrad_vel;
	if (Qdens) {
	  colnr[3] = 4;  coldat[3] = verr_vel;
	  n_vel = get_atable(velstr,4,colnr,coldat,n_vel);
	  dprintf(0,"[Found %d points in velocity field table (X,Y,V,DV)]\n",n_vel);
	} else {
	  n_vel = get_atable(velstr,3,colnr,coldat,n_vel);
	  verr_vel = NULL;
	  dprintf(0,"[Found %d points in velocity field table (X,Y,V)]\n",n_vel);
	}
	for (i=0; i<n_vel; i++) {
	  vsig_vel[i] = 1.0;   /* not used yet */
	}      
      }
      for (i=0; i<n_vel; i++) {
	dprintf(8,"%g %g %g %g\n",xpos_vel[i],ypos_vel[i],vrad_vel[i],vsig_vel[i]);
      }
      velptr = NULL;
    }
    strclose(velstr);                           /* and close file */
    
    inputs = burststring(getparam("units"),",");
    ninputs = xstrlen(inputs,sizeof(string))-1;
    if (ninputs > 0) {
      if (streq(inputs[0],"deg"))
	toarcsec = 3600.0;
      else if (streq(inputs[0],"arcmin") || streq(inputs[0],"min"))
	toarcsec = 60.0;
      else if (streq(inputs[0],"arcsec") || streq(inputs[0],"sec"))
	toarcsec = 1.0;
      else if (streq(inputs[0],"rad"))
	toarcsec = 3600.0 * 180/PI;
      else {
	toarcsec = natof(inputs[0]);
	printf("Conversion factor %g used to get arcsec\n",toarcsec);
      }
    } else if (ninputs == 0)
      toarcsec = 1.0;
    else if (ninputs < 0)
      error("Bad units");

    if (ninputs > 1)
      tokms = natof(inputs[1]);
    else
      tokms = 1.0;

    if (velptr) {
      lmin=0;                         /* coordinates of map (xlo) */
      lmax=Nx(velptr)-1;              /*                    (xhi) */
      mmin=0;                         /*                    (ylo) */
      mmax=Ny(velptr)-1;              /*                    (yhi) */
      dx=toarcsec*Dx(velptr);         /* separation in X (in arcsec.) */
      dy=toarcsec*Dy(velptr);         /* separation in Y (in arcsec.) */
      if (dx<0) {
        warning("Repairing negative dx");
        dx = -dx;
      }
      if (dy<0) {
        warning("Repairing negative dy");
        dy = -dy;
      }

      undf=getdparam("blank");        /* the undefined value */
      n = 0;                          /* count # blank values in velmap */
      for (j=0; j<Ny(velptr); j++)
        for (i=0; i<Nx(velptr); i++) 
	  if (MapValue(velptr,i,j)==undf)
	    n++;
	  else
	    MapValue(velptr,i,j) *= tokms;
      
      printf("Mapsize is %g * %g arcsec , pixel=%g*%g; Found %d/%d undefined map values\n",
	     ABS(dx*(lmax-lmin+1.0)), ABS(dy*(mmax-mmin+1.0)), 
	     dx,dy,
	     n, Nx(velptr)*Ny(velptr));
      
      grid[0]=dx;        /* grid-separations in VELPAR (dx) */
      grid[1]=dy;        /*                            (dy) */
      pamp=0.0;          /* position angle of map */
    } else {
      undf=getdparam("blank");        /* the undefined value */
      grid[0]=dx=1.0;
      grid[1]=dy=1.0;
      pamp=0.0;
      printf("Mapsize unkown for point data, but toarcsec = %g, dx=%g\n",
	     toarcsec,dx);
    }

    /*** WARNING: beam= and dens= are currently tied, if beam= used,
	 the dens= is **also** used for beam smearing corrections,
	 which is a largely untested feature of the code 
    ***/
    
    n = nemoinpr(getparam("beam"),beam,2);   /* get size of beam from user */
    if (n==2 || n==1) {
      Qbeam = TRUE;
      printf("With beam correction\n");
      if (n==1) beam[1] = beam[0];               /* make circular beam */
      if (lunpri) fprintf(lunpri,"  beam: %g %g\n",beam[0],beam[1]);
    } else {
      Qbeam = FALSE;
      printf("No beam correction\n");
      beam[1] = beam[0] = 0.0;
      if (n!=0) warning("Parsing error beam=%s",getparam("beam"));
    }
    if (Qimage) {
      if (Qdens) {
	input = getparam("dens");
	if (lunpri) fprintf(lunpri," density file        : %s  beam: %g %g\n",
			    input,beam[0],beam[1]);	
	denstr = stropen(input,"r");
	read_image(denstr,&denptr);
	strclose(denstr);
	warning("Using density map for weights now");
      } else {
	denptr = NULL;
      }
    }

    *nring = nemoinpr(getparam("radii"),rad,ring+1);
    if (*nring != 2) error("radii=: Need two radii for a disk");
    *vsys = getdparam("vsys");
    n = nemoinpr(getparam("pa"),pan,ring);
    if (n<1) error("pa=: need at least one position angle (%d)",n);
    for (i=n;i<*nring;i++)
      pan[i] = pan[n-1];
    n = nemoinpr(getparam("inc"),inc,ring);
    if (n<1) error("inc=: need at least one inclincation (%d)",n);
    for (i=n;i<*nring;i++)
      inc[i] = inc[n-1];
    if (inc[0] == 0.0) {
      warning("Special Rachel mode to try and get circles, better have fixed=inc also");
      Qrachel = TRUE;
    }
    n = nemoinpr(getparam("center"),center,2);
    if (n==2) {                     /* if two value supplied */
      *x0 = center[0];            /* this will be the center of rotation */
      *y0 = center[1];
    } else if (n==1 && Qrotcur) {   /* rotcurmode: slit forced along Y */
      *x0 = 0;
      *y0 = center[0];            /* this will be the center of rotation */
    } else if (n==0) {                   /* if nothing supplied: */
      *x0 = 0.5*(lmin+lmax);      /* use center of map */
      *y0 = 0.5*(mmin+mmax);
    } else                          /* if all fails - barf */
      error("Need two numbers to define the center (%d)",n);
    *thf = getdparam("frang");

    printf("ROTCURSHAPE: free angle %4.1f (degrees)\n", *thf);
    if (lunpri) fprintf(lunpri," free angle          : %4.1f (degrees)\n",*thf);

    input = getparam("side");
    if (*input=='r') {
        *side = 1;
        printf("ROTCURSHAPE: rotation curve of receding half\n");
        if (lunpri) fprintf(lunpri," velocity field      : receding half\n");
    } else if (*input=='a') {
        printf("ROTCURSHAPE: rotation curve of approaching half\n");
        if (lunpri) fprintf(lunpri," velocity field      : approaching side\n");
        *side = 2;
    } else if (*input=='w') {
        error("side=wedge not supported yet");
    } else {
        printf("ROTCURSHAPE: rotation curve of both halves\n");
        if (lunpri) fprintf(lunpri," velocity field      : both halves\n");
        *side = 3;
    }
    if (*side != 3) {  /* should we keep some parameters fixed ? */
         mask[0]=0;         /* fixed: systemic velocity */
         mask[1]=0;         /* fixed: x-position of center */
         mask[2]=0;         /* fixed: y-position of center */
    }

    input = getparam("weight");
    if (*input == 'u') {
        printf("ROTCURSHAPE: uniform weighting\n");
        if (lunpri) fprintf(lunpri," used weighting      : uniform\n");
        *wpow = 0;
    } else if (streq(input,"cos-squared")) {
        printf("ROTCURSHAPE: cosine-squared weighting\n");
        if (lunpri) fprintf(lunpri," used weighting      : cos-squared\n");
        *wpow = 2;
    } else {            /* default is cosine weighting */
        printf("ROTCURSHAPE: cosine weighting\n");
        if (lunpri) fprintf(lunpri," used weighting      : cosine\n");
        *wpow = 1;
    }

    fixed = 0;
    iret=match(getparam("fixed"),"vsys,xpos,ypos,pa,inc,all",&fixed);
    if (iret<0) error("Illegal option in fixed=%s",getparam("fixed"));
    dprintf(1,"MASK: 0x%x ",fixed);
    if (Qrotcur) {                  /* for rotcurmode need to fix xpos,pa,inc */
      if (mask[1]==1) {
	warning("rotcurmode: XPOS fixed at 0");
	*x0 = 0.0;
	mask[1]=0;
      }
      if (mask[3]==1) {
	warning("rotcurmode: PA fixed at 0");
	pan[0] = 0.0;
	mask[3]=0;
      }
      if (mask[4]==1) {
	warning("rotcurmode: INC fixed at 30 (NOTE:: vel/=2 !!)");
	inc[0] = 30.0;
	mask[4]=0;
      }
    } 
    /* BUG: some of this is messing up some of the Qrotcur cases */
    /*      solution is to use fixed= explicitly via cmdline */
    if (fixed & (1<<GPARAM)) {
      for (i=0; i<GPARAM; i++) mask[i] = 0;
    } else {
      for (i=0; i<GPARAM; i++) {
	mask[i] = (fixed & (1<<i)) ? 0 : 1;
      }
    }
    for (i=0; i<GPARAM; i++)
      dprintf(1,"%d ",mask[i]);
    dprintf(1,"\n");

    for (i=0, k=0; i<nmod; i++) {
      for (j=0; j<npar[i]; j++, k++)
	mask[k+GPARAM] = mmsk[i][j];
    }

    for (nfixed=0,i=0; i<nparams; i++)      /* count number of fixed par's */
        if (mask[i]==0) nfixed++; 

    printf("ROTCURSHAPE: will fit the following parameter(s)\n");
    if (lunpri) fprintf(lunpri," parameters to fit   :");
    if (mask[0]==1) {
        printf("       - Systemic velocity \n");
        if (lunpri) fprintf(lunpri," vsys");
    }
    if (mask[1]==1) {
        printf("       - X position of center \n");
        if (lunpri) fprintf(lunpri," xpos");
    }
    if (mask[2]==1) {
        printf("       - Y position of center \n");
        if (lunpri) fprintf(lunpri," ypos");
    }
    if (mask[3]==1) {
        printf("       - Position angle \n");
        if (lunpri) fprintf(lunpri," pa");
    }
    if (mask[4]==1) {
        printf("       - Inclination \n");
        if (lunpri) fprintf(lunpri," inc");
    }
    for (i=GPARAM; i<nparams; i++) {
      if (mask[i]==1) {
	printf("       - P%d \n",i-GPARAM+1);
	if (lunpri) fprintf(lunpri," P%d",i-GPARAM+1);
      }
    }
    if (lunpri) fprintf(lunpri,"\n");

    /*******  FIX THIS OLD ROTCUR CODE **********/

    input = getparam("ellips");
    n=0;
    cor[0] = cor[1] = 0;        /* default: no ellipses */
    if (scanopt(input,"vsys")) if (n++ < 2) cor[n-1] = 1;
    if (scanopt(input,"vrot")) if (n++ < 2) cor[n-1] = 2;
    if (scanopt(input,"pa"))   if (n++ < 2) cor[n-1] = 3;
    if (scanopt(input,"inc"))  if (n++ < 2) cor[n-1] = 4;
    if (scanopt(input,"xpos")) if (n++ < 2) cor[n-1] = 5;
    if (scanopt(input,"ypos")) if (n++ < 2) cor[n-1] = 6;
    if (n>2) warning("can only plot two-dimensional ellips");

    tol = getdparam("tol");
    lab = getdparam("lab");
    itmax = getiparam("itmax");
    printf("ROTCURSHAPE: tol=%g lab=%g iTmax=%d\n",tol,lab,itmax);
    if (lunpri) fprintf(lunpri," tol, lab, itmax     : %g %g %d\n",
                                    tol,lab,itmax);

    if (lunpri) fprintf(lunpri,"\n\n\n\n"); /* space space space... */

    vobs  = vobs_c1;    
    vobsd = vobsd_c1;   
    vcor  = vcor_c1;    

    *nsigma = getdparam("nsigma");
}

/*
 *    ROTFIT does a  least  squares fit to the radial velocity
 *    field.
 *
 *    RI       real            inner radius of ring
 *    RO       real            outer radius of ring
 *    P        real array      estimated/fitted parameters
 *    E        real array      errors in parameters
 *    MASK     integer array   contains mask (fixed or free parameter)
 *    WPOW     integer         which kind of weight
 *    SIDE     integer         which part of galaxy
 *    THT      real            free angle around minor axis
 *    COR      integer         plot error ellipses ?
 *    ELP4     real array      for ellipse parameters
 *
 *  returns:
 *    IER      integer         result of fit (good or bad)
 *    NPT      integer         number of points in ring
 */

int rotfit(ri, ro, p, e, mask, wpow, side, thf, elp4, cor, npt, nsigma, lunres)
real ri,ro;      /* inner and outer radius of ring */
int mask[];      /* mask for free/fixed parameters */
int wpow;        /* weighting function */
int side;        /* half of galaxy */
int cor[];       /* print error ellipses ? */
real elp4[];     /* array for ellipse parameters */
real p[],e[];    /* initial estimates (I), results and erros (O) */
real thf;        /* free angle around minor axis (I) */
int *npt;	 /* number of points in ring (0) */
real nsigma;     /* if positive, remove outliers and fit again, also iff lunres */
stream lunres;   /* file for residuals */
{
    int ier;                                             /* error return code */
    bool  stop,full;           /* booleans for stop fitting and data overflow */
    int   nfr;                                   /* number of free parameters */
    int   h, k;                                              /* loop counters */
    int   nrt;                       /* error return code from subroutine FIT */
    real  pf[PARAMS];               /* array for storing intermediate results */
    real  df[PARAMS];                 /* array which stores difference vector */
    real  b[PARAMS];                         /* array for partial derivatives */
    real  eps[PARAMS];                             /* contains stop criterium */
    real  flip;                               /* factor for difference vector */
    real  chi,q;                                   /* old and new chi-squared */
    real  r;                                           /* mean radius of ring */
    real  a11,a12,a22,sigma2;                              /* matrix elements */
    real  sinp, cosp, cosi, xc1, xc2, yc1, yc2;    	     /* ring elements */
    real  resmean, ressig, ratio;
    int   nblank;
    int   i, j, n;                            /* n=number of points in a ring */
#if 0
    real  x[2*MAXPTS],y[MAXPTS],w[MAXPTS];/* arrays for coords, vels and wgts */
    real  wb[MAXPTS];
    int   idx[2*MAXPTS];
    int   iblank[MAXPTS];
    real  res[MAXPTS];                        /* array for y-x, the residuals */
#else
    /* if MAXPTS too large, most OS won't allow this amount on stack/heap     */
    /* so lets allocate() it                                                  */
    static real  *x, *y, *w, *wb, *res;
    static int   *idx, *iblank;
    static int npts = 0;

    if (npts==0) {
      if (Qimage) 
	npts = Nx(velptr)*Ny(velptr);
      else
	npts = MAXPTS;
      dprintf(0,"rotfit: scratch arrays with npts=%d\n",npts);
      x = (real *) allocate(2*npts*sizeof(real));
      y = (real *) allocate(npts*sizeof(real));
      w = (real *) allocate(npts*sizeof(real));
      wb = (real *) allocate(npts*sizeof(real));
      res = (real *) allocate(npts*sizeof(real));
      idx = (int *) allocate(2*npts*sizeof(int));
      iblank = (int *) allocate(npts*sizeof(int));
    } 

#endif

#if 0
    /* Qrotcur */
    for (i=0; i<GPARAM; i++) mask[i] = 0;    
#endif
   
    nfr=0;                                 /* reset number of free parameters */
    for (i=0; i<nparams; i++) {
        nfr += mask[i];                /* calculate number of free parameters */
        eps[i] = 0.1;                           /* and init the eps parameter */
    }
    r=0.5*(ri+ro);                                     /* mean radius of ring */

    printf(" Disk radii range: %g %g\n",ri,ro); 
    /* printf(" %4d  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f", */
    printf(" iter   VSYS     XPOS     YPOS       PA      INC  ");
    /*      xdddd__xxxx.xx__xxxx.xx__xxxx.xx__xxxx.xx__xxxx.xx__ */
    /*    0       0.00    30.00    60.00    63.50    63.50    5.00    1.00    2.00    1.00     */
    for (i=GPARAM; i<nparams; i++)
      printf("    P%d  ",i+1-GPARAM);
    printf("  points   sigma\n");

    printf("  \n");
#if 1
    printf("MASK vector: ");
    for (i=0; i<nparams; i++)
      printf(" %d",mask[i]);
    printf("  (nfr=%d)\n",nfr);
#endif

    getdat(x,y,w,idx,res,&n,npts,p,ri,ro,thf,wpow,&q,side,&full,nfr,0);
    for (i=0;i<n;i++) iblank[i] = 1;  /* 1 means pixel is ok, 0 means flagged */

    h=0;                                          /* reset itegration counter */
    nblank=0;

    perform_out(h,p,n,q);                             /* show first iteration */
    do {                                   /* and start the outer REPEAT loop */
         h++;                                               /* next iteration */
         chi=q;                                           /* save chi-squared */
         for(k=0;k<nparams;k++)              /* loop to save initial estimates */
            pf[k]=p[k];                                       /* save p in pf */
         nrt = nllsqfit(x,2,y,w,res,n,pf,e,mask,nparams,tol,itmax,lab,vobs,vobsd);
	 if (nsigma > 0)  {              /* do another fit w/ outliers marked */
	   nblank=set_blank(n,res,w,iblank,nfr,nsigma,&q);
	   if (nblank)
	     nrt = nllsqfit(x,2,y,w,res,n,p,e,mask,nparams,tol,itmax,lab,
			    vobs,vobsd);
	   if (Qimage) {
	     if (h==1)
	       warning("new feature: patching the input maps for nsigma rejects");
	     for (i=0; i<n; i++) 
	       if (iblank[i] == 0)
		 MapValue(velptr,idx[2*i],idx[2*i+1]) = undf;
   	   } else {
	     warning("sadly we don't have a good way to iterate w/ nsigma on tables");
	     break;
	   }
	 }
         if (nrt<0) {
            warning("nllsqfit=%d: must find better solution (n=%d)",nrt,n);
            break;    /* ???? could also try and continue here (like KGB did) */
	 }
         for (k=0; k<nparams; k++) 
            df[k]=pf[k]-p[k];                  /* calculate difference vector */
         flip=1.0;                                   /* factor for inner loop */
         for(;;) {                                        /* inner WHILE loop */
            stop = FALSE;
            for(k=0; k<nparams; k++)       /* loop to calculate new parameters */
               pf[k]=flip*df[k]+p[k];                       /* new parameters */
            pf[4]=MIN(pf[4],180.0-pf[4]);         /* in case inclination > 90 */
            getdat(x,y,w,idx,res,&n,npts,pf,ri,ro,thf,wpow,&q,side,&full,nfr,0);
#if 0
	    for (i=0;i<n;i++) iblank[i] = 1; 
	    for (i=0;i<n;i++) w[i] *= iblank[i];            /* apply blanking */
#endif
            if (q < chi) {                                     /* better fit ?*/
               perform_out(h,pf,n,q);                   /* show the iteration */
               for(k=0;k<nparams;k++)            /* loop to save new estimates */
                  p[k]=pf[k];
                stop=FALSE;  /* but make sure it doesn't quit from outer loop */
                dprintf(1,"Found a solution q=%g chi=%g\n",q,chi);
                break;                             /* leave inner XWHILE loop */
            } else {                                      /* not a better fit */
               stop=TRUE;                                    /* reset logical */
               for (k=0; k<nparams; k++)            /* loop through parameters */
                  stop = stop && (ABS(flip*df[k]) < eps[k]);
               stop = (ABS(q-chi)/chi < tol);      /* difference small enough */
               if (stop) {                                      /* is it so ? */
		 dprintf(1,"Chi^2 small enough anyhow....%g\n",(q-chi)/chi);
		 break;                             /* we found it ! XREPEAT */
               } else
		 dprintf(1,"Chi^2 not small...%g\n",(q-chi)/chi);
            }
            if (flip > 0.0)                               /* make it negative */
               flip *= -1.0;
            else                                          /* make it positive */
               flip *= -0.5;
         } /* end of inner loop */
         if (stop) break;   /* C-kludge to leave outer loop from within inner */
    } while(h < itmax) ;                                /*  end of outer loop */

    if (nrt<0)                                              /* error from fit */
        ier = nrt;                 /* error equals error return from nllsqfit */
    else if (h==itmax)                        /* maximum number of iterations */
        ier=-10;
    else
        ier=1;                                          /* signal all is well */

    /* THESE VERBAL MESSAGES SEEM OUT OF DATE WITH NLLSQFIT MESSAGES */

    switch(ier) {           /* process some error codes, see also nllsqfit(3) */
      case -1:                                       /* no degrees of freedom */
         error("ROTCURSHAPE: no degrees of freedom");
         break;
      case -2:
         error("ROTCURSHAPE: you have no free parameters");
         break;
      case -3:                                        /* problems with matrix */
         error("ROTCURSHAPE: not enough degrees of freedom");
         break;
      case -4:
         error("ROTCURSHAPE: problems with matrix inversion");
         break;
      case -5:                                        /* floating zero divide */
         error("ROTCURSHAPE: floating zero divide");
         break;
      case -10:                     /* maximum number of iterations too small */
         error("ROTCURSHAPE: max. number of iterations %d to small",itmax);
         break;
      default:
         perform_out(h,p,n-nblank,q);                  /* write final results */
         if (full)
            warning("not all points inside ring %g were used",r);
    }

    if (lunres) {
      if (Qimage) {
	if (resptr == NULL) 
	  copy_image(velptr,&resptr);
	/* should this be done each call, or only when resptr = NULL ? */
	for (i=0;i<Nx(resptr);i++)
	  for (j=0;j<Ny(resptr);j++)
	    MapValue(resptr,i,j) = 0.0;
      } else {
        fprintf(lunres,"#  %d : New ring %g - %g\n",n,ri,ro);
        fprintf(lunres,"#  Xsky Ysky Vobs Vobs-Vmod Xgal Ygal Rgal THETAgal\n");
      }
      if (Qfit) {
	warning("One more fit iteration for full fit field map, old map=%d points",n);
	getdat(x,y,w,idx,res,&n,npts,p,ri,ro,thf,wpow,&q,side,&full,nfr,1);
	warning("Found %d points",n);
      } else
	warning("Old map Found %d points",n);
      cosp = cos((p[3]+90)*F);
      sinp = sin((p[3]+90)*F);
      cosi = cos(p[4]*F);
      for (i=0; i<n; i++) {
	xc1 = x[2*i]-dx*p[1];                /* X and Y (sky) w.r.t. rotation center */
	yc1 = x[2*i+1]-dy*p[2];
	xc2 =   xc1*cosp + yc1*sinp;         /* X and Y (galaxy) w.r.t. rotation center */
	yc2 = (-xc1*sinp + yc1*cosp)/cosi;   /* Qrotcur */
	if (Qimage) {
	  MapValue(resptr,idx[2*i],idx[2*i+1]) = res[i];
	  if (i==0)
	    MapMin(resptr) = MapMax(resptr) = res[i];
	  else {
	    MapMin(resptr) = MIN(MapMin(resptr), res[i]);
	    MapMax(resptr) = MAX(MapMax(resptr), res[i]);
	  }
	} else {
	  fprintf(lunres,"%g %g %g %g %g %g %g %g\n",
		xc1, yc1, y[i], res[i],
		xc2, yc2, sqrt(xc2*xc2+yc2*yc2),atan2(yc2,xc2)/F);
	}
      }
      if (Qimage) write_image(lunres,resptr);
    }

    if (ier==1 && cor[0]>0 && cor[1]>0) {   /* calculate ellipse parameters ? */
         sigma2=0.0;                                      /* reset parameters */
         a11=0.0;
         a12=0.0;
         a22=0.0;
         for (i=0; i<n; i++) {                    /* loop through data points */
            vobsd(&x[2*i],p,b,nparams);            /* calculate derivatives */
            a11=a11+w[i]*b[cor[0]-1]*b[cor[0]-1];   /* calc ellips parameters */
            a22=a22+w[i]*b[cor[1]-1]*b[cor[1]-1]; 
            a12=a12+w[i]*b[cor[0]-1]*b[cor[1]-1]; 
            sigma2=sigma2+w[i]*sqr(y[i]-vobs(&x[2*i],p,nparams));
         }
         sigma2=sigma2/(real)n;                        /* calculate new sigma */
         elp4[0]=a11;                              /* save ellipse parameters */
         elp4[1]=a12;
         elp4[2]=a22;
         elp4[3]=sigma2;
         printf("  ===> Ellipse error: (elp4=%g %g %g %g)\n",
                        elp4[0], elp4[1], elp4[2], elp4[3]);
    }
    *npt = n;
    return ier;
} /* rotfit */

int perform_out(int h,real *p,int n,real q)
{
  int i;
  printf(" %4d  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f",
	 h,    p[0],   p[1],  p[2],  p[3],  p[4]);
  for (i=GPARAM; i<nparams; i++)  
    printf(" %7.2f",p[i]);
  printf("  %5d %8.3f\n", n,   q);
}


int set_blank(int n, real *res, real *w, int *iblank, int nfr, real nsigma, real *q)
{
  real resmean, ressig, ratio;
  int i, nblank;

  stat2(res,n,nfr,&resmean,&ressig);
  dprintf(1,"Refit: mean=%g sigma=%g nsigma=%g\n",resmean,ressig,nsigma);
  if (ABS(resmean) > ressig) 
    warning("Residuals not gaussian? mean=%g/sigma=%g => %g",
	    resmean,ressig,ABS(resmean)/ressig);
  nblank = 0;
  /* note that this loop ensures iblank[] is completely reset for all elements */
  for (i=0;i<n;i++) {
    ratio = ABS(res[i]-resmean)/ressig;
    if (ratio > nsigma) {
      iblank[i] = 0;
      w[i] = 0.0;
      nblank++;
    } else {
      iblank[i] = 1;
    }
  }
  if (nblank) {
    dprintf(0,"Rejecting %d/%d points for nsigma=%g\n",nblank,n,nsigma);
    stat2b(res,iblank,n,nfr,&resmean,&ressig);
    dprintf(1,"Refit-blanking: mean=%g sigma=%g nsigma=%g\n",resmean,ressig,nsigma);
    *q = ressig;
  }
  return nblank;
}

/******************************************************************************/

/*
 * ROTPLT originally plots, but now only produces a table if lunpri
 *  is set. Otherwise it exits.
 */

void rotplt(rad,vsy,evs,pan,epa,inc,ein,xce,exc,yce,eyc,p,e,
       mask,ifit,elp,lunpri,cor,npt,factor)
real rad[],vsy[],evs[],pan[],epa[],inc[],ein[],p[],e[],
     xce[],exc[],yce[],eyc[],elp[][4], factor;
int  mask[],ifit,cor[],npt[];
stream lunpri;
{
  int i,j;
  real mean, sig;

  if (lunpri==NULL) return;         /* determine if work to be done */
  if (ifit<0) return;               /* else return right now */

  fprintf(lunpri,"\n");
  for (i=0; i<ifit; i++) {
    fprintf(lunpri,"VSYS: %g %g\n",vsy[i],evs[i]);
    fprintf(lunpri,"XPOS: %g %g\n",xce[i],exc[i]);
    fprintf(lunpri,"YPOS: %g %g\n",yce[i],eyc[i]);
    fprintf(lunpri,"PA:   %g %g\n",pan[i],epa[i]);
    fprintf(lunpri,"INC:  %g %g\n",inc[i],ein[i]);
    for (j=GPARAM; j<nparams;j++)
      fprintf(lunpri,"P%d:  %g %g\n",j-GPARAM+1,p[j],e[j]*factor);
    fprintf(lunpri,"NPT:  %d\n",npt[i]);
  }
  fprintf(lunpri,"\n");
  fprintf(lunpri," Error correction factor: : %g\n",factor);

  if (ifit < 2) return;
   
  stat2(inc,ifit,1,&mean,&sig);
  fprintf(lunpri," average inclination      : %8.2f  (%8.3f)  degrees\n",
                                              mean,    sig);
  stat2(pan,ifit,1,&mean,&sig);
  fprintf(lunpri," average position angle   : %8.2f  (%8.3f)  degrees\n",
                                              mean,    sig);
  stat2(vsy,ifit,1,&mean,&sig);
  fprintf(lunpri," average systemic velocity: %8.2f  (%8.3f)  km/s\n",
                                              mean,    sig);
  stat2(xce,ifit,1,&mean,&sig);
  fprintf(lunpri," average x-position       : %8.2f  (%8.3f)  grids\n",
                                              mean,    sig);
  stat2(yce,ifit,1,&mean,&sig);
  fprintf(lunpri," average y-position       : %8.2f  (%8.3f)  grids\n",
                                              mean,    sig);
}

/* stat2: for an input array a[] of length n, get the mean and dispersion */

void stat2(real *a,int n,int nf, real *mean,real *sig)
{
    real s=0, sx=0, sxx=0;
    int i;

    if (n<=0) return;
    for (i=0; i<n; i++) {
        s += 1.0;
        sx += *a;
        sxx += sqr(*a++);
    }
    *mean = sx/s;
    if (n > 1) {
      s -= nf;
      *sig = sqrt( (sxx-s*sqr(*mean)) / MAX(1.0,s) );
    } else
      *sig = 0.0;
}

void stat2b(real *a,int *b,int n,int nf, real *mean,real *sig)
{
    real s=0, sx=0, sxx=0;
    int i;

    if (n<=0) return;
    for (i=0; i<n; i++) {
      if (b[i]) {
        s += 1.0;
        sx += *a;
        sxx += sqr(*a++);
      }
    }
    *mean = sx/s;
    if (s > 0) {
      s -= nf;
      *sig = sqrt( (sxx-s*sqr(*mean)) / MAX(1.0,s) );
    } else
      *sig = 0.0;
}

/* 
 *    GETDAT gets data from disk and calculates differences.
 *
 *    SUBROUTINE GETDAT(X,Y,W,IDX,RES,N,NMAX,P,RI,RO,THF,WPOW,Q,SIDE,FULL,NFR,MODE)
 *
 *    X        real array       sky coordinates of pixels inside ring
 *    Y        real array       radial velocities
 *    W        real array       weights of radial velocities
 *    IDX      integer array    sky coordinates in integer pixels
 *    RES      real array       residual velocities
 *    N        integer          number of pixels inside ring
 *    NMAX     integer          maximum number of pixels
 *    P        real array       parameters of ring
 *    RI       real             inner radius of ring
 *    RO       real             outer radius of ring
 *    THF      real             free angle around minor axis
 *    WPOW     integer          weighting mode
 *    Q        real             chi-squared
 *    SIDE     integer          receding or approaching side
 *    FULL     logical          too many points in ring
 *    NFR      integer          number of degrees of freedom
 *    MODE     integer          fitmode (0=fitting w/ masking 1=fit, no mask)
 */

void getdat(x,y,w,idx,res,n,nmax,p,ri,ro,thf,wpow,q,side,full,nfr,mode)
int   *n,nmax;       /* number of pixels loaded (O), max. number (I) */
int   nfr;           /* number of degrees of freedom */
int   wpow;          /* weigthing function */
int   side;          /* part of galaxy */
bool  *full;         /* output if data overflow */
real  x[],y[],w[];   /* arrays to store coords., vels. and weights */
real  res[];         /* array to store residuals "as is" (no fitting) */
int   idx[];         /* coordinates in pixel units (0,0 is lower left) */
real  p[];           /* parameters of ring */
real  ri,ro;         /* inner and outer radius of ring */
real  thf;           /* free angle */
real  *q;            /* output sigma */
int   mode;          /* fit mode */
{
/******************************************************************************/
    int   nlt,nmt;                                                /* counters */
    int   mdone,mstep,mleft,m,m1,m2,l,i,j;                        /* counters */
    bool  use;                                    /* boolean (for data point) */
    real  phi,inc,x0,y0,sinp,cosp,sini,cosi;        /* parameters for ellipse */
    real  a,b,s;                                 /* couple of dummy variables */
    real  xx[2],dn[2];       /* arrays store coordinates and dN/dx/N, dN/dy/N */
    real  free;                     /* relative free angle (for simple check) */
    real  v;                                                          /* vel. */
    real  wi;                                         /* weight of data point */
    real  theta,costh,xr,yr,r,rx,ry;        /* coordinates in plane of galaxy */
    int   llo,lhi,mlo,mhi;

    *full = FALSE;        /* reset logical for output */
    *q=0.0;               /* reset sigma for output */
    *n=0;                 /* reset number of pixels for output */

    phi=p[3]+pamp;        /* position angle plus map p.a. */
    inc=p[4];             /* inclination */
    x0=p[1];              /* x-position of center */
    y0=p[2];              /* y-position of center */
    free=ABS(sin(F*thf)); /* free angle in terms of sine */
    sinp=sin(F*phi);       /* sine of pa. */
    cosp=cos(F*phi);       /* cosine of pa. */
    if (Qrachel) {
      sini=cosi=1.0;
    } else {
      sini=sin(F*inc);       /* sine of inc. */
      cosi=cos(F*inc);       /* cosine of inc. */
    }
    a=sqrt(1.0-cosp*cosp*sini*sini);       /* find square around ellipse */
    b=sqrt(1.0-sinp*sinp*sini*sini);

    if (Qimage) {                         /* Image input data */
      /* BUG:  need to fix this, there is a rounding problem near center */
      /*       although this is where rotcur isn't quite all that valid  */
      if (mode) {
	llo = lmin;                             /* the whole map */
	mlo = mmin;
	lhi = lmax;
	mhi = mmax;
      } else {
	llo=MAX(lmin,(int)(x0-a*ro/dx));        /* square around ellipse */
	mlo=MAX(mmin,(int)(y0-b*ro/dy));
	lhi=MIN(lmax,(int)(x0+a*ro/dx));
	mhi=MIN(mmax,(int)(y0+b*ro/dy));
      }
      if (llo > lhi || mlo > mhi) {
	error("ring not inside map [%d %d %d %d]",llo,mlo,lhi,mhi);
	*n = 0;                        /* continue here ??? */
        return;
      } else {
        dprintf(1,"getdat: GEOM-Par=%g %g %g %g %g RC-Par=%g\n",
		p[0],p[1],p[2],p[3],p[4],p[5]);
        dprintf(1,"getdat: box [%d %d %d %d]\n",llo,mlo,lhi,mhi);
        dprintf(1,"getdat: box [%g %g %g %g]\n",
		x0-a*ro/dx,y0-b*ro/dy,x0+a*ro/dx,y0+b*ro/dy);
      }
      
      nlt=lmax-lmin+1;       /* number of pixels in X (not used) */
      nmt= mhi- mlo+1;       /* number of pixels in Y */
      
      mdone=0;      /* number of lines done */
      mstep=1;      /* number of lines to read in one call (old GIPSY relic) */
      mleft=nmt;    /* number of lines left to do */
      do {                      /* big loop */
	mstep=MIN(mstep,mleft);     /* number of lines to do in this run */
	m1=mlo+mdone;       /* first line in this run */
	m2=m1+mstep-1;      /* last line in this run */
	m=m1-1;             /* line counter */
	while (m < m2) {       /* loop */
	  m=m+1;                 /* increment line counter */
	  ry=dy*(real)(m);       /* Y position in plane of galaxy */
	  l=llo-1;               /* line position counter */
	  while (l < lhi) {   /* loop */
	    l++;        /* increment line position counter */
	    rx=dx*(real)(l);       /* X position in plane of galaxy */
	    /***** i=(m-m1)*nlt+l-lmin+1;        --> array pointer */
	    v = MapValue(velptr,l,m);        /* velocity at (l,m) */
	    if (mode || v != undf) {       /* undefined value ? */
	      xr=(-(rx-dx*x0)*sinp+(ry-dy*y0)*cosp);     /* X position in galplane */
	      yr=(-(rx-dx*x0)*cosp-(ry-dy*y0)*sinp)/cosi;
	      r=sqrt(xr*xr+yr*yr);                       /* distance from center */
	      if (r < 0.01)                              /* radius too small ? */
		theta=0.0;
	      else
		theta=atan2(yr,xr)/F;
	      costh=ABS(cos(F*theta));       /* calculate |cos(theta)| */
	      dprintf(5,"@ %d,%d : r=%g cost=%g xr=%g yr=%g\n",l,m,r,costh,xr,yr);
	      if (mode || (r>ri && r<ro && costh>free)) {     /* point inside ring ? */
		dprintf(5," ** adding this point\n");
		if (denptr) 
		  wi = MapValue(denptr,l,m);
		else
		  wi=1.0;                /* calculate weight of this point */
		
		for (i=0; i<wpow; i++) 
		  wi *= costh;
		xx[0]=rx;       /* x position */
		xx[1]=ry;       /* y position */
		v -= bmcorr(xx,p,l,m);  /* beam-correction factor */
		use=FALSE;        /* reset logical */
		switch(side) {    /* which side of galaxy */
		case 1:             /* receding half */
		  use=(ABS(theta)<=90.0);        /* use this data point ? */
		  break;
		case 2:             /* approaching half */
		  use=(ABS(theta)>=90.0);        /* use this data point ? */
		  break;
		case 3:             /* both halves */
		  use=TRUE;         /* use this data point */
		  break;
		default:
		  error("wrong side (%d) of galaxy",side);
		}
		if (use) {
		  *full = (*n==(nmax-1));    /* buffers full */
		  if (! *full) {         /* save data ? */
		    x[*n*2]=rx;        /* load X-coordinate */
		    x[*n*2+1]=ry;      /* load Y-coordinate */
		    y[*n]=v;           /* load radial velocity */
		    w[*n]=wi;          /* load weight */
		    idx[*n*2]=l;       /* load integer X-coordinate */
		    idx[*n*2+1]=m;     /* load integer Y-coordinate */
		    if (mode)
		      s=vobs(xx,p,nparams);      /* fitted map */
		    else 
		      s=(v-vobs(xx,p,nparams));  /* corrected difference */
		    res[*n] = s;               
		    *q += s*s*wi;       /* calculate chi-squared */
		    *n += 1;           /* increment number of pixels */
		  }
		}
	      }
	      if (*full) break;      /* if buffers are filled - quit */
	    } /* v != undf */
	    if (*full) break;       /* if buffers are filled - quit */
	  }  /* l-loop */
	  if (*full) break;          /* if buffers are filled - quit */
	} /*  mloop */
	mdone += mstep;      /* number of lines done */
	mleft -= mstep;      /* number of lines left */
      } while (mleft>0);        /* big loop */
    } else {                                               /* tabular input data */
      for (i=0; i<n_vel; i++) {
	rx = xpos_vel[i];
	ry = ypos_vel[i];
	v  = vrad_vel[i];
	if (v == undf) continue;
	xr=(-(rx-x0)*sinp+(ry-y0)*cosp);     /* X position in galplane */
	yr=(-(rx-x0)*cosp-(ry-y0)*sinp)/cosi;
	r=sqrt(xr*xr+yr*yr);                       /* distance from center */
	if (r < 0.01)                              /* radius too small ? */
	  theta=0.0;
	else
	  theta=atan2(yr,xr)/F;
	costh=ABS(cos(F*theta));       /* calculate |cos(theta)| */
	if (r>ri && r<ro && costh>free) {     /* point inside ring ? */
	  dprintf(5,"@ (%g,%g,%g) r=%g cost=%g xr=%g yr=%g\n",rx,ry,v,r,costh,xr,yr);
	  dprintf(5," ** adding this point\n");

	  if (verr_vel)
	    wi = 1.0/sqr(verr_vel[i]);
	  else
	    wi=1.0;                /* calculate weight of this point */
	  for (j=0; j<wpow; j++) 
	    wi *= costh;
	  xx[0]=rx;       /* x position */
	  xx[1]=ry;       /* y position */
	  /* no beam correction here */
	  use=FALSE;        /* reset logical */
	  switch(side) {    /* which side of galaxy */
	  case 1:             /* receding half */
	    use=(ABS(theta)<=90.0);        /* use this data point ? */
	    break;
	  case 2:             /* approaching half */
	    use=(ABS(theta)>=90.0);        /* use this data point ? */
	    break;
	  case 3:             /* both halves */
	    use=TRUE;         /* use this data point */
	    break;
	  default:
	    error("wrong side (%d) of galaxy",side);
	  }
	  if (use) {
	    *full = (*n==(nmax-1));    /* buffers full */
	    if (! *full) {         /* save data ? */
	      x[*n*2]=rx;        /* load X-coordinate */
	      x[*n*2+1]=ry;      /* load Y-coordinate */
	      y[*n]=v;           /* load radial velocity */
	      w[*n]=wi;          /* load weight */
	      s=(v-vobs(xx,p,nparams));  /* corrected difference */
	      res[*n] = s;               
	      *q += s*s*wi;       /* calculate chi-squared */
	      *n += 1;           /* increment number of pixels */
	    }
	  }
	  if (*full) break;      /* if buffers are filled - quit */
	} /* inside ring ? */
      } /* loop over all points */
    }

    if (*n > nfr)  /* enough data points ? */
      *q=sqrt(*q/(real)(*n-nfr));     /* calculate sigma */
    else
      *q=1.0e+30;      /* some extreme value */
} /* getdat */

real bmcorr(xx,p,l,m)
real xx[2];     /* real x and y in galaxy plane */
real p[];       /* velocity field parameters */
int l,m;        /* grid coordinates w.r.t. 0,0 */
{
    real d, vdif, dn[2];

    if (Qbeam && denptr) {       /* beam-correction wanted ? */
        d = MapValue(denptr,l,m);
        dn[0] = (MapValue(denptr,l+1,m)-MapValue(denptr,l-1,m))/d/grid[0];
        dn[1] = (MapValue(denptr,l,m+1)-MapValue(denptr,l,m-1))/d/grid[1];
        vcor(xx,p,&vdif,dn); /*  calculate correction */
        return vdif;
    } else
        return 0.0;
}

/*  External functions needed by the nllsqfit routine */


/* global parameters for the local nllsqfit routines */

int  p_init = 1;
int  i,j;                                                    /* loop counters */
real vs,vc,phi,inc;                           /* parameters of velocity field */
real cosp1,cosp2,sinp1,sinp2,cosi1,cosi2,sini1,sini2;     /* different angles */
real x,y;                                                  /* sky coordinates */
real cost1,cost2,sint1,sint2,xx1,yy1,r,r1;  /* coordinates in plane of galaxy */
real fc,t[5],q[5],bx1,bx2,by1,by2;    /* vars for calculating beam-correction */
real vn,v2;                                  /* correction to radial velocity */
real rcderv[MAXPAR*MAXMOD];
int nrcpar = 0;

/*
 *
 *    VOBS calculates radial velocity from rotation curve.
 *
 *      real vobs(C,P,M)
 *
 *    C       real array       grid position in plane of galaxy
 *    P       real array       parameters of ring
 *    M       integer          dummy
 * ---------------------------------------------------------------------------
 *    VOBSD calculates the  partial  derivatives with respect
 *    to the parameters.
 *
 *      void vobsd(C,P,D,M)
 *
 *    C        real array       grid position in plane of galaxy
 *    P        real array       parameters of ring
 *    D        real array       for partial derivatives
 *    M        integer          dummy variable
 * ---------------------------------------------------------------------------
 *    VCOR calculates the beam smearing correction.
 *
 *      void vcor(C,P,VD,DN)
 *
 *    C        real array       grid position in plane of galaxy
 *    P        real array       parameters of ring
 *    VD       real             correction term
 *    DN       real array       contains dN/dx/N and dN/dy/N
 */

real vobs_c1(real *c,real *p,int m)
{
      perform_init(p,c);
      return vs+vc*sini1*cost1;          /* circular motion */
}


void vobsd_c1(real *c,real *p,real *d,int m)
{
  int i;
  perform_init(p,c);

    /* disk geometry parameters (Vsys,X0,Y0,PA,INC) */

  d[0]=1.0;                    /* partial derivative with respect to Vsys */
  d[1]= grid[0]*vc*(sint1*sinp1-cost1*cosp1/cosi1)*sint1*sini1*r1;  /* X0 */
  d[2]=-grid[1]*vc*(sint1*cosp1+cost1*sinp1/cosi1)*sint1*sini1*r1;  /* Y0 */
  d[3]=F*vc*(cosi2+sini2*cost2)*sint1*sini1/cosi1;                  /* PA */
  d[4]=F*vc*(cosi2-sini2*sint2)*cost1/cosi1;                       /* Inc */
  
  /* d[5] and up are rotation curve parameters */

  for (i=GPARAM; i<m; i++)
    d[i] = rcderv[i-GPARAM] *sini1*cost1;
}

void vcor_c1(real *c,real *p,real *vd,real *dn)
{
    perform_init(p,c);

    bx1=G*beam[0]*r1;       /* modified beam */
    bx2=bx1*bx1;
    by1=G*beam[1]*r1;       /* modified beam */
    by2=by1*by1;

    /* perform the vncorr */

    vn=0.0;        /* reset variables */

    t[0]=sint2;                         /*  fill T  */
    t[1]=-sint1*cost1;

    /* warning("Why does this loop i<1 ?? "); */
    for (i=0; i<1; i++){        /* loop for calculating correction  */
        if (i==0) {
            fc=dn[i]*vc*sini1*bx2*r;     /*  factor  */
            q[0]=-sinp1;                   /*   fill Q  */
            q[1]=-cosp1/cosi1;
        } else {
            fc=dn[i]*vc*sini1*by2*r;     /*  factor  */
            q[0] = cosp1;                   /*  fill Q  */
            q[1] =-sinp1/cosi1;
        }
        for (j=0; j<2; j++)         /* loop for calculating correction  */
            vn += fc*q[j]*t[j];      /* correction term  */
    }

    /* perform the v2corr   */

    v2=0.0;                 /* reset variables */

    t[0]=-3.0*cost1*sint2;                /* fill T */
    t[1]=sint1*(2.0*cost2-sint2);
    t[2]=cost1*(2.0*sint2-cost2);
    for (i=0; i<2; i++) {      
        if (i==0) {
            fc=0.50*vc*sini1*bx2;
            q[0]=sinp2;
            q[1]= 2.0*sinp1*cosp1/cosi1;
            q[2]=cosp2/cosi2;
        } else {
            fc=0.50*vc*sini1*by2;
            q[0]=cosp2;
            q[1]=-2.0*sinp1*cosp1/cosi1;
            q[2]=sinp2/cosi2;
        }
        for (j=0; j<3; j++)  /* loop to calculate correction */
            v2 += fc*q[j]*t[j];
    }

    *vd=vn+v2;      /* corrected velocity        */
}

/*
 *  The following functions are used to fit the sin(theta) 
 *  part - useful for pure radial out(in) flow instead or 
 *  rotation
 */

real vobs_s1(real *c,real *p,int m)
{
      error("cannot do inflow yet");
      perform_init(p,c);
      return vs+vc*sini1*sint1;      /* radial outflow motion */
}

void vobsd_s1(real *c,real *p,real *d,int m)
{
    perform_init(p,c);
    error("cannot do inflow yet");
    d[0]=1.0;                    /* partial derivative with respect to Vsys */
    d[1]= grid[0]*vc*(0);
    d[2]=-grid[1]*vc*(0);
    d[4]=F*vc*(sint2*sini2-1)*cost1*sini1/cosi1;                      /* PA */
    d[4]=F*vc*(1/cosi2-sini2*sint2)*sint1/cosi1;                     /* INC */
    d[5]=sini1*sint1;                                               /* Vrot */
}


void vcor_s1(real *c,real *p,real *vd,real *dn)
{
    error("vcor_s1 NOT IMPLEMENTED");

    perform_init(p,c);

    bx1=G*beam[0]*r1;       /* modified beam */
    bx2=bx1*bx1;
    by1=G*beam[1]*r1;       /* modified beam */
    by2=by1*by1;

    /* perform the vncorr */

    vn=0.0;        /* reset variables */

    t[0]=sint2;                         /*  fill T  */
    t[1]=-sint1*cost1;

    for (i=0; i<1; i++){        /* loop for calculating correction  */
        if (i==0) {
            fc=dn[i]*vc*sini1*bx2*r;     /*  factor  */
            q[0]=-sinp1;                   /*   fill Q  */
            q[1]=-cosp1/cosi1;
        } else {
            fc=dn[i]*vc*sini1*by2*r;     /*  factor  */
            q[0] = cosp1;                   /*  fill Q  */
            q[1] =-sinp1/cosi1;
        }
        for (j=0; j<2; j++)         /* loop for calculating correction  */
            vn += fc*q[j]*t[j];      /* correction term  */
    }

    /* perform the v2corr   */

    v2=0.0;                 /* reset variables */

    t[0]=-3.0*cost1*sint2;                /* fill T */
    t[1]=sint1*(2.0*cost2-sint2);
    t[2]=cost1*(2.0*sint2-cost2);
    for (i=0; i<2; i++) {      
        if (i==0) {
            fc=0.50*vc*sini1*bx2;
            q[0]=sinp2;
            q[1]= 2.0*sinp1*cosp1/cosi1;
            q[2]=cosp2/cosi2;
        } else {
            fc=0.50*vc*sini1*by2;
            q[0]=cosp2;
            q[1]=-2.0*sinp1*cosp1/cosi1;
            q[2]=sinp2/cosi2;
        }
        for (j=0; j<3; j++)  /* loop to calculate correction */
            v2 += fc*q[j]*t[j];
    }

    *vd=vn+v2;      /* corrected velocity        */
}


int perform_init(real *p,real *c)
{
  int i;

  vs=p[0];                 /* systemic velocity */
  if (p[3] != phi || p_init) {   /* new position angle ? */
    phi=p[3];               /* position angle */
    cosp1=cos(F*phi);       /* cosine */
    cosp2=cosp1*cosp1;      /* cosine squared */
    sinp1=sin(F*phi);       /* sine */
    sinp2=sinp1*sinp1;      /* sine squared */
  }
  if (p[4] != inc || p_init) {   /* new inclination ?  */
    inc=p[4];               /* inclination */
    cosi1=cos(F*inc);       /* cosine */
    cosi2=cosi1*cosi1;      /* cosine squared */
    if (Qrachel)
      sini1=sini2=1.0;
    else {
      sini1=sin(F*inc);       /* sine */
      sini2=sini1*sini1;      /* sine squared */
    }
  }
  p_init = 0;                   /* academic case :-) */
  x=c[0]-p[1]*grid[0];          /* calculate x */
  y=c[1]-p[2]*grid[1];          /* calcualte y */
  xx1=(-x*sinp1+y*cosp1);        /* x in plane of galaxy */
  yy1=(-x*cosp1-y*sinp1)/cosi1;  /* y in plane of galaxy */
  r=sqrt(xx1*xx1+yy1*yy1);          /* distance from center */
  r1=1/r;                         /* save inverse for speed */
  cost1=xx1*r1;         /* cosine of angle in plane of galaxy */
  sint1=yy1*r1;         /* sine of angle in plane of galaxy */
  cost2=cost1*cost1;     /* cosine squared */
  sint2=sint1*sint1;     /* sine squared */

  if (nmod > 1) {
    vc = 0;
    for (i=0; i<nmod; i++)
      vc += sqr( (*rcfn[0])(r,npar[i],&p[GPARAM+ipar[i]],&rcderv[ipar[i]]) );
    vc = sqrt(vc);
  } else
    vc = (*rcfn[0])(r,npar[0],&p[GPARAM],&rcderv[0]);
}


