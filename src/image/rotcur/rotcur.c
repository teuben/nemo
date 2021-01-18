/*******************************************************************************
 *    ROTCUR derives the kinematical parameters from the observed
 *    velocity field by fitting tilted-rings to the velocity field.  It
 *    does a least-squares-fitting to the function:
 *
 *                v(x,y) = VSYS + VROT * cos(theta) * sin(INC)
 *
 *                           - (x-XPOS) * sin(PA) + (y-YPOS) * cos(PA)
 *     with:    cos(theta) = -----------------------------------------
 *                                            r
 *
 *                           - (x-XPOS) * cos(PA) - (y-YPOS) * sin(PA)
 *     and:     sin(theta) = -------------------------------------
 *                                        r * cos(INC)
 *
 *    Here  v(x,y)  denotes  the  radial  velocity  at  rectangular sky
 *    coordinates  x and y ,  VSYS  the  systemic  velocity,  VROT  the
 *    rotational velocity,  INC  the  inclination angle  and  theta the
 *    azimuthal  distance from  the major  axis  in the  plane  of  the
 *    galaxy.  Theta is a function  of the  inclination  (INC)  and the
 *    position angle (PA) of the major axis.  XPOS and  YPOS denote the
 *    position of  the rotation centre in grids.  This program will fit
 *    for each ring the parameters  VSYS, VROT, INC, PA, XPOS and YPOS.
 *    The position  angle PA of the major axis is defined as the angle,
 *    taken in anti-clockwise  direction between the north direction on
 *    the sky and the major axis of the receding half of the galaxy.
 *
 *
 *    Author:  K. Begeman (S77) P.J. Teuben (C)
 *
 *     History: 19/jul/83 : original program                       (KGB)
 *               9/mar/85 : revision of program                    (KGB)
 *              23/may/86 : migrated to VAX-VMS                    (KGB)
 *              27/nov/88 : UNIX version                            pjt
 *               8-feb-91 : flushed buffers ROTCUR.DAT each write   pjt
 *              30-apr-91 : moved to NEMO in C                      pjt
 *               8-may-91 : working on beam smearing corrections    pjt
 *              16-may-91 : using match(*) for minimal match        pjt
 *                              <<<<bug in match>>>>
 *              12-jul-91 : fixed bug when match() not used         pjt
 *		17-oct-91 : added Npts/ring in output table (tab=)  pjt
 *		30-oct-91 : added blank= for Indians		    pjt
 *		21-may-92 : extra blank before " for formatting     pjt
 *		12-jun-92 : use inherit=t|f to take initial cond.  
 *	                    from previous ring instead              pjt
 *		17-jun-92 : try fitmode with different SIN/COS      pjt
 *				(but not implemented yet)
 *		23-jul-92 : 2.3a allow radii array to be decreasing pjt
 *		12-aug-92 : 2.3b optimized 1/r as r1		    pjt
 *                          sin(theta) fitting option
 *		20-aug-92 : 2.3e some code cleaning
 *		22-nov-94 : 2.3f fixed -DSINGLEPREC declarations    pjt
 *              15-oct-99 : 2.4  added residual computations        
 *                               to nllsqfit using resid=           pjt
 *              14-mar-01 : 2.5  added nsigma=                      pjt 
 *              23-apr-01 : 2.6  added central dv/dr estimator 
 *                               -- WORK NOT COMPLETED --     
 *                               see C-script 'rotcurcen' instead
 *               9-may-01 : 2.6a fixed error correction factor      pjt
 *               5-jun-01 : 2.7  allow density map also used as weight     PJT
 *               8-aug-01 :    a allow error correction factor 1.0  pjt
 *              26-jan-02      b allow scale factor for velocity    pjt
 *              26-jun-02 : 2.8  allow input to be an ascii table instead of image
 *              19-jul-02      a optional compilation for numrec
 *              11-sep-02   2.9  compute residuals as in rotcurshape       pjt
 *              13-nov-02        reuse= allow used points to be flagged as used (unfinished) PJT
 *              30-jan-03   2.10 allow vel.error  in tabular input if dens=t    pjt
 *              12-feb-03   2.10a  fix residual velfie as done in rotcurshape   PJT
 *
 *               4-oct-03   2.11 added option to use the WWB73 method     PJT
 *              25-may-04       a    fixed sqrt(N) problem in sigma estimate for nsigma    PJT
 *               2-jun-04   2.12 finally implemented the reuse= option     PJT
 *               5-jun-20   2.13 add wtmap= keyword                        PJT
 *              19-jan-21   2.14 revert back to old beam factor error calculation     PJT
 *
 *
 ******************************************************************************
 *
 * TODO:
 *	- keep track of pixel usage. debug output a map how many times
 *	  a pixel has been used by different rings
 *        In V2.12 this became reuse= and in essense pixel values in a used
 *        ring were set to 'undefined' so they would not be reused
 *
 *      - weights for tabular input (DONE 30-jan)
 *
 *      - clarify that our X,Y input coordinate system (even for images)
 *        is a  right-handed mathematical system, "X to the right, Y up"
 *        and that we screw the astronomical one, and force dx and dy > 0
 *        for images (mostly they will have dx < 0 , dy > 0)
 *
 *      - contrain the geometry (VSYS,XPOS,YPOS,PA,INC) to be the same
 *        over a set of rings, but still find their errors.
 *        This is a big change, and it will need new fitting functions
 *
 *      - allow expansion term, and elliptical like streaming
 *
 *      - allow a single radius (rmax) to fit a linearly rising rotation curve
 *        to match output for rotation curve plotting.
 *        Could also use rotcurves for this.
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <image.h>

/*     Set this appropriate if you want to use NumRec's mrqmin() based engine */
#if 0
#define nllsqfit nr_nllsqfit
#endif

#define PARAMS  6           /* number of parameters */

#define RING        500     /* maximum number of rings (17 arrays) */
#define MAXPTS    10000     /* maximum number of pixels per ring (4 arrays) */

#define DEF_TOL   0.001     /* tolerance for fit */
#define DEF_LAB   0.001     /* mixing parameter */
#define DEF_T     50        /* maximum number of iterations for fit */

#define F 0.0174532925  /* deg to rad conversion */
#define G 0.4246609001  /* FWHM to sigma conversion */

string defv[] = {
    "in=???\n        Input image velocity field",
    "radii=\n        Radii of rings (arcsec)",
    "vrot=\n         Rotation velocity",
    "pa=\n           Position angle (degrees)",
    "inc=\n          Inclination (degrees)",
    "vsys=\n         Systemic velocity",
    "center=\n       Rotation center (grids w.r.t. 0,0) [center of map]",
    "frang=20\n      Free angle around minor axis (degrees)",
    "side=\n         Side to fit: receding, approaching or [both]",
    "weight=\n       Weighting function: {uniform,[cosine],cos-squared}",
    "fixed=\n        Parameters to be kept fixed {vsys,vrot,pa,inc,xpos,ypos}",
    "ellips=\n       Parameters for which to plot error ellips",
    "beam=\n         Beam (arcsec) for beam correction [no correction]",
    "dens=\n         Image containing containing density map",
    "wtmap=\n        Weight map (should be mom0*mom0/mom2)",
    "tab=\n          If specified, this output table is used in append mode",
    "resid=\n        Output of residuals in a complicated plot",
    "tol=0.001\n     Tolerance for convergence of nllsqfit",
    "lab=0.001\n     Mixing parameter for nllsqfit",
    "itmax=50\n      Maximum number of allowed nllsqfit iterations",
    "units=deg,1\n   Units of input {deg, arcmin, arcsec, rad, #},{#} for length and velocity",
    "blank=0.0\n     Value of the blank pixel to be ignored",
    "inherit=t\n     Inherit initial conditions from previous ring",
    "reuse=t\n       Reuse points from previous rings if used before?",
    "fitmode=cos,1\n Basic Fitmode: cos(n*theta) or sin(n*theta)",
    "nsigma=-1\n     Iterate once by rejecting points more than nsigma resid",
    "imagemode=t\n   Input image mode? (false means ascii table)",
    "wwb73=f\n       Use simpler WWB73 linear method of fitting",
    "VERSION=2.14a\n 18-jan-2021 PJT",
    NULL,
};

string usage="nonlinear fit of kinematical parameters to a velocity field";

string cvsid="$Id$";


imageptr denptr, velptr, resptr, wtmapptr;    /* pointers to Images, if applicable */
image_maskptr   maskptr;
int   lmin, mmin, lmax, mmax;              /* boundaries of map */
real  grid[2];    /* grid separations in x and y (arcsec.) */
real  beam[2];    /* size of beam in arcseconds            */
real  dx,dy;      /* grid separation in x and y */
real  undf;       /* undefined value in map */
real  pamp;       /* position angle of map */

real  *xpos_vel, *ypos_vel, *vrad_vel, *verr_vel;
int  n_vel = 0;

real tol,lab;     /* parameters that go into nllsqfit */
int  itmax;
bool Qimage;      /* input mode */   
bool Qfirstring;
bool Qreuse;      /* reuse points from other rings ? */
bool Qwwb73;      /* use WWB73 method of linear fitting? */

/* Compute engine, derivative and beam smearing correction functions:
 * _c1 = cos(theta)	
 * _s1 = sin(theta)
 * _c2 = cos(2*theta)
 * _s2 = sin(2*theta)  etc.
 */

real vobs_c1(real *c, real *p, int m);
void vobsd_c1(real *c, real *p, real *d, int m);
void vcor_c1(real *c, real *p, real *vd, real *dn);
real vobs_s1(real *c, real *p, int m);
void vobsd_s1(real *c, real *p, real *d, int m);
void vcor_s1(real *c, real *p, real *vd, real *dn);




int rotinp(real *rad, real pan[], real inc[], real vro[], int *nring, int ring, real *vsys, 
	   real *x0, real *y0, real *thf, int *wpow, int mask[], int *side, int cor[], 
	   int *inh, int *fitmode, real *nsigma, stream lunpri);
int rotfit(real ri, real ro, real p[], real e[], int mask[], int wpow, int side, real thf, 
	   real elp4[], int cor[], int *npt, real *rms, int fitmode, real nsigma, bool useflag, stream lunres);
int perform_out(int h, real p[6], int n, real q);
void rotplt(real rad[], real vsy[], real evs[], real vro[], real evr[], real pan[], real epa[], 
	   real inc[], real ein[], real xce[], real exc[], real yce[], real eyc[], 
	   int mask[], int ifit, real elp[][4], stream lunpri, int cor[], real res[], int npt[], real factor);
void stat2(real a[], int n, real *mean, real *sig);
void getdat(real x[], real y[], real w[], int idx[], real res[], int *n, int nmax, real p[], 
	   real ri, real ro, real thf, int wpow, real *q, int side, bool *full, int nfr, bool useflag);
real bmcorr(real xx[2], real p[], int l, int m);
int perform_init(real *p, real *c);

typedef real (*my_proc1)(real *, real *, int);
typedef void (*my_proc2)(real *, real *, real *, int);
typedef void (*my_proc3)(real *, real *, real *, real *);
 
extern int nllsqfit(real *xdat, int xdim, real *ydat, real *wdat, real *ddat, int ndat, real *fpar, real *epar, int *mpar, int npar, real tol, int its, real lab, my_proc1 f, my_proc2 df);


my_proc1 vobs;
my_proc2 vobsd;
my_proc3 vcor;			/* pointers to the correct functions */


/******************************************************************************/
nemo_main()
{
    int  ier;            /* error return code */
    int  ifit=0;         /* counter for number of succesful fits */
    int  n;		 /* number of points in a ring */
    int  irng;   /* loop-counter */
    stream lunpri;       /* file for table output */
    stream lunres;       /* file for residual output */
    int  mask[PARAMS];/* mask to define the free(1) or fixed(0) parameters */
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
    real res[RING];           /* array containing rms vel in ring           */
    int  npt[RING];	      /* array containing number of points in ring  */
    real p[PARAMS],e[PARAMS]; /* arrays containing resp. pars and the errors */
    real ri,ro,r;     /* vars denoting inner, outer and mean radius */
    int  inherit;
    int  fitmode;
    real x0,y0,vsys;  /* vars for init. estim. of xpos, ypos and vsys */
    real thf;     /* var  denoting free angle around minor axis */
    real nsigma;
    real old_factor, factor;    /* factor > 1, by which errors need be multiplied */
    real rms;                  /* rms vel in a ring */

    if (hasvalue("tab"))
      lunpri=stropen(getparam("tab"),"a");  /* pointer to table stream output */
    else
      lunpri=NULL;                    /* no table output */
    Qimage = getbparam("imagemode");

    if (hasvalue("resid"))
      if (Qimage)
        lunres=stropen(getparam("resid"),"w");  /* pointer to table stream output */
      else
        lunres=stropen(getparam("resid"),"a");  /* pointer to table stream output */
    else
        lunres=NULL;                    /* no residual table output */

    rotinp(rad,pan,inc,vro,&nring,RING,     /* get input parameters */
           &vsys,&x0,&y0,&thf,
           &wpow,mask,&side,cor,&inherit,&fitmode,&nsigma,lunpri);

    old_factor = sqrt((beam[0]+grid[0])*(beam[1]+grid[1])/(grid[0]*grid[1]));
    if (beam[0] > 0 && beam[1] > 0)
      // factor = sqrt(FOUR_PI*beam[0]*beam[1]/(grid[0]*grid[1]));  /* Sicking 1997 !!! */
      factor = sqrt(1.13309*beam[0]*beam[1]/(grid[0]*grid[1]));  /* Sicking 1997 !!! */      
      //  should this not be pi/(4ln2) = 1.13309 ???
    else
      factor = 1.0;
    dprintf(0,"Sicking (1997)'s error multiplication factor=%g  (old_factor=%g)\n",
	    factor,old_factor);
    Qfirstring = TRUE;           /* for rotfit residual calc */
    for (irng=0; irng<nring-1; irng++) {  /* loop for each ring */
         ri=rad[irng];          /* inner radius of ring */
         ro=rad[irng+1];        /* outer radius of ring */
         if (ri > ro) {         /* check if need to be swapped */
            r = ri;  ri = ro;  ro = r;
         }
         r=0.5*(ri+ro);         /* mean radius of ring */

         p[0] = (mask[0] && ifit>0 && inherit) ? vsy[ifit-1] : vsys;
         p[1] = (mask[1] && ifit>0 && inherit) ? vro[ifit-1] : vro[irng];
         p[2] = (mask[2] && ifit>0 && inherit) ? pan[ifit-1] : pan[irng];
         p[3] = (mask[3] && ifit>0 && inherit) ? inc[ifit-1] : inc[irng];
         p[4] = (mask[4] && ifit>0 && inherit) ? xce[ifit-1] : x0;
         p[5] = (mask[5] && ifit>0 && inherit) ? yce[ifit-1] : y0;

         ier = rotfit(ri,ro,p,e,mask,wpow,side,thf,elp4,cor,&n,&rms,fitmode,-1.0,FALSE,lunres);
	 if (ier > 0 && nsigma > 0)
	   ier = rotfit(ri,ro,p,e,mask,wpow,side,thf,elp4,cor,&n,&rms,fitmode,nsigma,FALSE,lunres);
         if (ier>0) {           /* only if fit OK, store fit */
	   if (!Qreuse)
	     (void)rotfit(ri,ro,p,e,mask,wpow,side,thf,elp4,cor,&n,&rms,fitmode,nsigma,TRUE,lunres);
	   rad[ifit]=r;            /*  radius of ring */
	   vsy[ifit]=p[0];         /*  systemic velocity */
	   evs[ifit]=e[0]*factor;  /*  error in systemic velocity */
	   vro[ifit]=p[1];         /*  circular velocity */
	   evr[ifit]=e[1]*factor;  /*  error in circular velocity */
	   pan[ifit]=p[2];         /*  position angle */
	   epa[ifit]=e[2]*factor;  /*  error in position angle */
	   inc[ifit]=p[3];         /*  inclination */
	   ein[ifit]=e[3]*factor;  /*  error in inclination */
	   xce[ifit]=p[4];         /*  x-position */
	   exc[ifit]=e[4]*factor;  /*  error in x-position */
	   yce[ifit]=p[5];         /*  y-position */
	   eyc[ifit]=e[5]*factor;  /*  error in y-position */
	   elp[ifit][0]=elp4[0];   /*  save ellipse parameters */
	   elp[ifit][1]=elp4[1];   /* NOT corrected by 'factor' */
	   elp[ifit][2]=elp4[2];
	   elp[ifit][3]=elp4[3];
	   res[ifit] = rms;
	   npt[ifit] = n;
	   ifit++;
         }
    } /* end of loop through rings */
    if (lunres && Qimage) write_image(lunres,resptr);

    rotplt(rad,vsy,evs,vro,evr,pan,epa,         /* output the results */
           inc,ein,xce,exc,yce,eyc,
           mask,ifit,elp,lunpri,cor,res,npt,factor);
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
 *    INH      integer          inherit initial cond from previous ring?
 *    FITMODE  integer          fitmode: >0 cos-factor; <0 sin-factor
 *    LUNPRI   integer          LUN for print output
 */

rotinp(rad,pan,inc,vro,nring,ring,vsys,x0,y0,thf,wpow,mask,side,cor,inh,fitmode,nsigma,lunpri)
int  *wpow;           /* weighting function to be applied */
int  ring, *nring;    /* max. number of rings and wanted number of rings */
int  mask[];          /* mask to define free and fixed parameters */
int  *side;           /* variable denotes which part of galaxy to fit */
int  cor[];           /* plot error ellipses ? */
int  *inh;            /* inherit initial conditions from old ring? */
int  *fitmode;        /* fitmode: cos/sin and multiplicity */
real *rad;            /* user defined radii of rings */
real *vsys,vro[],pan[],inc[],*x0,*y0;  /* initial estimates (pars. 1 <--> PARAMS ) */
real *thf;            /* user defined free angle around minor axis */
real *nsigma;
stream  lunpri;       /* LUN for print output */
{
    char *input;
    string *inputs;
    int iret, i, j, n, nfixed, fixed, ninputs;
    real center[2], toarcsec, tokms;
    stream velstr, denstr, wtmapstr;
    string *burststring(), *fmode;
    bool Qdens, scanopt();
    bool Qwtmap;
    real *coldat[4];
    int colnr[4];

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
    Qwtmap = hasvalue("wtmap");
    Qdens = hasvalue("dens");                /* check if density given */
    if (Qimage) {                          /* get data from an image */
      read_image(velstr,&velptr);                 /* get data */
      copy_image(velptr,&resptr);
      if (Qreuse) {
	create_image_mask(velptr,&maskptr);
	for (i=0; i<Nx(maskptr); i++)
	for (j=0; j<Ny(maskptr); j++)
	  MapValue(maskptr,i,j) = TRUE;
      }
    } else {                                    /* get data from a table */
      n_vel = nemo_file_lines(input,100000);
      xpos_vel = (real *) allocate(n_vel * sizeof(real));
      ypos_vel = (real *) allocate(n_vel * sizeof(real));
      vrad_vel = (real *) allocate(n_vel * sizeof(real));
      if (Qdens) verr_vel = (real *) allocate(n_vel * sizeof(real));
      colnr[0] = 1;    coldat[0] = xpos_vel;
      colnr[1] = 2;    coldat[1] = ypos_vel;
      colnr[2] = 3;    coldat[2] = vrad_vel;
      if (Qdens) {
	colnr[3] = 4;    coldat[3] = verr_vel;
	n_vel = get_atable(velstr,4,colnr,coldat,n_vel);
	dprintf(1,"Found %d entries X,Y,V,DV in %s\n",input);
	for (i=0; i<n_vel; i++) {
	  dprintf(5,"%g %g %g %g\n",xpos_vel[i],ypos_vel[i],vrad_vel[i],verr_vel[i]);
	  verr_vel[i] = 1/sqr(verr_vel[i]);
	}
      } else {
	verr_vel = NULL;
	n_vel = get_atable(velstr,3,colnr,coldat,n_vel);
	dprintf(1,"Found %d entries X,Y,V in %s\n",input);
	for (i=0; i<n_vel; i++)
	  dprintf(5,"%g %g %g\n",xpos_vel[i],ypos_vel[i],vrad_vel[i]);
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
      
      printf("Mapsize is %g * %g arcsec; pixel %g arcsec; found %d/%d undefined map values\n",
	     ABS(dx*(lmax-lmin+1.0)), ABS(dy*(mmax-mmin+1.0)), ABS(dx),
	     n, Nx(velptr)*Ny(velptr));
      
      grid[0]=dx;        /* grid-separations in VELPAR (dx) */
      grid[1]=dy;        /*                            (dy) */
      pamp=0.0;          /* position angle of map */
    } else {
      printf("Mapsize unkown for point data, but toarcsec = %g\n",toarcsec);
      undf=getdparam("blank");        /* the undefined value */
      grid[0]=dx=1.0;
      grid[1]=dy=1.0;
      pamp=0.0;
    }
    if (Qwtmap) {
      input = getparam("wtmap");
      if (lunpri) fprintf(lunpri," weight map file      : %s\n", input);
      wtmapstr = stropen(input,"r");
      read_image(wtmapstr,&wtmapptr);
      strclose(wtmapstr);
      warning("Using special weight map for weights now");      
    } else
      wtmapstr = NULL;

    n = nemoinpr(getparam("beam"),beam,2);   /* get size of beam from user */
    if (n==2 || n==1) {       /* OK, got a beam, now get density map ... */
         if (n==1) beam[1] = beam[0];
         if (Qdens) {
            input = getparam("dens");
            if (lunpri) fprintf(lunpri," density file        : %s  beam: %g %g\n",
                                    input,beam[0],beam[1]);	
            denstr = stropen(input,"r");
    	    read_image(denstr,&denptr);
	    strclose(denstr);
	    if (!Qwtmap) 
	      warning("Using density map for weights now");
         } else {
            warning("beam defined, but no real beam correction used");
            if (lunpri) fprintf(lunpri,"  beam: %g %g\n",beam[0],beam[1]);
	    denptr = NULL;
         }
    } else {        /* no beam correction */
         beam[1] = beam[0] = 0.0;
         if (n!=0) warning("Parsing error beam=%s",getparam("beam"));
         printf("No beam correction\n");
         denstr = NULL;
	 denptr = NULL;
    }

    *nring = nemoinpr(getparam("radii"),rad,ring+1);
    if (*nring<2) error("radii=: Need at least two radii for one ring");
    *vsys = getdparam("vsys");
    n = nemoinpr(getparam("vrot"),vro,ring);
    if (n<1) error("vrot=: need at least one velocity (%d)",n);
    for (i=n;i<*nring;i++)
        vro[i] = vro[n-1];
    n = nemoinpr(getparam("pa"),pan,ring);
    if (n<1) error("vrot=: need at least one position angle (%d)",n);
    for (i=n;i<*nring;i++)
        pan[i] = pan[n-1];
    n = nemoinpr(getparam("inc"),inc,ring);
    if (n<1) error("vrot=: need at least one inclincation (%d)",n);
    for (i=n;i<*nring;i++)
        inc[i] = inc[n-1];
    n = nemoinpr(getparam("center"),center,2);
    if (n==2) {                     /* if two value supplied */
        *x0 = center[0];            /* this will be the center of rotation */
        *y0 = center[1];
    } else if (n==0) {                   /* if nothing supplied: */
        *x0 = 0.5*(lmin+lmax);      /* use center of map */
        *y0 = 0.5*(mmin+mmax);
    } else                          /* if all fails - barf */
        error("Need two numbers to define the center (%d)",n);
    *thf = getdparam("frang");

    printf("ROTCUR: free angle %4.1f (degrees)\n", *thf);
    if (lunpri) fprintf(lunpri," free angle          : %4.1f (degrees)\n",*thf);

#if 1
    input = getparam("side");
    if (*input=='r') {
        *side = 1;
        printf("ROTCUR: rotation curve of receding half\n");
        if (lunpri) fprintf(lunpri," velocity field      : receding half\n");
    } else if (*input=='a') {
        printf("ROTCUR: rotation curve of approaching half\n");
        if (lunpri) fprintf(lunpri," velocity field      : approaching side\n");
        *side = 2;
    } else if (*input=='w') {
        error("side=wedge not supported yet");
    } else {
        printf("ROTCUR: rotation curve of both halves\n");
        if (lunpri) fprintf(lunpri," velocity field      : both halves\n");
        *side = 3;
    }
#else
    iret = match(getparam("side"),"approaching,receding,both,wedge",side);
    switch (*side) {
     case 0x01:         /* approaching */
     case 0x02:         /* receding */
        break;
     case 0x04:         /* both */
        *side = 3;
        break;
     case 0x08:
        error("side=wedge not supported yet");
     default:
        error("Illegal side");
    }
#endif

    if (*side != 3) {  /* should we keep some parameters fixed ? */
         mask[0]=0;         /* fixed: systemic velocity */
         mask[4]=0;         /* fixed: x-position of center */
         mask[5]=0;         /* fixed: y-position of center */
    }
#if 1
    input = getparam("weight");
    if (*input == 'u') {
        printf("ROTCUR: uniform weighting\n");
        if (lunpri) fprintf(lunpri," used weighting      : uniform\n");
        *wpow = 0;
    } else if (streq(input,"cos-squared")) {
        printf("ROTCUR: cosine-squared weighting\n");
        if (lunpri) fprintf(lunpri," used weighting      : cos-squared\n");
        *wpow = 2;
    } else {            /* default is cosine weighting */
        printf("ROTCUR: cosine weighting\n");
        if (lunpri) fprintf(lunpri," used weighting      : cosine\n");
        *wpow = 1;
    }
#else
    iret=match(getparam("weight"),"uniform,cosine,cos-squared",wpow);
    select (*wpow) {
     case 0x01:     
        printf("ROTCUR: uniform weighting\n");
        if (lunpri) fprintf(lunpri," used weighting      : uniform\n");
        *wpow = 0;
        break;
     case 0x02:
        printf("ROTCUR: cosine weighting\n");
        if (lunpri) fprintf(lunpri," used weighting      : cosine\n");
        *wpow = 1;
        break;
     case 0x04:
        printf("ROTCUR: cosine-squared weighting\n");
        if (lunpri) fprintf(lunpri," used weighting      : cos-squared\n");
        *wpow = 2;
        break;
     default:       
        error("Illegal weight");
    }
#endif

#if 0
    input = getparam("fixed");
    if (scanopt(input,"vsys")) mask[0] = 0;
    if (scanopt(input,"vrot")) mask[1] = 0;
    if (scanopt(input,"pa"))   mask[2] = 0;
    if (scanopt(input,"inc"))  mask[3] = 0;        
    if (scanopt(input,"xpos")) mask[4] = 0;
    if (scanopt(input,"ypos")) mask[5] = 0;        
#else
    fixed = 0;
    iret=match(getparam("fixed"),"vsys,vrot,pa,inc,xpos,ypos",&fixed);
    if (iret<0) error("Illegal option in fixed=%s",getparam("fixed"));
    dprintf(1,"MASK: 0x%x ",fixed);
    for (i=0; i<6; i++) {
      mask[i] = (fixed & (1<<i)) ? 0 : 1;
      dprintf(1,"%d ",mask[i]);
    }
    dprintf(1,"\n");
#endif

    for (nfixed=0,i=0; i<PARAMS; i++)      /* count number of fixed par's */
        if (mask[i]==0) nfixed++; 

    printf("ROTCUR: will fit the following parameter(s)\n");
    if (lunpri) fprintf(lunpri," parameters to fit   :");
    if (mask[0]==1) {
        printf("       - Systemic velocity \n");
        if (lunpri) fprintf(lunpri," vsys");
    }
    if (mask[1]==1) {
        printf("       - Rotational velocity \n");
        if (lunpri) fprintf(lunpri," vrot");
    }
    if (mask[2]==1) {
        printf("       - Position angle \n");
        if (lunpri) fprintf(lunpri," pa");
    }
    if (mask[3]==1) {
        printf("       - Inclination \n");
        if (lunpri) fprintf(lunpri," inc");
    }
    if (mask[4]==1) {
        printf("       - X position of center \n");
        if (lunpri) fprintf(lunpri," xpos");
    }
    if (mask[5]==1) {
        printf("       - Y position of center \n");
        if (lunpri) fprintf(lunpri," ypos");
    }
    if (lunpri) fprintf(lunpri,"\n");
    if (nfixed==PARAMS) error("no free parameters left");

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
    printf("ROTCUR: tol=%g lab=%g iTmax=%d\n",tol,lab,itmax);
    if (lunpri) fprintf(lunpri," tol, lab, itmax     : %g %g %d\n",
                                    tol,lab,itmax);

    if (lunpri) fprintf(lunpri,"\n\n\n\n"); /* space space space... */

    *inh = (getbparam("inherit") ? 1 : 0);

    fmode = burststring(getparam("fitmode")," ,");
    if (fmode[0]==NULL || fmode[1]==NULL) 
        error("fitmode: %s needs {cos|sin},{#}",getparam("fitmode"));
    *fitmode = atoi(fmode[1]);
    switch (*fmode[0]) {
      case 's': *fitmode *= -1; break;
      case 'c': break;
      default: warning("Illegal fitmode: cos(%d*theta) assumed\n",*fitmode);
               break;
    }

    switch (*fitmode) {
      case 1:
        vobs  = vobs_c1;    vobsd = vobsd_c1;   vcor  = vcor_c1;    
        break;
      case -1:
        vobs  = vobs_s1;    vobsd = vobsd_s1;   vcor  = vcor_s1;
        break;
      default:
        error("fitmode %d not supported",*fitmode);
    }

    *nsigma = getdparam("nsigma");

    Qreuse = getbparam("reuse");
    Qwwb73 = getbparam("wwb73");

    if (!Qreuse)
      warning("New feature: fitted points will not be reused: reuse=f");
}

/*
 *    ROTFIT does a  least  squares fit to the radial velocity
 *    field. It can be called in two modes: nsigma=0 and nsigma>0
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
 *    FITMODE  integer         fitmode cos/sin(#)
 *
 *  returns:
 *    IER      integer         result of fit (good or bad)
 *    NPT      integer         number of points in ring
 *    RMS      real            residual rms velocities
 */

int rotfit(ri, ro, p, e, mask, wpow, side, thf, elp4, cor, npt, rms, fitmode, nsigma, useflag, lunres)
real ri,ro;      /* inner and outer radius of ring */
int mask[];      /* mask for free/fixed parameters */
int wpow;        /* weighting function */
int side;        /* half of galaxy */
int cor[];       /* print error ellipses ? */
real elp4[];     /* array for ellipse parameters */
real p[],e[];    /* initial estimates (I), results and erros (O) */
real thf;        /* free angle around minor axis (I) */
int *npt;	 /* number of points in ring (0) */
real *rms;       /* rms velocity in ring */
int fitmode;     /* basic fitmode (I) */
real nsigma;     /* if positive, remove outliers and fit again */
bool useflag;    /* flag: if TRUE, flag all ring points to undf, no fitting */
stream lunres;   /* file for residuals */
{
    int ier;                                             /* error return code */
    bool  stop,full;           /* booleans for stop fitting and data overflow */
    int   nfr;                                   /* number of free parameters */
    int   h, k;                                              /* loop counters */
    int   nrt;                       /* error return code from subroutine FIT */
    real  x[2*MAXPTS],y[MAXPTS],w[MAXPTS];/* arrays for coords, vels and wgts */
    int   iblank[MAXPTS];
    int   idx[2*MAXPTS];
    real  res[MAXPTS];                        /* array for y-x, the residuals */
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
   
    nfr=0;                                 /* reset number of free parameters */
    for (i=0; i<PARAMS; i++) {
        nfr += mask[i];                /* calculate number of free parameters */
        eps[i] = 0.1;                           /* and init the eps parameter */
    }
    r=0.5*(ri+ro);                                     /* mean radius of ring */

    printf(" radius of ring: %g \" \n",r); 
    printf("  iter.  systemic rotation position incli- ");
    printf("x-center y-center points  sigma\n");
    printf("  number velocity velocity   angle  nation ");
    printf("position position        velocity\n");

    getdat(x,y,w,idx,res,&n,MAXPTS,p,ri,ro,thf,wpow,&q,side,&full,nfr,useflag);  /* this ring */
    *rms = q;
    if (useflag)
      return -1;

    for (i=0;i<n;i++) iblank[i] = 1;

    h=0;                                           /* reset itegration counter */
    nblank=0;

    perform_out(h,p,n,q);                             /* show first iteration */
    if (Qwwb73) {
      return 1;
    }
    do {                                   /* and start the outer REPEAT loop */
         h++;                                               /* next iteration */
         chi=q;                                           /* save chi-squared */
         for(k=0;k<PARAMS;k++)              /* loop to save initial estimates */
            pf[k]=p[k];                                       /* save p in pf */
         nrt = nllsqfit(x,2,y,w,res,n,pf,e,mask,PARAMS,tol,itmax,lab,
                        vobs,vobsd);
	 
	 if (nsigma > 0 && h==1) {        /* first time only, with sigma > 0  */
	   stat2(res,n,&resmean,&ressig);
	   dprintf(1,"Refit: mean=%g sigma=%g nsigma=%g\n",resmean,ressig,nsigma);
	   if (ABS(resmean) > ressig) 
	     warning("Residuals not gaussian? mean=%g/sigma=%g => %g",
		     resmean,ressig,ABS(resmean)/ressig);
	   nblank = 0;
	   for (i=0;i<n;i++) {
	     ratio = ABS(res[i]-resmean)/ressig;
	     if (ratio > nsigma) {
	       iblank[i] = 0;
	       w[i] = 0.0;
	       nblank++;
	     } else {
	       iblank[i] = 1;
	     }
	     dprintf(4,"%d %g %g %g %g %s\n",i+1,y[i],w[i],res[i],
		     ratio,(ratio > nsigma) ? "***" : "" );
	   }
	   if (nblank) {
	     dprintf(0,"Rejecting %d/%d points for nsigma=%g\n",nblank,n,nsigma);
	     nrt = nllsqfit(x,2,y,w,res,n,p,e,mask,PARAMS,tol,itmax,lab,
                        vobs,vobsd);
	     if (nrt<0)
	       warning("nllsqfit=%d<0: PJT must find better solution (npoints=%d)",nrt,n);
	   }
	   break;
	 }
         if (nrt<0) {
            warning("nllsqfit=%d<0: KGB must find better solution (npoints=%d)",
		    nrt,n);
            break;    /* ???? could also try and continue here (like KGB did) */
	 }
         for (k=0; k<PARAMS; k++) 
            df[k]=pf[k]-p[k];                  /* calculate difference vector */
         flip=1.0;                                   /* factor for inner loop */
         for(;;) {                                        /* inner WHILE loop */
            stop = FALSE;
            for(k=0; k<PARAMS; k++)       /* loop to calculate new parameters */
               pf[k]=flip*df[k]+p[k];                       /* new parameters */
            pf[3]=MIN(pf[3],180.0-pf[3]);         /* in case inclination > 90 */
            getdat(x,y,w,idx,res,&n,MAXPTS,pf,ri,ro,thf,wpow,&q,side,&full,nfr,useflag);
           *rms = q;
	    for (i=0;i<n;i++) w[i] *= iblank[i];            /* apply blanking */
            if (q < chi) {                                     /* better fit ?*/
               perform_out(h,pf,n,q);                   /* show the iteration */
               for(k=0;k<PARAMS;k++)            /* loop to save new estimates */
                  p[k]=pf[k];
                stop=FALSE;  /* but make sure it doesn't quit from outer loop */
                dprintf(1,"Found a solution\n");
                break;                             /* leave inner XWHILE loop */
            } else {                                      /* not a better fit */
               stop=TRUE;                                    /* reset logical */
               for (k=0; k<PARAMS; k++)            /* loop through parameters */
                  stop = stop && (ABS(flip*df[k]) < eps[k]);
               stop = (ABS(q-chi)/chi < tol);      /* difference small enough */
               if (stop) {                                      /* is it so ? */
                  dprintf(1,"Chi^2 small enough anyhow....\n");
                  break;                             /* we found it ! XREPEAT */
               }
            }
	    /* with flip=1 initially the next if() makes this for() loop 
	     * go through flip=1,-1,1/2,-1/2,1/4,-1/4,......
	     * until convergence is reached 
	     */
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
         warning("ROTCUR: no degrees of freedom");
         break;
      case -2:
         warning("ROTCUR: you have no free parameters");
         break;
      case -3:                                        /* problems with matrix */
         warning("ROTCUR: not enough degrees of freedom");
         break;
      case -4:
         warning("ROTCUR: problems with matrix inversion");
         break;
      case -5:                                        /* floating zero divide */
         warning("ROTCUR: floating zero divide");
         break;
      case -10:                     /* maximum number of iterations too small */
         warning("ROTCUR: max. number of iterations %d to small",itmax);
         break;
      default:
         perform_out(h,p,n-nblank,q);                  /* write final results */
         if (full)
            warning("not all points inside ring %g were used",r);
    }

    if (lunres) {
      if (Qimage) {
	if (Qfirstring) {     /* for the first ring, reset the whole resid vel field */
	  for (i=0;i<Nx(resptr);i++)
	    for (j=0;j<Ny(resptr);j++)
	      MapValue(resptr,i,j) = 0.0;
	}
      } else {
        fprintf(lunres,"#  %d : New ring %g - %g\n",n,ri,ro);
        fprintf(lunres,"#  Xsky Ysky Vobs Vobs-Vmod Xgal Ygal Rgal THETAgal\n");
      }
      cosp = cos((p[2]+90)*F);
      sinp = sin((p[2]+90)*F);
      cosi = cos(p[3]*F);
      for (i=0; i<n; i++) {
	xc1 = x[2*i]-dx*p[4];
	yc1 = x[2*i+1]-dy*p[5];
	xc2 =   xc1*cosp + yc1*sinp;
	yc2 = (-xc1*sinp + yc1*cosp)/cosi;
	if (Qimage) {
	  MapValue(resptr,idx[2*i],idx[2*i+1]) = res[i];
	  if (i==0 && Qfirstring)
	    MapMin(resptr) = MapMax(resptr) = res[i];
	  else {
	    MapMin(resptr) = MIN(MapMin(resptr), res[i]);
	    MapMax(resptr) = MAX(MapMax(resptr), res[i]);
	  }
	  Qfirstring = FALSE;
	} else
	  fprintf(lunres,"%g %g %g %g %g %g %g %g\n",
		  xc1, yc1, y[i], res[i],
		  xc2, yc2, sqrt(xc2*xc2+yc2*yc2),atan2(yc2,xc2)/F);
      }
    }

    if (ier==1 && cor[0]>0 && cor[1]>0) {   /* calculate ellipse parameters ? */
         sigma2=0.0;                                      /* reset parameters */
         a11=0.0;
         a12=0.0;
         a22=0.0;
         for (i=0; i<n; i++) {                    /* loop through data points */
            vobsd(&x[2*i],p,b,PARAMS);            /* calculate derivatives */
            a11=a11+w[i]*b[cor[0]-1]*b[cor[0]-1];   /* calc ellips parameters */
            a22=a22+w[i]*b[cor[1]-1]*b[cor[1]-1]; 
            a12=a12+w[i]*b[cor[0]-1]*b[cor[1]-1]; 
            sigma2=sigma2+w[i]*sqr(y[i]-vobs(&x[2*i],p,PARAMS));
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

perform_out(int h,real *p,int n,real q)
{
/*  FORMAT(1H ,I4,4X,3(F7.2,2X),F5.2,2X,2(F7.2,2X),I5,2X,F8.3) eq.fortran */
    printf(" %4d    %7.2f  %7.2f  %7.2f  %5.2f  %7.2f  %7.2f  %5d  %8.3f\n",
              h,    p[0],   p[1],  p[2],  p[3],  p[4],  p[5],  n,   q);
}
/******************************************************************************/

/*
 * ROTPLT originally plots, but now only produces a table if lunpri
 *  is set. Otherwise it exits.
 */

void rotplt(rad,vsy,evs,vro,evr,pan,epa,inc,ein,xce,exc,yce,eyc,
       mask,ifit,elp,lunpri,cor,res,npt,factor)
real rad[],vsy[],evs[],vro[],evr[],pan[],epa[],inc[],ein[],
     xce[],exc[],yce[],eyc[],elp[][4], res[], factor;
int  mask[],ifit,cor[],npt[];
stream lunpri;
{
  int i;
  real mean, sig;

  if (lunpri==NULL) return;         /* determine if work to be done */
  if (ifit<0) return;               /* else return right now */

  fprintf(lunpri,"  radius   systemic   error  rotation   error");
  fprintf(lunpri," position    error   inclination  error   ");
  fprintf(lunpri,"x-position   error   y-position   error\n");

  fprintf(lunpri,"          velocity           velocity        ");
  fprintf(lunpri,"   angle               angle              ");
  fprintf(lunpri," of center            of center\n");

  fprintf(lunpri," (arcsec)   (km/s)   (km/s)   (km/s)    (km/s)");
  fprintf(lunpri,"(degrees)  (degrees)  (degrees)   (degrees) ");
  fprintf(lunpri,"(grids w.r.t. (0,0))  (grids w.r.t. (0,0))\n");
    
  for (i=0; i<ifit; i++) {
    fprintf(lunpri," %8.2f  %8.2f   %6.2f  %8.2f  %6.2f  %7.2f     ",
                         rad[i],vsy[i], evs[i],vro[i],evr[i],pan[i]);
    fprintf(lunpri," %5.2f     %6.2f     %6.2f   %8.2f  %8.2f   %8.2f  %8.2f",
                     epa[i],   inc[i],   ein[i], xce[i],exc[i], yce[i],eyc[i]);
    fprintf(lunpri," %d",npt[i]);
    fprintf(lunpri," %6.2f",res[i]);
    fprintf(lunpri,"\n");
  }
  fprintf(lunpri,"\n");
  fprintf(lunpri," Error correction factor: : %g\n",factor);
   
  stat2(inc,ifit,&mean,&sig);
  fprintf(lunpri," average inclination      : %8.2f  (%8.3f)  degrees\n",
                                              mean,    sig);
  stat2(pan,ifit,&mean,&sig);
  fprintf(lunpri," average position angle   : %8.2f  (%8.3f)  degrees\n",
                                              mean,    sig);
  stat2(vsy,ifit,&mean,&sig);
  fprintf(lunpri," average systemic velocity: %8.2f  (%8.3f)  km/s\n",
                                              mean,    sig);
  stat2(xce,ifit,&mean,&sig);
  fprintf(lunpri," average x-position       : %8.2f  (%8.3f)  grids\n",
                                              mean,    sig);
  stat2(yce,ifit,&mean,&sig);
  fprintf(lunpri," average y-position       : %8.2f  (%8.3f)  grids\n",
                                              mean,    sig);
}

/* stat2: for an input array a[] of length n, get the mean and dispersion */

void stat2(real *a,int n,real *mean,real *sig)
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
    *sig = n>1 ? sqrt((sxx-s*sqr(*mean)) / MAX(1.0,s-1))  :  0.0;
}

/* 
 *    GETDAT gets data from disk and calculates differences.
 *
 *    SUBROUTINE GETDAT(X,Y,W,IDX,RES,N,NMAX,P,RI,RO,THF,WPOW,Q,SIDE,FULL,NFR,UFLAG)
 *
 *    X        real array       sky coordinates of pixels inside ring
 *    Y        real array       radial velocities
 *    W        real array       weights of radial velocities
 *    IDX      integer array    back indexed into pixels (x,y)
 *    RES      real array       residuals
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
 *    UFLAG    logical          flag points in ring
 */

void getdat(x,y,w,idx,res,n,nmax,p,ri,ro,thf,wpow,q,side,full,nfr,useflag)
int   *n,nmax;       /* number of pixels loaded (O), max. number (I) */
int   nfr;           /* number of degrees of freedom */
int   wpow;          /* weigthing function */
int   side;          /* part of galaxy */
bool  *full;         /* output if data overflow */
real  x[],y[],w[];   /* arrays to store coords., vels. and weights */
int   idx[];
real  res[];
real  p[];           /* parameters of ring */
real  ri,ro;         /* inner and outer radius of ring */
real  thf;           /* free angle */
real  *q;             /* output sigma */
bool  useflag;
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

    phi=p[2]+pamp;        /* position angle plus map p.a. */
    inc=p[3];             /* inclination */
    x0=p[4];              /* x-position of center */
    y0=p[5];              /* y-position of center */
    free=ABS(sin(F*thf)); /* free angle in terms of sine */
    sinp=sin(F*phi);       /* sine of pa. */
    cosp=cos(F*phi);       /* cosine of pa. */
    sini=sin(F*inc);       /* sine of inc. */
    cosi=cos(F*inc);       /* cosine of inc. */
    a=sqrt(1.0-cosp*cosp*sini*sini);       /* find square around ellipse */
    b=sqrt(1.0-sinp*sinp*sini*sini);

    if (Qimage) {
      /* BUG:  need to fix this, there is a rounding problem near center */
      /*       although this is where rotcur isn't quite all that valid  */
      llo=MAX(lmin,(int)(x0-a*ro/dx));        /* square around ellipse */
      mlo=MAX(mmin,(int)(y0-b*ro/dy));
      lhi=MIN(lmax,(int)(x0+a*ro/dx));
      mhi=MIN(mmax,(int)(y0+b*ro/dy));
      if (llo > lhi || mlo > mhi) {
	error("ring not inside map [%d %d %d %d]",llo,mlo,lhi,mhi);
	*n = 0;                        /* continue here ??? */
        return;
      } else {
        dprintf(1,"getdat: P=%g %g %g %g %g %g\n",
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
	    /***** i=(m-m1)*nlt+l-lmin+1;        /* array pointer */
	    v = MapValue(velptr,l,m);        /* velocity at (l,m) */
	    if (v != undf) {       /* undefined value ? */
	      xr=(-(rx-dx*x0)*sinp+(ry-dy*y0)*cosp);     /* X position in galplane */
	      yr=(-(rx-dx*x0)*cosp-(ry-dy*y0)*sinp)/cosi;/* Y position in galplane */
	      r=sqrt(xr*xr+yr*yr);                       /* distance from center */
	      if (r < 0.01)                              /* radius too small ? */
		theta=0.0;
	      else
		theta=atan2(yr,xr)/F;
	      costh=ABS(cos(F*theta));       /* calculate |cos(theta)| */
	      dprintf(5,"@ %d,%d : r=%g cost=%g xr=%g yr=%g\n",l,m,r,costh,xr,yr);
	      if (r>ri && r<ro && costh>free) {     /* point inside ring ? */
		dprintf(5," ** adding this point\n");
		if (wtmapptr)
		  wi = MapValue(wtmapptr,l,m);
	        else if (denptr) 
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
		    idx[*n*2]=l;
		    idx[*n*2+1]=m;
		    s=(v-vobs(xx,p,PARAMS));  /* corrected difference */
		    res[*n] = s;
		    *q += s*s*wi;       /* calculate chi-squared */
		    *n += 1;           /* increment number of pixels */
		    if (useflag) {
		      MapValue(velptr,l,m) = undf;
		      continue;
		    }
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
    } else {                /* read from table instead of image */
      for (i=0; i<n_vel; i++) {
	rx = xpos_vel[i];
	ry = ypos_vel[i];
	v  = vrad_vel[i];
	wi = (verr_vel ? verr_vel[i] : 1.0);      /* weight */
	if (v == undf) continue;
	xr=(-(rx-x0)*sinp+(ry-y0)*cosp);     /* X position in galplane */
	yr=(-(rx-x0)*cosp-(ry-y0)*sinp)/cosi;/* Y position in galplane */
	r=sqrt(xr*xr+yr*yr);                       /* distance from center */
	if (r < 0.01)                              /* radius too small ? */
	  theta=0.0;
	else
	  theta=atan2(yr,xr)/F;
	costh=ABS(cos(F*theta));       /* calculate |cos(theta)| */
	if (r>ri && r<ro && costh>free) {     /* point inside ring ? */
	  dprintf(5,"@ r=%g cost=%g xr=%g yr=%g\n",r,costh,xr,yr);
	  dprintf(5," ** adding this point\n");
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
	      idx[*n*2] = i;
	      idx[*n*2+1] = -1;
	      s=(v-vobs(xx,p,PARAMS));  /* corrected difference */
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

    if (denptr) {       /* beam-correction wanted ? */
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

int  i,j;                                                    /* loop counters */
real vs,vc,phi=0.0,inc=0.0;                   /* parameters of velocity field */
real cosp1,cosp2,sinp1,sinp2,cosi1,cosi2,sini1,sini2;     /* different angles */
real x,y;                                                  /* sky coordinates */
real cost1,cost2,sint1,sint2,xx1,yy1,r,r1;  /* coordinates in plane of galaxy */
real fc,t[5],q[5],bx1,bx2,by1,by2;    /* vars for calculating beam-correction */
real vn,v2;                                  /* correction to radial velocity */

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
      return vs+vc*sini1*cost1;      /* circular motion */
}


void vobsd_c1(real *c,real *p,real *d,int m)
{
    perform_init(p,c);

    d[0]=1.0;                    /* partial derivative with respect to Vsys */
    d[1]=sini1*cost1;                                               /* Vrot */
    d[2]=F*vc*(cosi2+sini2*cost2)*sint1*sini1/cosi1;                  /* PA */
    d[3]=F*vc*(cosi2-sini2*sint2)*cost1/cosi1;                       /* Inc */
    d[4]= grid[0]*vc*(sint1*sinp1-cost1*cosp1/cosi1)*sint1*sini1*r1;  /* X0 */
    d[5]=-grid[1]*vc*(sint1*cosp1+cost1*sinp1/cosi1)*sint1*sini1*r1;  /* Y0 */
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
      perform_init(p,c);
      return vs+vc*sini1*sint1;      /* radial outflow motion */
}

void vobsd_s1(real *c,real *p,real *d,int m)
{
    perform_init(p,c);

    d[0]=1.0;                    /* partial derivative with respect to Vsys */
    d[1]=sini1*sint1;                                               /* Vrot */
    d[2]=F*vc*(sint2*sini2-1)*cost1*sini1/cosi1;                      /* PA */
    d[3]=F*vc*(1/cosi2-sini2*sint2)*sint1/cosi1;                     /* INC */
    d[4]= grid[0]*vc*(0);
    d[5]=-grid[1]*vc*(0);
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


perform_init(real *p,real *c)
{
    vs=p[0];                 /* systemic velocity */
    vc=p[1];                 /* circular velocity */
    if (p[2] != phi) {   /* new position angle ? */
       phi=p[2];               /* position angle */
       cosp1=cos(F*phi);       /* cosine */
       cosp2=cosp1*cosp1;      /* cosine squared */
       sinp1=sin(F*phi);       /* sine */
       sinp2=sinp1*sinp1;      /* sine squared */
    }
    if (p[3] != inc) {   /* new inclination ?  */
       inc=p[3];               /* inclination */
       cosi1=cos(F*inc);       /* cosine */
       cosi2=cosi1*cosi1;      /* cosine squared */
       sini1=sin(F*inc);       /* sine */
       sini2=sini1*sini1;      /* sine squared */
    }
    x=c[0]-p[4]*grid[0];          /* calculate x */
    y=c[1]-p[5]*grid[1];          /* calcualte y */
    xx1=(-x*sinp1+y*cosp1);        /* x in plane of galaxy */
    yy1=(-x*cosp1-y*sinp1)/cosi1;  /* y in plane of galaxy */
    r=sqrt(xx1*xx1+yy1*yy1);          /* distance from center */
    r1=1/r;                         /* save inverse for speed */
    cost1=xx1*r1;         /* cosine of angle in plane of galaxy */
    sint1=yy1*r1;         /* sine of angle in plane of galaxy */
    cost2=cost1*cost1;     /* cosine squared */
    sint2=sint1*sint1;     /* sine squared */
}


