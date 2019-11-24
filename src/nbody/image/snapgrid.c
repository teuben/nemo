/* 
 *  SNAPGRID:  program grids a snapshot into a 3D image-cube
 *             can do emission as well as absorption
 *             generalization of snapccd
 *
 *	19-jan-89  V1.0 -- derived from snapccd	PJT
 *      12-mar-89  V1.1 -- allow expression for emmisivity	PJT
 *                         allow taking mean instead of adding up
 *		       a  + warning when emissivity = 0
 *	 2-nov-90  V2.0 allow stacked snapshots
 *      21-oct-91  V3.0 allow moment=-1 or -2 for fast vel/sig gridding
 *                      also think about smoothing
 *	12-jun-92  V3.1 allow selection by time
 *      18-jul-92  V3.2 fixed bug for moment<0 and stack=t
 *                      floor() in gridding functions
 *      30-jul-93  V4.0 write out multiple images from multiple evar's
 *                      and using the new Unit's member of Image's
 *      26-jan-94  adding the new tvar= variable for optical depth (tau)
 *	22-jan-95  V4.2 fixed up the mean=t code
 *                      local variables
 *      11-jul-96  V4.3 finite size particles (svar) or 'variable smoothing'
 *      23-nov-96      a   prototypes, now works on linux
 *	 4-mar-97      b   support for SINGLEPREC
 *      18-jun-98   4.4 added optional overridden {x,y,z}lab's
 *      30-jul-98      a  fixed gridding bug introduced in 4.4
 *       7-oct-02      c  atof -> natof
 *       7-may-04   5.0 sky projections using proj= and new Image(5NEMO) format
 *       7-feb-06   5.1 integrate= option to compact the 3rd axis on the fly (ccdmom-like)
 *       9-feb-06   5.2 make it work for stack=t
 *       2-mar-11   5.3 implemented h3,h4 as moment -3 and -4
 *      18-may-12   5.4 added smoothing in VZ (szvar)
 *     13-feb-2013  6.0 units changed on a cube (now density instead of surface brightness?)
 *
 * Todo: - mean=t may not be correct for nz>1 
 *       - hermite h3 and h4 for proper kinemetry
 *         Gerhard(1993), v d Marel & Franx (1993)   
           0807.pdf?
 */

#include <stdinc.h>		/* nemo general */
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <strlib.h>

#include <snapshot/body.h>      /* snapshot's */
#include <snapshot/snapshot.h>
#include <snapshot/get_snap.c>

#include <image.h>              /* images */

string defv[] = {		/* keywords/default values/help */
	"in=???\n			  input filename (a snapshot)",
	"out=???\n			  output filename (an image)",
	"times=all\n			  standard selection of times",
	"xrange=-2:2\n			  x-range in gridding",
	"yrange=-2:2\n	  	 	  y-range in gridding",
	"zrange=-infinity:infinity\n	  z-range in gridding",
	"xvar=x\n			  x-variable to grid",
	"yvar=y\n			  y-variable to grid",
	"zvar=-vz\n			  z-variable to grid",
        "evar=m\n                         emission variable(s)",
        "tvar=0\n                         absorbtion variable",
        "dvar=z\n                         depth variable w/ tvar=",
        "svar=\n                          smoothing size variable in XY",
	"szvar=\n                         smoothing size variable in ZVAR",
	"nx=64\n			  x size of image",
	"ny=64\n			  y size of image",
	"nz=1\n			  	  z size of image/cube",
	"xlab=\n                          Optional X label [xvar]",
	"ylab=\n                          Optional Y label [yvar]",
	"zlab=\n                          Optional Z label [zvar]",
	"moment=0\n			  moment in zvar (-2,-1,0,1,2...)",
	"mean=f\n			  mean (moment=0) or sum per cell",
	"stack=f\n			  Stack all selected snapshots?",
	"integrate=f\n                    Sum or Integrate along 'dvar'?",
	"proj=\n                          Sky projection (SIN, TAN, ARC, NCP, GLS, MER, AIT)",
	"VERSION=6.0a\n			  27-jul-2019 PJT",
	NULL,
};

string usage="grid a snapshot into a 2/3D image-cube";

string cvsid="$Id$";

#define HUGE      1.0e20        /* don't use INF, ccdfits writes bad headers */
#define TIMEFUZZ  0.000001
#define CUTOFF    4.0		/* cutoff of gaussian in terms of sigma */
#define MAXVAR	  16		/* max evar's */

local stream  instr, outstr;				/* file streams */

		/* SNAPSHOT INTERFACE */
local int    nobj, bits, nbody=0, noutxy=0, noutz=0, nzero=0;
local real   tnow;
local Body   *btab = NULL;
local string times;

		/* IMAGE INTERFACE */
local imageptr  iptr=NULL, iptr0=NULL, iptr1=NULL, iptr2=NULL, iptr3=NULL, iptr4=NULL;
local int    nx,ny,nz;			     /* map-size */
local real   xrange[3], yrange[3], zrange[3];      /* range of cube */
local real   xbeam, ybeam, zbeam;                  /* >0 if convolution beams */
local int    xedge, yedge, zedge;                  /* check if infinite edges */

local string xvar, yvar, zvar;  	/* expression for axes */
local string xlab, ylab, zlab;          /* labels for output */
local rproc  xfunc, yfunc, zfunc;	/* bodytrans expression evaluator for axes */
local string *evar;
local rproc  efunc[MAXVAR];
local int    nvar;			/* number of evar's present */
local string dvar, tvar, svar, szvar;
local rproc  dfunc, tfunc, sfunc, szfunc;

local int    moment;	                /* moment to take in velocity */
local real   zsig;			/* positive if convolution in Z used */
local real   zmin, zmax;

local bool   Qmean;			/* adding or taking mean for image ?? */
local bool   Qint;                      /* integrate, instead of sum, along the 3rd dimension */
local bool   Qstack;                    /* stacking snapshots ?? */
local bool   Qdepth;                    /* need dfunc/tfunc for depth integration */
local bool   Qsmooth;                   /* (variable) smoothing */
local bool   Qzsmooth;                  /* (variable) smoothing */

local bool   Qwcs;                      /* use a real astronomical WCS in "fits" degrees */
local string proj;         
local double xref, yref, xrefpix, yrefpix, xinc, yinc, rot;

extern string  *burststring(string,string);
extern rproc   btrtrans(string);


local void setparams(void);
local void compfuncs(void);
local int read_snap(void);
local int allocate_image(void);
local int clear_image(void);
local int bin_data(int ivar);
local int free_snap(void);
local void los_data(void);
local int rescale_data(int ivar);
local int xbox(real x);
local int ybox(real y);
local int zbox(real z);
local real odepth(real tau);
local int setaxis(string rexp, real range[3], int n, int *edge, real *beam);


void nemo_main (void)
{
    int i;
    
    setparams();                /* set from user interface */
    compfuncs();                /* get expression functions */
    allocate_image();		/* make space for image(s) */
    if (Qstack) clear_image();	/* clear the images */
    if (Qstack && Qdepth) 
      error("stack=t and depth analysis cannot be used together, use snapmerge first");
    while (read_snap())	{                   /* read next N-body snapshot */
        for (i=0; i<nvar; i++) {
            if (!Qstack) {
            	if (nvar>1) dprintf(0,"Gridding evar=%s\n",evar[i]);
		clear_image();
	    }
            bin_data(i);	            /* bin and accumulate */
            if (!Qstack) {                  /* if multiple images: */
	      if (Qdepth||Qint)
		los_data();
	      rescale_data(i);            /* rescale */
	      write_image(outstr,iptr);   /* and write them out */
	      if (i==0) reset_history();  /* clean history */
            }
        }
        free_snap();
    }
    if (Qstack) {
      if (Qdepth||Qint)
	los_data();
      rescale_data(0);                    /* and rescale before ... */
      write_image (outstr,iptr);	    /* write the image */
    }

    strclose(instr);
    strclose(outstr);
}


void wcs(real *x, real *y)
{
  double xpix, ypix;
  double xpos = *x,
         ypos = *y;
  if (xypix(xpos, ypos, xref, yref, xrefpix, yrefpix, xinc, yinc, rot, proj, &xpix, &ypix)) {
    dprintf(0,"xypix %s %15.10g %15.10g %15.10g %15.10g\n",proj,xpos,ypos,xpix,ypix);
    error("problem-1 with WCS conversion");
  }
  *x = xpix;
  *y = ypix;
  /* dprintf(0,"xypix %s %15.10g %15.10g %15.10g %15.10g\n",proj,xpos,ypos,xpix,ypix); */
}


void to_dms(double dval, int *dd, int *mm, double *ss)
{
  int sign = SGN(dval);
  dval = ABS(dval);
  *dd = (int) floor(dval);
  dval = (dval-(*dd))*60.0;
  *mm = (int) floor(dval);
  *ss = (dval-(*mm))*60.0;
  *dd *= sign;
}

void show_wcs(string id, real longitude, real latitude)
{
    int d1,m1,d2,m2;
    double s1,s2;
    double xpos = longitude, ypos = latitude;
    double xpix, ypix;

    if (worldpos(xpos, ypos, xref, yref, xrefpix, yrefpix, xinc, yinc, rot, proj, &xpix, &ypix))
      warning("problem-2 with WCS conversion");

    to_dms(xpix,&d1,&m1,&s1);
    to_dms(ypix,&d2,&m2,&s2);
    dprintf(0,"%03d:%02d:%06.3f %03d:%02d:%06.3f   %s\n",d1,m1,s1,d2,m2,s2,id);
}




void setparams()
{
    char  *cp;

    times = getparam("times");
    Qmean = getbparam("mean");
    Qstack = getbparam("stack");
    Qint = getbparam("integrate");

    if (Qint) {
      warning("new version with integrate=t, only for MWP");
      warning("only use moment=0 stack=f mean=f etc.");
    }

      
    nx = getiparam("nx");
    ny = getiparam("ny");
    nz = getiparam("nz");
    if (nx<1 || ny<1 || nz<1) error("Bad grid size (%d %d %d)",nx,ny,nz);

    setaxis(getparam("xrange"), xrange, nx, &xedge, &xbeam);
    if (xbeam > 0.0) error("No convolve in X allowed yet");

    setaxis(getparam("yrange"), yrange, ny, &yedge, &ybeam);
    if (ybeam > 0.0) error("No convolve in Y allowed yet");

    Qwcs = hasvalue("proj");
    if (Qwcs) {  /* X and Y are now interpreted as a LONGITUDE / LATITUDE paid, in degrees */
      proj = getparam("proj");
      if (*proj != '-')
	error("proj=%s should be 4 uppercase letters starting with -, e.g. proj=-TAN",proj);
      rot = 0.0;                   /* we can only handle non-rotated frames */

      /* first save the astro WCS */
      /* note we're forcing the center pixel as the reference pixel */

      xrefpix = 0.5*(nx-1);
      yrefpix = 0.5*(ny-1);
      xref = 0.5*(xrange[0] + xrange[1]);
      yref = 0.5*(yrange[0] + yrange[1]);
      xinc = (xrange[1] - xrange[0])/nx;
      yinc = (yrange[1] - yrange[0])/ny;

      dprintf(0,"WCS: ref: %g %g refpix: %g %g  inc: %g %g proj: %s\n",
	      xref,yref,xrefpix,yrefpix,xinc,yinc,proj);
      

      /* and now replace the pixel coordinates */

      xrange[0] = 0.0;
      yrange[0] = 0.0;
      xrange[1] = nx-1.0;
      yrange[1] = ny-1.0;
      xrange[2] = (xrange[1] - xrange[0])/nx;
      yrange[2] = (yrange[1] - yrange[0])/ny;

      show_wcs("center pixel",xrefpix,yrefpix);
      show_wcs("lower left", xrange[0],yrange[0]);
      show_wcs("lower right",xrange[1],yrange[0]);
      show_wcs("upper right",xrange[1],yrange[1]);
      show_wcs("upper left", xrange[0],yrange[1]);


    }

    setaxis(getparam("zrange"), zrange, nz, &zedge, &zbeam);
    if (zbeam > 0.0) {                          /* with convolve */
        zsig = zbeam;                           /* beam */
        zmin = zrange[0] - CUTOFF*zsig;         /* make edges wider */
        zmax = zrange[1] + CUTOFF*zsig;
    } else {                                    /* straight gridding */
        zsig = 0.0;                             /* no beam */
        zmin = zrange[0];                       /* exact edges */
        zmax = zrange[1];
    }
    if (zedge && zsig>0.0)
        error("Cannot convolve in Z with an infinite edge");
    if (zedge && nz>1)
        error("Cannot use multiple Z planes with an infinite edge");

    zrange[2] = (zrange[1]-zrange[0])/nz;   /* reset grid spacing */
    dprintf (1,"size of IMAGE cube = %d * %d * %d\n",nx,ny,nz);

    xvar = getparam("xvar");
    yvar = getparam("yvar");
    zvar = getparam("zvar");
    evar = burststring(getparam("evar"),",");
    nvar = xstrlen(evar,sizeof(string)) - 1;
    dvar = getparam("dvar");
    tvar = getparam("tvar");
    Qsmooth = hasvalue("svar");
    if (Qsmooth) svar = getparam("svar");
    Qzsmooth = hasvalue("szvar");
    if (Qzsmooth) szvar = getparam("szvar");
    if (nvar < 1) error("Need evar=");
    if (nvar > MAXVAR) error("Too many evar's (%d > MAXVAR = %d)",nvar,MAXVAR);
    if (Qstack && nvar>1) error("stack=t with multiple (%d) evar=",nvar);
    xlab = hasvalue("xlab") ? getparam("xlab") : xvar;
    ylab = hasvalue("ylab") ? getparam("ylab") : yvar;
    zlab = hasvalue("zlab") ? getparam("zlab") : zvar;


    moment=getiparam("moment");
    if (moment<-4) error("Illegal moment=%d",moment);
    instr = stropen (getparam ("in"), "r");
    outstr = stropen (getparam("out"),"w");
}


void compfuncs()
{
    int i;
    
    xfunc = btrtrans(xvar);
    yfunc = btrtrans(yvar);
    zfunc = btrtrans(zvar);
    for (i=0; i<nvar; i++)
        efunc[i] = btrtrans(evar[i]);
    Qdepth = !streq(tvar,"0");
    if (Qdepth || Qint) {
        dfunc = btrtrans(dvar);
        tfunc = btrtrans(tvar);
    }
    if (Qsmooth)
        sfunc = btrtrans(svar);
    if (Qzsmooth)
        szfunc = btrtrans(szvar);
    
}

int read_snap()
{		
    for(;;) {		
        get_history(instr);
        get_snap(instr,&btab,&nobj,&tnow,&bits);
        if (bits==0) 
            break;           /* no snapshot at all */
        if ( (bits&PhaseSpaceBit)==0 )
            continue;       /* skip, no data, probably diagnostics */
        if (streq(times,"all") || within(tnow, times, TIMEFUZZ))
            break;          /* take it, if time in timerange */
    }
    if (bits) {
    	if (Qstack)
	    dprintf(1,"Adding data for times=%g; bits=0x%x\n",tnow,bits);
	else
	    dprintf(1,"New data for times=%g; bits=0x%x\n",tnow,bits);
    }
    return bits;
}

#define CV(i)  (CubeValue(i,ix,iy,iz))

allocate_image()
{
    create_cube (&iptr,nx,ny,nz);
    if (iptr==NULL) error("No memory to allocate first image");

    if (Qmean || moment <= -1) {
    	if (moment)
    	    warning("%d: mean=t requires moment=0",moment);
        create_cube(&iptr0,nx,ny,nz);
        if (iptr0==NULL) 
            error("No memory to allocate normalization image - try mean=f");
    }

    if (moment <= -1) {
        create_cube(&iptr1,nx,ny,nz);
        if (iptr1==NULL) 
            error("No memory to allocate 2nd image - try moment>0");
    }

    if (moment <= -2) {
        create_cube(&iptr2,nx,ny,nz);
        if (iptr2==NULL) 
            error("No memory to allocate 3rd image - try moment>0");
    }

    if (moment <= -3) {
        create_cube(&iptr3,nx,ny,nz);
        if (iptr3==NULL) 
            error("No memory to allocate 4th image - try moment>0");
    }

    if (moment <= -4) {
        create_cube(&iptr4,nx,ny,nz);
        if (iptr4==NULL) 
            error("No memory to allocate 5th image - try moment>0");
    }

    if (Qwcs) {
      Axis(iptr) = 1;

      Xmin(iptr) = xref;
      Ymin(iptr) = yref;
      Zmin(iptr) = zrange[0] + 0.5*zrange[2];
      
      Dx(iptr) = xinc;
      Dy(iptr) = yinc;
      Dz(iptr) = zrange[2];

      Xref(iptr) = xrefpix;
      Yref(iptr) = yrefpix;
      Zref(iptr) = 0.0;
      
      
      Namex(iptr) = hasvalue("xlab") ? xlab : sconc("RA--",proj);
      Namey(iptr) = hasvalue("ylab") ? ylab : sconc("DEC-",proj);
      Namez(iptr) = zlab;
      
    } else {
      Xmin(iptr) = xrange[0] + 0.5*xrange[2];
      Ymin(iptr) = yrange[0] + 0.5*yrange[2];
      Zmin(iptr) = zrange[0] + 0.5*zrange[2];
      
      Dx(iptr) = xrange[2];
      Dy(iptr) = yrange[2];
      Dz(iptr) = zrange[2];
      
      Namex(iptr) = xlab;
      Namey(iptr) = ylab;
      Namez(iptr) = zlab;
      
    }
    Beamx(iptr) = 0.0;  /* we don't allow smooth in the image plane for now */
    Beamy(iptr) = 0.0;
    Beamz(iptr) = (zsig>0.0 ? zsig : 0.0);
}

clear_image()
{
    int ix,iy,iz;

    for (ix=0; ix<nx; ix++)		        /* initialize all cubes to 0 */
    for (iy=0; iy<ny; iy++)
    for (iz=0; iz<nz; iz++) {
    	CV(iptr) = 0.0;
        if(iptr0) CV(iptr0) = 0.0;        
        if(iptr1) CV(iptr1) = 0.0;        
        if(iptr2) CV(iptr2) = 0.0;
        if(iptr3) CV(iptr3) = 0.0;
        if(iptr4) CV(iptr4) = 0.0;
    }
}


typedef struct point {
  real em, ab, z, depth;          /* emit, absorb, moment var, depth var */
  int i;                          /* particle id */
  struct point *next, *last;      /* pointers to aid */
} Point;

#define AddIndex(i,j,em,ab,z)  /**/

local Point **map =NULL;

local int pcomp(Point **a, Point **b);

bin_data(int ivar)
{
    real brightness, cell_factor, x, y, z, z0, t,sum;
    real expfac, fac, sfac, flux, b, emtau, depth;
    real e, emax, twosqs;
    int    i, j, k, ix, iy, iz, n, nneg, ioff;
    int    ix0, iy0, ix1, iy1, m, mmax;
    Body   *bp;
    Point  *pp, *pf,*pl, **ptab;
    bool   done;
    
    if (Qdepth || Qint) {
      /* first time around allocate a map[] of pointers to Point's */
        if (map==NULL)
            map = (Point **) allocate(Nx(iptr)*Ny(iptr)*sizeof(Point *));
        if (map==NULL || !Qstack) {
            for (iy=0; iy<Ny(iptr); iy++)
            for (ix=0; ix<Nx(iptr); ix++)
                map[ix+Nx(iptr)*iy] = NULL;
        }
    }
    expfac = (zsig > 0.0 ? 1.0/(sqrt(TWO_PI)*zsig) : 0.0);
    if (Qmean)
        cell_factor = 1.0;
    else {
        cell_factor = 1.0 / ABS(Dx(iptr)*Dy(iptr));
	if (Nz(iptr) > 1) cell_factor /= ABS(Dz(iptr));
    }

    if (Qint)
        cell_factor = 1.0;   

    nbody += nobj;
    if (Qsmooth) 
        mmax = MAX(Nx(iptr),Ny(iptr));
    else
        mmax = 1;
    emax = 10.0;

		/* big loop: walk through all particles and accumulate ccd data */
    for (i=0, bp=btab; i<nobj; i++, bp++) {
        x = xfunc(bp,tnow,i);            /* transform */
	y = yfunc(bp,tnow,i);
	if (Qwcs) wcs(&x,&y);            /* convert to an astronomical WCS, if requested */
        z = zfunc(bp,tnow,i);
        flux = efunc[ivar](bp,tnow,i);
        if (Qdepth || Qint) {
            emtau = odepth( tfunc(bp,tnow,i) );
            depth = dfunc(bp,tnow,i);
	}
        if (Qsmooth) {
            twosqs = sfunc(bp,tnow,i);
            twosqs = 2.0 * sqr(twosqs);
        }

	ix0 = xbox(x);                  /* direct gridding in X and Y */
	iy0 = ybox(y);


	if (ix0<0 || iy0<0) {           /* outside area (>= nx,ny never occurs */
	    noutxy++;
	    continue;
	}
        if (z<zmin || z>zmax) {         /* initial check in Z */
            noutz++;
            continue;
        }
        if (flux == 0.0) {              /* discard zero flux cases */
            nzero++;
            continue;
        }

      	dprintf(4,"%d @ (%d,%d) from (%g,%g)\n",
      	        i+1,ix0,iy0,x,y);

        for (m=0; m<mmax; m++) {        /* loop over smoothing area */
            done = TRUE;
            for (iy1=-m; iy1<=m; iy1++)
            for (ix1=-m; ix1<=m; ix1++) {       /* current smoothing edge */
                ix = ix0 + ix1;
                iy = iy0 + iy1;
        	if (ix<0 || iy<0 || ix >= Nx(iptr) || iy >= Ny(iptr))
        	    continue;

                if (m>0 && ABS(ix1) != m && ABS(iy1) != m) continue;
                if (m>0)
                    if (twosqs > 0)
                        e = (sqr(ix1*Dx(iptr))+sqr(iy1*Dy(iptr)))/twosqs;
                    else 
                        e = 2 * emax;
                else 
                    e = 0.0;
                if (e < emax) {
                    sfac = exp(-e);
                    done = FALSE;
                } else
                    sfac = 0.0;

                brightness =   sfac * flux * cell_factor;	/* normalize */
                b = brightness;
        	for (k=0; k<ABS(moment); k++) brightness *= z;  /* moments in Z */
                if (brightness == 0.0) continue;

                if (Qdepth || Qint) {	   /* stack away relevant particle info */
                    pp = (Point *) allocate(sizeof(Point));
                    pp->em = brightness;
                    pp->ab = emtau;
                    pp->z  = z;
		    pp->i  = i;
		    pp->depth = depth;
                    pp->next = NULL;
                    ioff = ix + Nx(iptr)*iy; /* location in grid map[] */
                    pf = map[ioff];
                    if (pf==NULL) {
                        map[ioff] = pp;
                        pp->last = pp;
                    } else {
                        pl = pf->last;
                        pl->next = pp;
                        pf->last = pp;
                    }

                    continue;       /* goto next accumulation now    CHECK */
                }


                if (zsig > 0.0) {           /* with Gaussian convolution in Z */
                    for (iz=0, z0=Zmin(iptr); iz<nz; iz++, z0 += Dz(iptr)) {
        	        fac = (z-z0)/zsig;
	        	if (ABS(fac)>CUTOFF)    /* if too far from gaussian center */
		            continue;           /* no contribution added */
        		fac = expfac*exp(-0.5*fac*fac);
                        CV(iptr) += brightness*fac;     /* moment */
			if(iptr0) CV(iptr0) += fac;     /* do we need that even here? */
                        if(iptr1) CV(iptr1) += b*fac;   /* moment -1,-2 */
                        if(iptr2) CV(iptr2) += b*z*fac; /* moment -2 */
                        if(iptr3) CV(iptr3) += b*z*z*fac; /* moment -3 */
                        if(iptr4) CV(iptr4) += b*z*z*z*fac; /* moment -4 */
        	    }
        	} else {
                    iz = zbox(z);
                    CV(iptr) +=   brightness;   /* moment */
                    if(iptr0) CV(iptr0) += 1.0; /* for mean */
                    if(iptr1) CV(iptr1) += b;   /* moment -1,-2 */
                    if(iptr2) CV(iptr2) += b*z; /* moment -2 */
                    if(iptr3) CV(iptr3) += b*z*z;   /* moment -3 */
                    if(iptr4) CV(iptr4) += b*z*z*z; /* moment -4 */
                }

            } /* for (iy1/ix1) */
            if (done) break;
        } /* m */
    }  /*-- end particles loop --*/
}


void los_data(void)
{
  real brightness, cell_factor, x, y, z, z0, t, dz, sum;
  real expfac, fac, sfac, flux, b, emtau, depth;
  real e, emax, twosqs;
  int    i, j, k, ix, iy, iz, n, nneg, ioff;
  int    ix0, iy0, ix1, iy1, m, mmax;
  Body   *bp;
  Point  *pp, *pf,*pl, **ptab;
  bool   done;
  int    nzero = 0, none = 0;

  if (Qdepth || Qint) {
    
    /* a hack, need to fix checking before we do this */
    if (iptr0) {
      iptr0 = 0; warning("iptr0 not 0, conflicting args");
    }
    if (iptr1) {
      iptr1 = 0; warning("iptr1 not 0, conflicting args");
    }
    if (iptr2) {
      iptr2 = 0; warning("iptr2 not 0, conflicting args");
    }
    if (iptr3) {
      iptr3 = 0; warning("iptr3 not 0, conflicting args");
    }
    if (iptr4) {
      iptr4 = 0; warning("iptr4 not 0, conflicting args");
    }

    
    /* for (i=0; i<Nx(iptr)*Ny(iptr); i++)  */
    /* accumulate */
    for (iy=0, i=0; iy<ny; iy++) {
      for (ix=0; ix<nx; ix++,i++) {
	if (map[i]==NULL) {
	  nzero++;
	  continue;
	}
	dprintf(1,"Id's along z:");
	
	n=1;                    /* count number along line of sight */
	pp = map[i];
	dprintf(1," %d (%g)",pp->i,pp->depth);
	while (pp->next) {
	  n++;
	  pp = pp->next;
	  dprintf(1," %d (%g)",pp->i,pp->depth);
	}
	dprintf(1,"\n");
	dprintf(1,"cell %d: %d\n",i,n);

	if (n == 1) {
	  none++;
	  continue;
	}
	
	/* allocate a huge array of pointers that we can mess (sort) with */
	ptab = (Point **) allocate (n * sizeof(Point *));
        
	n=1;
	pp = ptab[0] = map[i];
	while (pp->next) {
	  ptab[n++] = pp->next;
	  pp = pp->next;
	}

	/* sort the particles by decreasing 'depth' */
	qsort(ptab, n, sizeof(Point *), pcomp);
	
	for(j=0; j<n; j++)
	  dprintf(1,"depth(%d) = %d (%g)\n",j,ptab[j]->i,ptab[j]->depth);
	
	
	sum = 0.0;
	for (j=1; j<n; j++) {
	  dz = ptab[j]->depth - ptab[j-1]->depth;
	  if (dz == 0.0) continue;
	  sum += 0.5*(ptab[j]->em + ptab[j-1]->em) * dz;
	}
	dprintf(1,"sum=%g\n",sum);
	CubeValue(iptr,ix,iy,0) = sum;
      } /* ix */
    } /* iy */
    if (!Qstack) {          /* free the chain of visible particles */
      for (i=0; i<Nx(iptr)*Ny(iptr); i++) {
	pf = map[i];
	if (pf==NULL) continue;
	while(pf->next) {
	  pp = pf->next;
	  free(pf);
	  pf = pp;
	}
      } /* i */
    }
  }
  if (nzero)
    warning("%d pixels with no emission to integrate",nzero);
  if (none)
    warning("%d pixels with only 1 sample, set to 0, not enough to integrate",none);
    
}

int pcomp(Point **a, Point **b)
{
  return (*a)->depth > (*b)->depth;
}

free_snap()
{
    free(btab);         /* free snapshot */
    btab = NULL;        /* and make sure it can realloc at next get_snap() */
}    

rescale_data(int ivar)
{
  real m_min, m_max, brightness, total, x, y, z, b, mean, sigma, h_factor;
  int    i, j, k, ix, iy, iz, nneg, ndata;

    dprintf(1,"rescale(%d)\n",ivar);
    m_max = -HUGE;
    m_min = HUGE;
    total = 0.0;
    ndata = 0;

    /* Add the variable 4th and 5th dimension coordinates */
    Unit(iptr) = evar[ivar]; /* for Qmean=t should use proper mean cell units */
    Time(iptr) = tnow;
    
    /* handle special cases when mean=t or moment=-1 or moment=-2 */
    if (iptr4) {                   /* moment = -4 : h4 */
      dprintf(1,"rescale(%d) iptr4\n",ivar);
      nneg=0;
      h_factor = 1.0/(8*sqrt(6));
      for (ix=0; ix<nx; ix++)         /* loop over whole cube */
        for (iy=0; iy<ny; iy++)
	  for (iz=0; iz<nz; iz++) {
            if(CV(iptr) == 0.0) continue;
	    mean = CV(iptr2)/CV(iptr1);
	    sigma = CV(iptr3)/CV(iptr1) - mean*mean;
	    if (sigma < 0.0) {
	      nneg++;
	      b = 0.0;
	    } else if (CV(iptr0) == 1) {
	      b = 0.0;
	    } else {
	      b = -3.0 + 
		((CV(iptr)-4*CV(iptr4)*mean+6*CV(iptr3)*mean*mean)/CV(iptr1) -
		 3*mean*mean*mean*mean) / (sigma*sigma);
	      dprintf(1,"ixyz4: %d %d  %g %g %g  %g\n",ix,iy,b,mean,sigma,CV(iptr0));
	    }
            CV(iptr) = b * h_factor;
	  }
      if(nneg>0) warning("%d pixels with sig^2<0",nneg);
#if 0
/* kurt: zeta_2 = mu_4 / mu_2^2      = 3 for a gaussian    */
/*       zeta_2 = 3 + 8 sqrt(6) h_4  approx */
    mean = sum1/sum0;
    sigma = sum2/sum0 - mean*mean;

    tmp = -3.0 +
           ((sum4-4*sum3*mean+6*sum2*mean*mean)/sum0 - 3*mean*mean*mean*mean) /
           (sigma*sigma);
    return tmp;

#endif
    } else if (iptr3) {                    /* moment = -3 : h4 */
      dprintf(1,"rescale(%d) iptr3\n",ivar);
      nneg=0;
      h_factor = 1.0/(4*sqrt(3));
      for (ix=0; ix<nx; ix++)         /* loop over whole cube */
        for (iy=0; iy<ny; iy++)
	  for (iz=0; iz<nz; iz++) {
            if(CV(iptr) == 0.0) continue;
	    mean = CV(iptr2)/CV(iptr1);
	    sigma = CV(iptr3)/CV(iptr1) - mean*mean;
            if(sigma<0.0) {
	      nneg++;
	      b = 0;
	    } else if (CV(iptr0) == 1) {
	      b = 0.0;
            } else {
	      sigma = sqrt(sigma);
	      b = ((CV(iptr)-3*CV(iptr3)*mean)/CV(iptr1) + 2*mean*mean*mean) /
		(sigma*sigma*sigma);
	      dprintf(1,"ixyz4: %d %d  %g %g %g  %g\n",ix,iy,b,mean,sigma,CV(iptr0));
	    }
	    CV(iptr) = b * h_factor;
	  }
      if(nneg>0) warning("%d pixels with sig^2<0",nneg);
#if 0
/* skew: zeta_1 = mu_3 / mu_2^(3/2)  = 0 for a gaussian    */
/*       zeta_1 = 4 sqrt(3) h_3      approx */
    mean = sum1/sum0;
    sigma = sum2/sum0 - mean*mean;
    if (sigma < 0.0) sigma = 0.0;
    sigma = sqrt(sigma);    
    tmp = ((sum3-3*sum2*mean)/sum0 + 2*mean*mean*mean) /
            (sigma*sigma*sigma);
    return tmp;

#endif
    } else if (iptr2) {                    /* moment = -2 : dispersion output */
        dprintf(1,"rescale(%d) iptr2\n",ivar);
        nneg=0;
        for (ix=0; ix<nx; ix++)         /* loop over whole cube */
        for (iy=0; iy<ny; iy++)
        for (iz=0; iz<nz; iz++) {
            if(CV(iptr) == 0.0) continue;
	    mean = CV(iptr2)/CV(iptr1);
            sigma = CV(iptr)/CV(iptr1)-mean*mean;
            if(sigma < 0.0) {
	      nneg++;
	      b = 0.0;
            } else {
	      b = sqrt(sigma);
	      dprintf(1,"ixyz4: %d %d  %g %g %g  %g\n",ix,iy,b,mean,sigma,CV(iptr));
	    }
            CV(iptr) = b;
        }
        if(nneg>0) warning("%d pixels with sig^2<0",nneg);
    } else if(iptr1) {
        dprintf(1,"rescale(%d) iptr1\n",ivar);
        for (ix=0; ix<nx; ix++)         /* loop over whole cube */
        for (iy=0; iy<ny; iy++)
        for (iz=0; iz<nz; iz++) {
            if(CV(iptr) == 0.0) continue;       /* no emission */
            CV(iptr) /= CV(iptr1);
        }
    } else if(iptr0) {
        dprintf(1,"rescale(%d) iptr0\n",ivar);
        for (ix=0; ix<nx; ix++)         /* loop over whole cube */
        for (iy=0; iy<ny; iy++)
        for (iz=0; iz<nz; iz++) {
            if(CV(iptr) == 0.0) continue;       /* no emission */
            CV(iptr) /= CV(iptr0);
        }
        
    }

    for (ix=0; ix<nx; ix++)         	/* determine maximum in picture */
    for (iy=0; iy<ny; iy++)
    for (iz=0; iz<nz; iz++) {
       	  brightness = CV(iptr);
	  total += brightness;
	  m_max =  MAX(m_max,brightness);
	  m_min =  MIN(m_min,brightness);
	  if (brightness!=0.0)
	  	ndata++;
    }
   
    MapMin(iptr) = m_min;               /* min and max of data */
    MapMax(iptr) = m_max;
    BeamType(iptr) = NONE;              /* no smoothing yet */

    dprintf (1,"Total %d particles within grid\n",nbody-noutxy-noutz);
    dprintf (1,"     (%d were outside XY range, %d were not in Z range)\n",
                    noutxy, noutz);
    dprintf (1,"%d cells contain non-zero data,  min and max in map are %f %f\n",
    		ndata, m_min, m_max);
    dprintf (1,"Total mass in map is %f\n",total*Dx(iptr)*Dy(iptr));
    if (nzero)
        warning("There were %d stars with zero emissivity in the grid",nzero);
}

/*	compute integerized coordinates in X Y and Z, return < 0 if
 *	outside the range's - handle Z a little differently, as it
 *	allows edges at infinity
 */

int xbox(real x)
{
    if (xrange[0] < xrange[1]) {
        if (xrange[0]<x && x<xrange[1])
            return (int)floor((x-xrange[0])/xrange[2]);
    } else if (xrange[1] < xrange[0]) {
        if (xrange[1]<x && x<xrange[0])
            return (int)floor((x-xrange[0])/xrange[2]);
    } 
    return -1;
}

int ybox(real y)
{
    if (yrange[0] < yrange[1]) {
        if (yrange[0]<y && y<yrange[1])
            return (int)floor((y-yrange[0])/yrange[2]);
    } else if (yrange[1] < yrange[0]) {
        if (yrange[1]<y && y<yrange[0])
            return (int)floor((y-yrange[0])/yrange[2]);
    } 
    return -1;
}

int zbox(real z)
{
    if (zedge==0x03)			/* Both edges at infinity */
        return 0;                       /* always inside */
    if (zedge==0x01 && z<zrange[1])	/* left edge at infinity */
        return 0;
    if (zedge==0x02 && z>zrange[0])	/* right side at infinity */
        return 0;
    if (zedge)
        return -1;                      /* remaining cases outside */
	
    return (int)floor((z-zrange[0])/zrange[2]);    /* simple gridding */
}

real odepth(real tau)
{
    return exp(-tau);
}

/*
 * parse an expression of the form beg:end[,sig] into
 * range[2]:   0 = beg
 *	       1 = end
 *	       2 = (end-beg)/n   if sig not present (end>beg)
 * beam:
 *	       sig          if sig (>0) present
 *             -1           if no beam
 * 
 */
setaxis (string rexp, real range[3], int n, int *edge, real *beam)
{
    char *cp;
    
    *edge = 0;                  /* no infinite edges yet */
    cp = strchr(rexp,':');
    if (cp==NULL)
        error("range=%s must be of form beg:end[,sig]",rexp);
    if (strncmp(rexp,"-inf",4)==0) {
        range[0] = -HUGE;
        *edge |= 0x01;              /* set left edge at inifinity */	
    } else {
        *cp = 0;
        range[0] = natof(rexp);
    }
    cp++;
    if (strncmp(cp,"inf",3)==0) {
        range[1] = HUGE;
	*edge |= 0x02;              /* set right edge at infinity */
    } else
        range[1] = natof(cp);
    range[2] = (range[1]-range[0])/(real)n;       /* step */
    if (range[2] < 0)
      warning("%s: This axis %d has negative step",rexp,n);
    cp = strchr(cp,',');
    if (cp)
        *beam = natof(++cp);                  /* convolution beam */
    else
        *beam = -1.0;                        /* any number < 0 */
}
