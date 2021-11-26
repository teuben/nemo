/* 
 *  SNAPMAP:   program grids body variable snapshot into a 2D image
 *
 *	20-jun-09  V1.0 -- derived from snapgrid - GalaxyMasses09 - PJT
 *                      quick and dirty: only mean and svar used
 *       7-may-11  V1.3 added cone smoothing
 *       8-may-11  V1.4 make evar=m the default
 *
 *   TODO:
 *     fix gaussian weighted (still has the old svar= code from snapgrid)
 *     implement linear, i.e. ccdintpol, code (see also miriad::varmaps)
 *     check multiple evar's
 *     check stacking
 *
 */

#include <stdinc.h>		/* nemo general */
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <strlib.h>
#include <extstring.h>

#include <snapshot/body.h>      /* snapshot's */
#include <snapshot/snapshot.h>
#include <snapshot/get_snap.c>

#include <image.h>              /* images */

string defv[] = {		/* keywords/default values/help */
  "in=???\n		  input filename (a snapshot)",
  "out=???\n		  output filename (an image)",
  "times=all\n		  standard selection of times",
  "xrange=-2:2\n	  x-range in gridding",
  "yrange=-2:2\n	  y-range in gridding",
  "xvar=x\n		  x-variable to grid",
  "yvar=y\n		  y-variable to grid",
  "evar=m\n               emission variable(s)",
  "svar=\n                smoothing size variable",
  "sfunc=gauss\n          smoothing function: gauss,cone",
  "nx=64\n		  x size of image",
  "ny=64\n		  y size of image",
  "xlab=\n                Optional X label [xvar]",
  "ylab=\n                Optional Y label [yvar]",
  "mode=mean\n            Mode: mean,linear",
  "stack=f\n		  Stack all selected snapshots?",
  "proj=\n                Sky projection (SIN,TAN,ARC,NCP,GLS,MER,AIT)",
  "emax=10\n              Maximum exp smoothing parameter",
  "VERSION=2.1\n	  23-sep-2021 PJT",
  NULL,
};

string usage="grid a snapshot into a 2D map";

string cvsid="$Id$";

#define TIMEFUZZ  0.000001
#define MAXVAR	  16		/* max evar's */

#define SMODE_GAUSS  1
#define SMODE_CONE   2

local stream  instr, outstr;				/* file streams */

		/* SNAPSHOT INTERFACE */
local int    nobj, bits, nbody=0, noutxy=0, nzero=0;
local real   tnow;
local Body   *btab = NULL;
local string times;

		/* IMAGE INTERFACE */
local imageptr  iptr=NULL, iptr0=NULL, iptr1=NULL, iptr2=NULL;
local int    nx,ny;			     /* map-size */
local real   xrange[3], yrange[3];           /* range of map */
local real   xbeam, ybeam;                   /* >0 if convolution beams */
local int    xedge, yedge;                   /* check if infinite edges */

local string xvar, yvar;          	/* expression for axes */
local string xlab, ylab;                /* labels for output */
local rproc  xfunc, yfunc;	        /* bodytrans expression evaluator for axes */
local string *evar;
local rproc  efunc[MAXVAR];
local int    nvar;			/* number of evar's present */
local string svar;
local rproc  sfunc;
local real   emax; 
local int    smode;                     /* 1=gauss 2=cone */

local int    moment;	                /* moment to take in velocity */

local bool   Qmean;			/* adding or taking mean for image ?? */
local bool   Qstack;                    /* stacking snapshots ?? */
local bool   Qsmooth;                   /* (variable) smoothing */
local bool   Qmap;                      /* make map of linked list to the points for interpolation */

local bool   Qwcs;                      /* use a real astronomical WCS in "fits" degrees */
local string proj;         
local double xref, yref, xrefpix, yrefpix, xinc, yinc, rot;

extern string  *burststring(string,string);
extern rproc   btrtrans(string);


local void setparams(void);
local void compfuncs(void);
local int read_snap(void);
local void allocate_image();
local void clear_image();
local void bin_data(int ivar);
local void free_snap();
local void rescale_data(int ivar);
local int xbox(real x);
local int ybox(real y);
local int zbox(real z);
local real odepth(real tau);
local void setaxis(string rexp, real range[3], int n, int *edge, real *beam);



void nemo_main()
{
    int i;
    
    setparams();                /* set from user interface */
    compfuncs();                /* get expression functions */
    allocate_image();		/* make space for image(s) */
    if (Qstack) clear_image();	/* clear the images */
    while (read_snap())	{                   /* read next N-body snapshot */
        for (i=0; i<nvar; i++) {
            if (!Qstack) {
            	if (nvar>1) dprintf(0,"Gridding evar=%s\n",evar[i]);
		clear_image();
	    }
            bin_data(i);	            /* bin and accumulate */
            if (!Qstack) {                  /* if multiple images: */
	      rescale_data(i);            /* rescale */
	      write_image(outstr,iptr);   /* and write them out */
	      if (i==0) reset_history();  /* clean history */
            }
        }
        free_snap();
    }
    if (Qstack) {
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
  if (xypix(xpos,ypos,xref,yref,xrefpix,yrefpix,xinc,yinc,rot,proj,&xpix,&ypix)) {
    dprintf(0,"xypix %s %15.10g %15.10g %15.10g %15.10g\n",proj,xpos,ypos,xpix,ypix);
    error("problem-1 with WCS conversion");
  }
  *x = xpix;
  *y = ypix;
  /* dprintf(2,"xypix %s %15.10g %15.10g %15.10g %15.10g\n",proj,xpos,ypos,xpix,ypix); */
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
  
  if (worldpos(xpos,ypos,xref,yref,xrefpix,yrefpix,xinc,yinc,rot,proj,&xpix,&ypix))
    warning("problem-2 with WCS conversion");

  to_dms(xpix,&d1,&m1,&s1);
  to_dms(ypix,&d2,&m2,&s2);
  dprintf(0,"%03d:%02d:%06.3f %03d:%02d:%06.3f   %s\n",d1,m1,s1,d2,m2,s2,id);
}




void setparams()
{
  char  *cp;

  times = getparam("times");
  cp = getparam("mode");
  Qmean = (*cp == 'm');   /* mean mode */
  Qmap  = (*cp == 'l');   /* linear mode */
  Qstack = getbparam("stack");
      
  nx = getiparam("nx");
  ny = getiparam("ny");
  if (nx<1 || ny<1) error("Bad grid size (%d %d)",nx,ny);
  
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

  dprintf (1,"size of IMAGE map = %d * %d\n",nx,ny);
  
  xvar = getparam("xvar");
  yvar = getparam("yvar");
  evar = burststring(getparam("evar"),",");
  nvar = xstrlen(evar,sizeof(string)) - 1;
  Qsmooth = hasvalue("svar");
  if (Qsmooth) {
      if (Qmap) {
	warning("Cannot smooth in mode=linear map mode, smooth disabled");
	Qsmooth = FALSE;
      } else 
	svar = getparam("svar");
  }
  emax = getdparam("emax");
  cp = getparam("sfunc");
  if (*cp == 'g')
    smode = SMODE_GAUSS;
  else if (*cp == 'c') {
    smode = SMODE_CONE;
    emax = 0.5;
  } else
    error("Bad sfunc=%s",cp);
  
  if (nvar < 1) error("Need evar=");
  if (nvar > MAXVAR) error("Too many evar's (%d > MAXVAR = %d)",nvar,MAXVAR);
  if (Qstack && nvar>1) error("stack=t with multiple (%d) evar=",nvar);
  xlab = hasvalue("xlab") ? getparam("xlab") : xvar;
  ylab = hasvalue("ylab") ? getparam("ylab") : yvar;
  
  instr = stropen (getparam ("in"), "r");
  outstr = stropen (getparam("out"),"w");
}


void compfuncs()
{
    int i;
    
    xfunc = btrtrans(xvar);
    yfunc = btrtrans(yvar);
    for (i=0; i<nvar; i++)
        efunc[i] = btrtrans(evar[i]);
    if (Qsmooth)
        sfunc = btrtrans(svar);
}

int read_snap()
{
  dprintf(0,"SnapRead\n");
  for(;;) {		
        get_history(instr);
        get_snap(instr,&btab,&nobj,&tnow,&bits);
        if (bits==0) 
            break;           /* no snapshot at all */
        if ( (bits&PhaseSpaceBit)==0 )
            continue;       /* skip, no data, probably diagnostics */
        if (streq(times,"all") || within(tnow, times, TIMEFUZZ)) {
	  dprintf(0,"Time: %g\n",tnow);
          break;          /* take it, if time in timerange */
	}
  }
  if (bits) {
    	if (Qstack)
	    dprintf(1,"Adding data for times=%g; bits=0x%x\n",tnow,bits);
	else
	    dprintf(1,"New data for times=%g; bits=0x%x\n",tnow,bits);
  }
  return bits;
}

#define CV(i)  (MapValue(i,ix,iy))

void allocate_image()
{
    create_image(&iptr,nx,ny);
    if (iptr==NULL) error("No memory to allocate first image");

    create_image(&iptr0,nx,ny);
    if (iptr0==NULL) 
      error("No memory to allocate normalization image - try mean=f");

    if (Qwcs) {
      Axis(iptr) = 1;

      Xmin(iptr) = xref;
      Ymin(iptr) = yref;
      
      Dx(iptr) = xinc;
      Dy(iptr) = yinc;

      Xref(iptr) = xrefpix;
      Yref(iptr) = yrefpix;
      
      
      Namex(iptr) = hasvalue("xlab") ? xlab : sconc("RA--",proj);
      Namey(iptr) = hasvalue("ylab") ? ylab : sconc("DEC-",proj);
      
    } else {
      Xmin(iptr) = xrange[0] + 0.5*xrange[2];
      Ymin(iptr) = yrange[0] + 0.5*yrange[2];
      
      Dx(iptr) = xrange[2];
      Dy(iptr) = yrange[2];
      
      Namex(iptr) = xlab;
      Namey(iptr) = ylab;
      
    }
    Beamx(iptr) = 0.0;  /* we don't allow smooth in the image plane for now */
    Beamy(iptr) = 0.0;
}

void clear_image()
{
  int ix,iy;

  for (ix=0; ix<nx; ix++)		        /* initialize all cubes to 0 */
    for (iy=0; iy<ny; iy++) {
      CV(iptr) = 0.0;
      if(iptr0) CV(iptr0) = 0.0;        
    }
}


typedef struct point {
  real x, y, e;                   /* point info */
  struct point *next, *last;      /* pointers to aid */
} Point;

#define AddIndex(i,j,em,ab,z)  /**/

local Point **map =NULL;


void bin_data(int ivar)
{
  real brightness, x, y, z, z0, t,sum;
  real expfac, fac, sfac, flux, emtau, depth;
  real e, twosqs, delta;
  int    i, j, k, ix, iy, iz, n, nneg, ioff;
  int    ix0, iy0, ix1, iy1, m, mmax;
  Body   *bp;
  Point  *pp, *pf,*pl;
  bool   done;
    
  if (Qmap) {
    /* first time around allocate a map[] of pointers to Point's */
    warning("new linear mapping mode");
    if (map==NULL)
      map = (Point **) allocate(Nx(iptr)*Ny(iptr)*sizeof(Point *));
    if (map==NULL || !Qstack) {
      for (iy=0; iy<Ny(iptr); iy++)
	for (ix=0; ix<Nx(iptr); ix++)
	  map[ix+Nx(iptr)*iy] = NULL;
    }
  }
  
  nbody += nobj;
  if (Qsmooth) 
    mmax = MAX(Nx(iptr),Ny(iptr));
  else
    mmax = 1;

		/* big loop: walk through all particles and accumulate data */
  for (i=0, bp=btab; i<nobj; i++, bp++) {
    x = xfunc(bp,tnow,i);            /* transform */
    y = yfunc(bp,tnow,i);
    if (Qwcs) wcs(&x,&y);            /* convert to an astronomical WCS, if requested */
    flux = efunc[ivar](bp,tnow,i);
    if (Qsmooth) {
      delta = sfunc(bp,tnow,i);
      twosqs = 2.0 * sqr(delta);
    }
    
    ix0 = xbox(x);                  /* direct gridding in X and Y */
    iy0 = ybox(y);
    
    
    if (ix0<0 || iy0<0) {           /* outside area (>= nx,ny never occurs) */
      noutxy++;
      continue;
    }
    
    if (flux == 0.0) {              /* discard zero flux cases */
      nzero++;
      continue;
    }
    
    dprintf(4,"%d @ (%d,%d) from (%g,%g)\n",i+1,ix0,iy0,x,y);

    if (Qmap) {
      dprintf(1,"i=%d %g %g %g\n",i,x,y,e);
      pp = (Point *) allocate(sizeof(Point));
      pp->x = x;
      pp->y = y;
      pp->e = flux;
      pp->next = NULL;
      ioff = ix0 + Nx(iptr)*iy0; /* location in grid map[] */
      pf = map[ioff];
      if (pf==NULL) {          /* first point in this cell */
	dprintf(1,"first point\n");
	map[ioff] = pp;
	pp->last = pp;
      } else {                 /* append to last point in this cell, adjust pointers */
	dprintf(1,"continuing point\n");
	pl = pf->last;
	pl->next = pp;
	pf->last = pp;
      }      
    } else {
      /*  loop around the point and see where it has to deposit emission
       * 'm' is counting the ring from the center
       *  m=0 center point, m=1 the 8 point ring around it, etc.etc.
       *  it loops until all points in a ring have not contributed anymore
       */
      for (m=0; m<mmax; m++) {        /* loop over smoothing area in m-rings */
	done = TRUE;                          /* claim no more rings needed */
	for (iy1=-m; iy1<=m; iy1++)           /* unless in this ring */
	  for (ix1=-m; ix1<=m; ix1++) {       /* it is proven else */
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
	    
	    if (e < emax) {                    /* here is one contributing */
	      if (smode == SMODE_GAUSS)
		sfac = exp(-e);
	      else if (smode == SMODE_CONE)
		sfac = 1-sqrt(e);
	      done = FALSE;                    /* set the loop to still continue */
	    } else
	      sfac = 0.0;

	    if (sfac == 0.0) continue;
	    
	    brightness =   sfac * flux;     /* normalize */
	    CV(iptr) +=   brightness;       /* mean */
	    if(iptr0) CV(iptr0) += sfac;    /* for svar smoothing */

	  } /* for (iy1/ix1) */
	if (done) break;                    /* if there were no contributing pixels */
      } /* m */                             /* in this m-ring, quit */
    }
  }  /*-- end particles loop --*/
}


void free_snap()
{
  free(btab);         /* free snapshot */
  btab = NULL;        /* and make sure it can realloc at next get_snap() */
}    

void rescale_data(int ivar)
{
  real m_min, m_max, brightness, total, x, y, z, b, sum0, sum1;
  int    i, j, k, ix, iy, iz, nneg, ndata;
  Point  *pp, *pf,*pl;

  dprintf(1,"rescale(%d)\n",ivar);
  m_max = -HUGE;
  m_min = HUGE;
  total = 0.0;
  ndata = 0;

  /* Add the variable 4th and 5th dimension coordinates */
  Unit(iptr) = evar[ivar]; /* for Qmean=t should use proper mean cell units */
  Time(iptr) = tnow;
  
  if (Qmap) {
    warning("rescale: Linear interpolating, but testing it as mean");
    for (ix=0; ix<nx; ix++)         /* loop over whole cube */
      for (iy=0; iy<ny; iy++) {
	pp = map[ix + Nx(iptr)*iy];
	if (pp) {
	  sum0 = 1.0;
	  sum1 = pp->e;
	  while (pp->next) {
	    pp = pp->next;
	    sum0 += 1.0;
	    sum1 += pp->e;
	  }
	  sum1 /= sum0;
	  m_max =  MAX(m_max,sum1);
	  m_min =  MIN(m_min,sum1);
	  CV(iptr) = sum1;
	  ndata++;
	} else
	  CV(iptr) = 0.0;
      }
  } else {
    if(iptr0) {
        dprintf(1,"rescale(%d) iptr0\n",ivar);
        for (ix=0; ix<nx; ix++)         /* loop over whole cube */
	for (iy=0; iy<ny; iy++) {
                if(CV(iptr) == 0.0) continue;       /* no emission */
            CV(iptr) /= CV(iptr0);
        }
    }

    for (ix=0; ix<nx; ix++)         	/* determine maximum in picture */
      for (iy=0; iy<ny; iy++) {
       	  brightness = CV(iptr);
	  total += brightness;
	  m_max =  MAX(m_max,brightness);
	  m_min =  MIN(m_min,brightness);
	  if (brightness!=0.0)
	  	ndata++;
      }
  }
   
  MapMin(iptr) = m_min;               /* min and max of data */
  MapMax(iptr) = m_max;
  BeamType(iptr) = NONE;              /* no smoothing yet */

  dprintf (1,"Total %d particles within grid\n",nbody-noutxy);
  dprintf (1,"     (%d were outside XY range)\n",
	   noutxy);
  dprintf (1,"%d cells contain non-zero data,  min and max in map are %f %f\n",
	   ndata, m_min, m_max);
  dprintf (1,"Total mass in map is %f\n",total*Dx(iptr)*Dy(iptr));
  if (nzero)
    warning("There were %d stars with zero emissivity in the grid",nzero);
}

/*	compute integerized coordinates in X Y and Z, return < 0 if
 *	outside the range's - handle Z a little differently, as it
 *	also allows edges at infinity
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
void setaxis (string rexp, real range[3], int n, int *edge, real *beam)
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
