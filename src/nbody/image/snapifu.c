/* 
 *  SNAPIFU:   generate spectra at a grid of points
 *
 *	 8-apr-09  V1.0 - cloned off snapgrid             PJT
 *
 * Todo: - mean=t may not be correct for nz>1 
 *       - xgrid,ygrid should be from file?
 *       - is it IFS or IFU  ??
 */

#include <stdinc.h>		/* nemo general */
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <strlib.h>

#include <snapshot/body.h>      /* snapshot's */
#include <snapshot/snapshot.h>
#include <snapshot/get_snap.c>

#include <image.h>              /* images */

string defv[] = {		/* keywords/default values/help */
  "in=???\n			  input filename (a snapshot)",
  "out=???\n			  output filename (an image)",
  "times=all\n			  standard selection of times",
  "xgrid=0\n                      X grid position(s)",
  "ygrid=0\n                      Y grid position(s)",
  "size=0.1\n                     diameter of IFU fiber",
  "xvar=x\n			  x-variable to grid",
  "yvar=y\n			  y-variable to grid",
  "zvar=-vz\n			  z-variable to grid",
  "evar=m\n                       emission variable(s)",
  "tvar=0\n                       absorbtion variable",
  "dvar=z\n                       depth variable w/ tvar=",
  "zrange=-infinity:infinity\n	  z-range in gridding",
  "nz=1\n			  z size of image",
  "moment=0\n			  moment in zvar (-2,-1,0,1,2...)",
  "mean=f\n			  mean (moment=0) or sum per cell",
  "stack=f\n			  Stack all selected snapshots?",
  "VERSION=1.0\n		  8-apr-09 PJT",
  NULL,
};

string usage="take spectra from a snapshot at a set of specified grid points";

string cvsid="$Id$"; 

#define HUGE      1.0e20        /* don't use INF, ccdfits writes bad headers */
#define TIMEFUZZ  0.000001
#define CUTOFF    4.0		/* cutoff of gaussian in terms of sigma */
#define MAXVAR	  16		/* max evar's */
#define MAXPOINTS 1024          /* max # gridpoints */

local stream  instr, outstr;				/* file streams */

		/* SNAPSHOT INTERFACE */
local int    nobj, bits, nbody=0, noutxy=0, noutz=0, nzero=0;
local real   tnow;
local Body   *btab = NULL;
local string times;


		/* IMAGE INTERFACE for spectra */
local imageptr  iptr=NULL, iptr0=NULL, iptr1=NULL, iptr2=NULL;
local int    nx,ny,nz;    		/* spectrum-size */
local real   zrange[3];                 /* range of spectra */
local real   zbeam;                     /* >0 if convolution beams */
local int    zedge;                     /* check if infinite edges */

local real   xgrid[MAXPOINTS], ygrid[MAXPOINTS];   /* IFU grid points */
local real   size;                                 /* size of each point */
local int    ngridx, ngridy;                       /* number of gridpoints */

local string xvar, yvar, zvar;  	/* expression for axes */
local string xlab, ylab, zlab;          /* labels for output */
local rproc  xfunc, yfunc, zfunc;	/* bodytrans expression evaluator for axes */
local string *evar;
local rproc  efunc[MAXVAR];
local int    nvar;			/* number of evar's present */
local string dvar, tvar;
local rproc  dfunc, tfunc;

local int    moment;	                /* moment to take in velocity */
local real   zsig;			/* positive if convolution in Z used */
local real   zmin, zmax;

local bool   Qmean;			/* adding or taking mean for image ?? */
local bool   Qstack;                    /* stacking snapshots ?? */
local bool   Qdepth;                    /* need dfunc/tfunc for depth integration */


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


nemo_main ()
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
		clear_image();
            	if (nvar>1) dprintf(0,"Gridding evar=%s\n",evar[i]);
	    }
            bin_data(i);	            /* bin and accumulate */
            if (!Qstack) {                  /* if multiple images: */
	      if (Qdepth)
		los_data();
	      rescale_data(i);            /* rescale */
	      write_image(outstr,iptr);   /* and write them out */
	      if (i==0) reset_history();  /* clean history */
            }
        }
        free_snap();
    }
    if (Qstack) {
      if (Qdepth)
	los_data();
      rescale_data(0);                    /* and rescale before ... */
      write_image (outstr,iptr);	    /* write the image */
    }

    strclose(instr);
    strclose(outstr);
}



void setparams()
{
    char  *cp;

    times = getparam("times");
    Qmean = getbparam("mean");
    Qstack = getbparam("stack");
      
    ny = 1;
    nz = getiparam("nz");
    if (nz<1) error("Bad grid size (%d)",nz);

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
    if (nvar < 1) error("Need evar=");
    if (nvar > MAXVAR) error("Too many evar's (%d > MAXVAR = %d)",nvar,MAXVAR);
    if (Qstack && nvar>1) error("stack=t with multiple (%d) evar=",nvar);
    xlab = xvar;
    ylab = yvar;
    zlab = zvar;

    ngridx = nemoinpd(getparam("xgrid"),xgrid,MAXPOINTS);
    ngridy = nemoinpd(getparam("ygrid"),ygrid,MAXPOINTS);
    if (ngridx < 0) error("xgrid parsing error");
    if (ngridy < 0) error("ygrid parsing error");
    if (ngridx != ngridy) error("xgrid and ygrid need equal number");

    nx = ngridx;
    
    size = getdparam("size");

    moment=getiparam("moment");
    if (moment<-2) error("Illegal moment=%d",moment);
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
  if (Qdepth) {
    dfunc = btrtrans(dvar);
    tfunc = btrtrans(tvar);
  }
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

    if (Qmean) {
    	if (moment)
    	    warning("%d: mean=t requires moment=0",moment);
        create_cube(&iptr0,nx,ny,nz);
        if (iptr0==NULL) 
            error("No memory to allocate normalization image - try mean=f");
    }

    if (moment == -1 || moment == -2) {
        create_cube(&iptr1,nx,ny,nz);
        if (iptr1==NULL) 
            error("No memory to allocate second image - try moment>0");
    }

    if (moment == -2) {
        create_cube(&iptr2,nx,ny,nz);
        if (iptr2==NULL) 
            error("No memory to allocate third image - try moment>0");
    }

    Xmin(iptr) = 0.0;
    Ymin(iptr) = 0.0;
    Zmin(iptr) = zrange[0] + 0.5*zrange[2];
      
    Dx(iptr) = 1.0;
    Dy(iptr) = 1.0;
    Dz(iptr) = zrange[2];
      
    Namex(iptr) = xlab;
    Namey(iptr) = ylab;
    Namez(iptr) = zlab;
      
    Beamx(iptr) = 0.0; 
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
    real r2, r2max;
    int    i, j, k, l, ix, iy, iz, n, nneg, ioff;
    Body   *bp;
    Point  *pp, *pf,*pl, **ptab;

    r2max = 0.25 * size * size;    /* size is the diameter, we need r^2 here */
    
    if (Qdepth) {
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
    else
        cell_factor = 1.0 / (r2max*PI);

    nbody += nobj;
    
    for (i=0, bp=btab; i<nobj; i++, bp++) {    /* loop over all particles */
      x = xfunc(bp,tnow,i);                    /* transform */
      y = yfunc(bp,tnow,i);
      z = zfunc(bp,tnow,i);
      flux = efunc[ivar](bp,tnow,i);
      if (Qdepth) {
	emtau = odepth( tfunc(bp,tnow,i) );
	depth = dfunc(bp,tnow,i);
      }
      for (l=0; l<ngridx; l++) {                   /* loop over all grid positions */
	ix = l;
	r2 = sqr(x-xgrid[l])+sqr(y-ygrid[l]);
	if (r2>r2max)
	  continue;

        if (z<zmin || z>zmax) {         /* initial check in Z */
            noutz++;
            continue;
        }
        if (flux == 0.0) {              /* discard zero flux cases */
            nzero++;
            continue;
        }

	brightness =   flux * cell_factor;	/* normalize */
	b = brightness;
	for (k=0; k<ABS(moment); k++) brightness *= z;  /* moments in Z */
	if (brightness == 0.0) continue;

	if (Qdepth) {	   /* stack away relevant particle info */
	  pp = (Point *) allocate(sizeof(Point));
	  pp->em = brightness;
	  pp->ab = emtau;
	  pp->z  = z;
	  pp->i  = i;
	  pp->depth = depth;
	  pp->next = NULL;
	  ioff = ix + Nx(iptr)*iy;      /* location in grid map[] */
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

	if (zsig > 0.0) {                   /* with Gaussian convolve in Z */
	  for (iz=0, z0=Zmin(iptr); iz<nz; iz++, z0 += Dz(iptr)) {
	    fac = (z-z0)/zsig;
	    if (ABS(fac)>CUTOFF)            /* if too far from gaussian center */
	      continue;                     /* no contribution added */
	    fac = expfac*exp(-0.5*fac*fac);
	    CV(iptr) += brightness*fac;     /* moment */
	    if(iptr1) CV(iptr1) += b*fac;   /* moment -1,-2 */
	    if(iptr2) CV(iptr2) += b*z*fac; /* moment -2 */
	  }
	} else {                            /* straight box gridding */
	  iz = zbox(z);
	  CV(iptr) +=   brightness;         /* moment */
	  if(iptr0) CV(iptr0) += 1.0;       /* for mean */
	  if(iptr1) CV(iptr1) += b;         /* moment -1,-2 */
	  if(iptr2) CV(iptr2) += b*z;       /* moment -2 */
	}
      } /*-- (l) end grid loop --*/
    }  /*-- (i) end particles loop --*/
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

  error("code not adapted yet for IFU");

  if (Qdepth) {
    
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
    real m_min, m_max, brightness, total, x, y, z, b;
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
    if (iptr2) {                    /* moment = -2 : dispersion output */
        dprintf(1,"rescale(%d) iptr2\n",ivar);
        nneg=0;
        for (ix=0; ix<nx; ix++)         /* loop over whole cube */
        for (iy=0; iy<ny; iy++)
        for (iz=0; iz<nz; iz++) {
            if(CV(iptr) == 0.0) continue;
            b = CV(iptr)/CV(iptr1)-sqr(CV(iptr2)/CV(iptr1));
            if(b<0.0) {
                nneg++;
                b=0.0;
            }
            CV(iptr) = sqrt(b);
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

/*	compute integerized coordinates in Z, return < 0 if
 *	outside the range's - handle Z a little differently, as it
 *	allows edges at infinity
 */

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
