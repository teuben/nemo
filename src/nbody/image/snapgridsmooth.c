/* 
 *  SNAPGRID:  program grids a snapshot into a 3D image-cube, with adaptive smoothing
 *             simplification of snapgrid
 *
 *	 1-nov-06  V0.1 -- derived from snapgrid	PJT/Ed shaya
 *
 */

#include <stdinc.h> 
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <strlib.h>

#include <snapshot/body.h>
#include <snapshot/snapshot.h>
#include <snapshot/get_snap.c>

#include <image.h>        

string defv[] = {		
  "in=???\n			  input filename (a snapshot)",
  "out=???\n			  output filename (an image cube)",
  "times=all\n			  standard selection of times",
  "xrange=-2:2\n		  x-range in gridding",
  "yrange=-2:2\n	  	  y-range in gridding",
  "zrange=-2:2\n	  	  z-range in gridding",
  "xvar=x\n			  x-variable to grid",
  "yvar=y\n			  y-variable to grid",
  "zvar=z\n			  z-variable to grid",
  "evar=m\n                       emission variable(s)",
  "svar=\n                        smoothing size variable",
  "nx=64\n			  x size of cube",
  "ny=64\n			  y size of cube",
  "nz=64\n		  	  z size of cube",
  "stack=f\n			  Stack all selected snapshots?",
  "VERSION=0.1\n		  1-nov-06 PJT",
  NULL,
};

string usage="grid a snapshot into a 3D image-cube with adaptive smoothing";

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
local imageptr  iptr=NULL, iptr0=NULL, iptr1=NULL, iptr2=NULL;
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
local string svar;
local rproc  dfunc, tfunc, sfunc;

local real   zsig;			/* positive if convolution in Z used */
local real   zmin, zmax;

local bool   Qstack;                    /* stacking snapshots ?? */
local bool   Qsmooth;                   /* (variable) smoothing */

/* local double xref, yref, xrefpix, yrefpix, xinc, yinc, rot; */

extern string  *burststring(string,string);
extern rproc   btrtrans(string);


local void setparams(void);
local void compfuncs(void);
local int read_snap(void);
local int allocate_image(void);
local int clear_image(void);
local int bin_data(int ivar);
local int free_snap(void);
local int rescale_data(int ivar);
local int xbox(real x);
local int ybox(real y);
local int zbox(real z);
local int setaxis(string rexp, real range[3], int n, int *edge, real *beam);


nemo_main ()
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



void setparams()
{
    char  *cp;

    times = getparam("times");
    Qstack = getbparam("stack");

    nx = getiparam("nx");
    ny = getiparam("ny");
    nz = getiparam("nz");
    if (nx<1 || ny<1 || nz<1) error("Bad grid size (%d %d %d)",nx,ny,nz);

    setaxis(getparam("xrange"), xrange, nx, &xedge, &xbeam);
    if (xbeam > 0.0) error("No convolve in X allowed yet");

    setaxis(getparam("yrange"), yrange, ny, &yedge, &ybeam);
    if (ybeam > 0.0) error("No convolve in Y allowed yet");

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
    Qsmooth = hasvalue("svar");
    if (Qsmooth) svar = getparam("svar");
    if (nvar < 1) error("Need evar=");
    if (nvar > MAXVAR) error("Too many evar's (%d > MAXVAR = %d)",nvar,MAXVAR);
    if (Qstack && nvar>1) error("stack=t with multiple (%d) evar=",nvar);
    xlab = xvar;
    ylab = yvar;
    zlab = zvar;

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
    if (Qsmooth)
        sfunc = btrtrans(svar);
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

    Xmin(iptr) = xrange[0] + 0.5*xrange[2];
    Ymin(iptr) = yrange[0] + 0.5*yrange[2];
    Zmin(iptr) = zrange[0] + 0.5*zrange[2];
    
    Dx(iptr) = xrange[2];
    Dy(iptr) = yrange[2];
    Dz(iptr) = zrange[2];
    
    Namex(iptr) = xlab;
    Namey(iptr) = ylab;
    Namez(iptr) = zlab;
    
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
    }
}


typedef struct point {
  real em, ab, z, depth;          /* emit, absorb, moment var, depth var */
  int i;                          /* particle id */
  struct point *next, *last;      /* pointers to aid */
} Point;

#define AddIndex(i,j,em,ab,z)  /**/

local Point **map =NULL;

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
    
    expfac = (zsig > 0.0 ? 1.0/(sqrt(TWO_PI)*zsig) : 0.0);
    cell_factor = 1.0 / ABS(Dx(iptr)*Dy(iptr));

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
        z = zfunc(bp,tnow,i);
        flux = efunc[ivar](bp,tnow,i);
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
                if (brightness == 0.0) continue;

                if (zsig > 0.0) {           /* with Gaussian convolve in Z */
                    for (iz=0, z0=Zmin(iptr); iz<nz; iz++, z0 += Dz(iptr)) {
        	        fac = (z-z0)/zsig;
	        	if (ABS(fac)>CUTOFF)    /* if too far from gaussian center */
		            continue;           /* no contribution added */
        		fac = expfac*exp(-0.5*fac*fac);
                        CV(iptr) += brightness*fac;     /* moment */
        	    }
        	} else {
                    iz = zbox(z);
                    CV(iptr) +=   brightness;   /* moment */
                }

            } /* for (iy1/ix1) */
            if (done) break;
        } /* m */
    }  /*-- end particles loop --*/
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
