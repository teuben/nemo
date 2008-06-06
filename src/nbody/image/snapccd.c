/* 
 *  SNAPCCD:  program makes a CCD frame from a snapshot
 *
 *	16-Jun-87  V1.0 adapted from slit : 1D version   PJT
 *	17-Jun-87       dynamic allocation of 'frame'	 PJT
 *	25-jun-87  V2.0 write out image for further processsing   PJT
 *	30-jun-87  V2.1 new filestructure for image(5) PJT
 *			oeps, bugje
 *	 6-jul-87  V2.3 roundoff controlled through floor()	PJT
 *	 8-jul-87  V2.4 correct implementation of cellposition   PJT
 *	20-aug-87  V3.0 better memory management in reading Nbody PJT
 *	 8-Mar-88  V3.1 added data history	PJT
 *	 1-jun-88  V4.0 rename name of file	PJT
 *			new filestruct; get/put_hist
 *	22-dec-88  V4.1 allow ranges in radial velocity space PJT
 *	18-jan-89  V4.2 compatible with 3D image structure PJT
 *	30-jan-89  V4.3 vel into Zmin	PJT
 *      16-mar-90  V4.4 made GCC happy, made helpvec  PJT
 *	13-nov-90  V4.5 new location of <snapshot.h>	PJT
 *       4-mar-97  V4.6 NEMO V2.x, fixes for SINGLEPREC pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <snapshot/snapshot.h>
#include <image.h>

string defv[] = {
	"in=???\n			input filename (a snapshot)",
	"out=???\n			output filename (an image)",
	"origin=0,0\n			origin of CCD on 'sky'",
	"size=4.0\n			size (will be square)",
	"cell=0.0625\n			cellsize : 64 pixels if size=4",
	"vrange=-infinity:infinity\n	range in velocity space",
	"moment=0\n			velocity moment to weigh with",
	"VERSION=4.6b\n			6-jun-08 PJT",
	NULL,
};

string usage = "simple conversion of snapshot to image";

#define HPI  1.5702
#define RPD (3.1415/360.0)
#define EPS  0.00001
#ifndef HUGE
# define HUGE 1.0e20
#endif

local string	infile, outfile;			/* file names */
local stream  instr, outstr;				/* file streams */

		/* SNAPSHOT INTERFACE */
local int    nobj;
local real   tnow;
local real   *mass=NULL;
local real   *phase=NULL;
local string headline;
		/* IMAGE INTERFACE */
local imageptr  iptr=NULL;		/* will be allocated dynamically */
local int    nx,ny,nsize;		/* map-size */
local real origin[2];			/* of map in sky coordinates */
local real size;			/* size of frame (square) */
local real cell;			/* cell or pixel size (square) */

local int    moment;	                /* moment to take in velocity */
local real vmin,vmax;			/* range in velocity space */
local real vmean, vsig;		        /* beam in velocity space */
local bool   vbeam;

local char *axname[] = {               /* axnames */
        "x", "y", "z", "vel-x", "vel-y", "vel-z"
    }; 


nemo_main ()
{
	setparams();                    /* stuff command line [pars] */

	instr = stropen (infile, "r");
	outstr = stropen (outfile,"w");

	read_snap();			/* read N-body data */
	allocate_image();		/* make space for image */
	bin_data();			/* do the heavy work */
	write_image (outstr,iptr);	/* write the image */

	strclose(instr);
	strclose(outstr);
}

setparams()
{
	string  tmpstr;
	int	tmpint;
	char    *cp;

	if (NDIM!=3)
		error ("This program requires 3D-data, 2D could work\n");

	infile =  getparam ("in");
	outfile = getparam("out");

        tmpint = nemoinpr(getparam("origin"),origin,2);
	if (tmpint!=2)
	  error ("Keyword 'origin' must have two values, e.g. origin=0,0\n");

	size = getdparam("size");
	cell = getdparam("cell");
	nsize = size/cell + EPS;    /* fudge factor to prevent unfavourable roundoff */
	dprintf (1,"size of CCD frame = %d square\n",nsize);

	moment=getiparam("moment");
	dprintf (1,"Vel_moment=%d\n",moment);
	if (moment<0)
		error("Vel-moments should be non-negative");

	tmpstr = getparam("vrange");
	cp = strchr(tmpstr,':');
	if (cp==NULL) {			/* beam specified:   vmean,vsig */
            cp = strchr(tmpstr,',');
            if (cp==NULL)
                error("syntax error vrange=vmin:vmax or vrange=vmean,vsig\n");
	    vmean = atof(tmpstr);
            vsig = atof(++cp);
            vmin = vmean - 3*vsig;
            vmax = vmean + 3*vsig;
            dprintf(2,"Gaussian beam specified: vmean=%f vsig=%f\n",vmean,vsig);
            vbeam = TRUE;
        } else {                    /* vrange specified: vmin:vmax */
            if (strncmp(tmpstr,"-infinity",9)==0)
                vmin = -HUGE;
            else
                vmin = atof(tmpstr);
            if (strncmp(++cp,"infinity",8)==0)
                vmax = HUGE;
            else
                vmax = atof(cp);
            dprintf(2,"vrange = %f : %f\n",vmin,vmax);
            vbeam = FALSE;
        }
}

read_snap()
{				
    int    i;
    real   imass;    

    get_history(instr);

    if (! get_tag_ok(instr, SnapShotTag))
    	return(0);				/* no snapshot found */
    	
    get_set(instr, SnapShotTag);
      get_set(instr, ParametersTag);
        get_data(instr, NobjTag, IntType, &nobj, 0);
	if (get_tag_ok(instr, TimeTag))
	   get_data(instr, TimeTag, RealType, &tnow, 0);
	else
	   tnow=0.0;
      get_tes(instr, ParametersTag);
      dprintf (1,"nobj=%d time=%f\n", nobj, tnow);
      get_set(instr, ParticlesTag);
         if (mass==NULL) {		/* first time around */
            mass = (real *) malloc(nobj * sizeof(real));
            phase = (real *) malloc(nobj * 2 * NDIM * sizeof(real));
            if (mass==NULL || phase==NULL)
                error ("read_snap: not enough memory for N-body snapshot\n");
         }
         if (get_tag_ok(instr,MassTag))
           get_data(instr, MassTag, RealType, mass, nobj, 0);
         else {
	   dprintf (0,"WARNING: no masses provided, will assume 1/%d\n",nobj);
	   imass = 1.0/(real)nobj;
	   for (i=0; i<nobj; i++)
	      mass[i] = imass;
	 }
         get_data(instr, PhaseSpaceTag, RealType, phase, nobj, 2, NDIM, 0);
      get_tes(instr, ParticlesTag);
    get_tes(instr,SnapShotTag);
    return(1);      /* make GCC happy */
}

allocate_image()
{
/*	iptr = (imageptr ) malloc(sizeof(image));   OLD STYLE */
        create_image (&iptr,nsize,nsize);       /* force square image */
	if (iptr==NULL) {
		dprintf (0,"Not enough memory for image\n");
		stop(0);
	}
	Nx(iptr) = Ny(iptr) = nsize;			/* image will be square */
	Xmin(iptr) = origin[0] - 0.5*size + 0.5*cell;   /* note definition of */
	Ymin(iptr) = origin[1] - 0.5*size + 0.5*cell;   /* cell position here */
	Dx(iptr) = Dy(iptr) = cell;
        Dz(iptr) = 0.0;
	Frame(iptr) = (real *) malloc(Nx(iptr)*Ny(iptr)*sizeof(real));
	if (Frame(iptr)==NULL) {
		dprintf (0,"Not enough memory to allocate %d*%d frame array\n",Nx(iptr),Ny(iptr));
		stop(1);
	}
        Namex(iptr) = axname[0];        /* may change */
        Namey(iptr) = axname[1];
        Namez(iptr) = axname[5];
}

bin_data()
{
    real xsky, ysky, vrad;
    real m_min, m_max, brightness, inv_surden, total;
    int  i, k, ix, iy, nx, ny, cnt, noutside, noutvel, ndata;
    real *pptr;
    
	/* (re)initialize CCD and some other local  variables */
    nx=Nx(iptr);
    ny=Ny(iptr);
    for (ix=0; ix<nx; ix++)
       for (iy=0; iy<ny; iy++)
	  MapValue(iptr,ix,iy) = 0.0;
    m_max = -HUGE;
    m_min = HUGE;
    inv_surden = 1.0 / (cell * cell);		/* scaling factor */
    noutside=noutvel=ndata=0;
    total=0.0;
		/* walk through all particles and accumulate ccd data */
    for (i=0, pptr=phase; i<nobj; i++) {
        xsky = *pptr++;			/* x */
        ysky = *pptr++;			/* y */
        pptr += NDIM;
        vrad = -(*pptr++);		/* v_z */
	ix = floor((xsky-Xmin(iptr))/cell+0.5+EPS);	/* integer coords in CCD */
	iy = floor((ysky-Ymin(iptr))/cell+0.5+EPS);
	if (ix<0 || iy<0 || ix>=nx || iy>=ny) {
	    noutside++;			/* not in CCD frame */
	    continue;
	}
        if (vrad<vmin || vrad>vmax) {
            noutvel++;
            continue;
        }
	brightness =   mass[i] * inv_surden;
        if (vbeam)
            brightness *= exp( -0.5*sqr((vrad-vmean)/vsig) );
	for (k=0; k<moment; k++)
	    brightness *= vrad;
	MapValue(iptr,ix,iy) +=   brightness;
     }  /*-- end particles loop --*/

	/* determine maximum in picture */
    for (ix=0; ix<nx; ix++)
       for (iy=0; iy<ny; iy++) {
       	  brightness = MapValue(iptr,ix,iy);
	  total += brightness;
	  if (brightness > m_max)
	  	m_max = brightness;
	  else if (brightness < m_min)
	  	m_min = brightness;
	  if (brightness!=0.0)
	  	ndata++;
	}

   
    MapMin(iptr) = m_min;               /* min and max of data */
    MapMax(iptr) = m_max;
    Beamx(iptr)=Beamy(iptr)= 0.0;       /* beam in position space */
    BeamType(iptr) = NONE;              /* no smoothing yet */
    if (vbeam) {
       Zmin(iptr) = vmean;               /* velocity of data */
       Beamz(iptr) = vsig;
    } else
       Zmin(iptr) = 0.5*(vmin+vmax);

    dprintf (1,"%d particles in CCD frame\n",nobj-noutside-noutvel);
    dprintf (1,"     (%d were outside frame, %d were not in velocity range)\n",
                    noutside, noutvel);
    dprintf (1,"%d cells contain non-zero data,  min and max in map are %f %f\n",
    		ndata, m_min, m_max);
    dprintf (1,"Total mass in map is %f\n",total*Dx(iptr)*Dy(iptr));

}


