/* 
 * CCDMATH: manipulate a series of maps to generate a new map,
 *		can also start from scratch
 *
 *	1-Jul-87	Original version PJT
 *	3-Jul-87	V1.1: order of keywords changed for future expansion
 *	1-jun-88	V2.0: renamed programname from combine to ccdmath
 *				new filestruct, although code still same
 *      18-dec-88	V2.1: Gipsy's fie() functions: order and name
 *                      of keywords adapted to fie()-stuff
 *	11-jan-89	V2.2: be able to create maps from scratch
 *                            needs the size=nx,ny keyword in this case
 *	22-jan-89	V2.3a: map can also be full 3D
 *	19-jun-89	V2.3b: fie -> fien for F2C interface
 *	 7-jul-89	V2.3c: nemoinpi instead of nemoinp, back to fie again
 *	13-nov-90       V2.4 get_nan is now in library in NEMO V2.x
 *			     also forget to return value in fie_remap
 *	 7-mar-92	V2.5 modern style - gcc 2.0 happy		PJT
 *				fie() is now in real, not float
 *	24-feb-98	V2.6 added seed=				PJT
 *	 8-sep-01	a    init_xrandom
 *      24-jul-02       2.7  allow nz=0, so it does not make a 3rd axis  PJT
 *       1-mar-03       3.0  new keywords to set/change WCS              PJT 
 *      19-jun-03       3.1  allow %w, %r for 2d/3d radii; allow crpix offset     PJT
 *      26-aug-04       3.2  fix bad error in setting crpix for cube generation   PJT
 *      10-may-05       3.2a use the wcs routines that have moved to wcsio.c      PJT
 *
 *       because of the float/real conversions and
 *       to eliminate excessive memory usage, operations 'fie' are
 *       done on a column by column basis.
 *                      
 */

#include <stdinc.h>
#include <ctype.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <strlib.h>
#include <image.h>

string defv[] = {
  "in=\n           Input file(s), separated by comma's (optional)",
  "out=???\n       Output file",
  "fie=???\n       Expression %1,%2,.. for input maps; %x,%y,%z, %w,%r for new map/cube",
  "size=10,10,1\n  2- or 3D dimensions of map/cube",
  "crpix=\n        Override/Set crpix (1,1,1)",
  "crval=\n        Override/Set crval (0,0,0)",
  "cdelt=\n        Override/Set cdelt (1,1,1)",
  "seed=0\n        Random seed",
  "VERSION=3.2b\n  24-nov-2019 PJT",
  NULL,
};

string usage = "image arithmetic and wcs";

string cvsid="$Id$";

#ifndef HUGE
# define HUGE 1.0e20
#endif

#define MAXIMAGE 20

imageptr iptr[MAXIMAGE];	/* pointers to (input) images */
int      nimage;                /* actual number of input images */
bool     mapgen = FALSE;	/* no input files: create from scratch ? */

#define MAXNAX 3

double crval[MAXNAX], crpix[MAXNAX], cdelt[MAXNAX];
int nwcs = 0;



local int set_axis(string var, int n, double *xvar, double defvar);
local int fie_remap(char *fie, bool map_create);
local void do_create(int nx, int ny, int nz);
local void do_combine(void);

extern  int debug_level;		/* see initparam() */
extern    void    dmpfien();
extern    int     inifien();
extern    void    dofien();

extern string *burststring(string,string);


void nemo_main ()
{
    string *fnames, fie;
    stream  instr[MAXIMAGE];            /* input files */
    stream  outstr;                     /* output file */
    int     size[3], nx, ny, nz;        /* size of scratch map */
    int     noper;                      /* number of images needed from oper */

    init_xrandom(getparam("seed"));
    nwcs += set_axis(getparam("crpix"),MAXNAX,crpix,1.0);
    nwcs += set_axis(getparam("crval"),MAXNAX,crval,0.0);
    nwcs += set_axis(getparam("cdelt"),MAXNAX,cdelt,1.0);
    fie = getparam("in");
    if (fie==NULL || *fie==0)
        mapgen = TRUE;      /* no input files: create maps from scratch */
    else
        fnames = burststring(getparam("in"), ", ");  /* input file names */
    if(mapgen) {            /* map from scratch */
        dprintf(1,"Generating a map from scratch\n");
        switch (nemoinpi(getparam("size"),size,3)) {
            case 1:			/*  nx[,nx,1] */
                size[1] = size[0];
		size[2] = 1;
                break;
            case 2:			/*  nx,ny[,1] */
		size[2] = 1;
                break;
	    case 3:			/*  nx,ny,nz  */
                break;
            case 0:			/*  [10,10,1] */
                dprintf(0,"Cannot have no size, default 10 assumed\n");
                size[0] = size[1] = 10;
                size[2] = 1;
                break;
            default:			/* --- some error --- */
                error("Syntax error in size keyword\n");
        }
        nx = size[0];
        ny = size[1];
        nz = size[2];
    }
    fie = scopy(getparam("fie"));
    if (fie_remap(fie,mapgen) < 0)      /* remap %x,%y,%z to %1,%2,%3 */
        error("syntax error in fie = %s",fie);
    noper = inifien(fie);                   /* highest parameter # needed */
    if (debug_level >= 5)  dmpfien();        /* debug output from inifie */
    if (noper < 0)
         error ("Error in parsing fie expression %s\n inifie returned %d",
                 getparam("fie"),noper);

    outstr = stropen (getparam("out"),"w");  /* open output file first ... */

    nimage = 0;         /* count number of images/files in in= keyword */
    while (!mapgen){         /* .. then open input files one by one */
        if (nimage==MAXIMAGE)
            error("Too many input files, maximum is %d\n",MAXIMAGE);
        if (fnames[nimage]==NULL)
            break;                      /* done with file names */
        instr[nimage] = stropen(fnames[nimage],"r");    /* open file */
        iptr[nimage] = NULL;        /* make sure to init it right */
        read_image (instr[nimage], &iptr[nimage]);
        dprintf (2,"Image %d read in\n",nimage);
        if (nimage) {                   /* check size consistency */
            if (Nx(iptr[nimage]) != Nx(iptr[nimage-1]))
                error ("Input map %d does have different Nx\n",nimage);
            if (Ny(iptr[nimage]) != Ny(iptr[nimage-1]))
                error ("Input map %d does have different Ny\n",nimage);
            if (Nz(iptr[nimage]) != Nz(iptr[nimage-1]))
                error ("Input map %d does have different Nz\n",nimage);
        }
        strclose(instr[nimage]);        /* close input file */
        nimage++;
    }
    if (!mapgen) {
        dprintf (1,"%d image files read\n",nimage);
        if (nimage < noper)
            error("Not enough input maps (%d) for %d operations\n",
                nimage, noper);
    }

    if (mapgen)
        do_create(nx,ny,nz);
    else
        do_combine();

    write_image (outstr,iptr[0]);         /* write image to file */
    strclose(outstr);
}


local int set_axis(string var, int n, double *xvar, double defvar)
{
  int i, nret;
  if (var == 0 || *var == 0) {
    for (i=0; i<n; i++) 
      xvar[i] = defvar;
    return 0;
  } else {
    nret = nemoinpd(var,xvar,n);
    if (nret < 0) error("Parsing error %d in %s",nret,var);
    for (i=nret; i<n; i++)
      xvar[i] = defvar;
    return 1;
  }
  return 999;
}



/* 
 *  Checks if string 'fie' has nothing but numbers
 *  after a %, i.e. %1, %2... when map_create==false
 *  In case map_create==true % can only be followed by
 *  an 'x', 'y' or 'z' - in this case x -> 1, y -> 2, z -> 3
 *  'w' is sqrt(x*x+y*y),   'r' is sqrt(x*x+y*y+z*z).
 *  Returns 0 on success, -1 on some failure
 */

local int fie_remap(char *fie, bool map_create)
{
    for(;*fie;fie++) {
        fie = strpbrk(fie,"%$");                /* look for a $ or % */
        if (fie==NULL)
            return(0);                          /* done if End of String */
        fie++;
        if (*fie=='x' && map_create)            /* x */
            *fie = '1';
        else if (*fie=='y' && map_create)       /* y */
            *fie = '2';
        else if (*fie=='z' && map_create)       /* z */
            *fie = '3';
        else if (*fie=='w' && map_create)       /* w (2d) */
            *fie = '4';
        else if (*fie=='r' && map_create)       /* r (3D) */
            *fie = '5';
        else if (isdigit(*fie) && !map_create)
            ;                                   /* ignore */
        else {
            dprintf(1,"Warning: syntax error starting at fie=%s\n",--fie);
            return(-1);
        }
    }
    return 0;
}

/*
 *  create new map from scratch, using %x and %y as position parameters 
 *		0..nx-1 and 0..ny-1
 */
local void do_create(int nx, int ny,int nz)
{
    double m_min, m_max, total;
    real   fin[5], fout;
    int    ix, iy, iz;
    int    badvalues;
    
    m_min = HUGE; m_max = -HUGE;
    total = 0.0;		/* count total intensity in new map */
    badvalues = 0;		/* count number of bad operations */

    if (nz > 0) {
      if (!create_cube (&iptr[0], nx, ny, nz))	/* create default empty image */
        error("Could not create 3D image from scratch");
      wcs_f2i(3,crpix,crval,cdelt,iptr[0]);

      for (iz=0; iz<nz; iz++) {
        fin[2] = iz-crpix[2]+1;   /* crpix is 1 for first pixel (FITS convention) */
        for (iy=0; iy<ny; iy++) {
	  fin[1] = iy-crpix[1]+1;
	  for (ix=0; ix<nx; ix++) {
	    fin[0] = ix-crpix[0]+1;
	    fin[3] = sqrt(sqr(fin[0])+sqr(fin[1]));              /* w */
	    fin[4] = sqrt(sqr(fin[0])+sqr(fin[1])+sqr(fin[2]));  /* r */
	    dofien (fin, 1, &fout, 0.0);   /* do the work --- see: fie.3 */
	    CubeValue(iptr[0],ix,iy,iz) = fout;
	    m_min = MIN(m_min,fout);         /* and check for new minmax */
	    m_max = MAX(m_max,fout);
	    total += fout;                   /* add up totals */
	  }
        }
      }
    } else {
      if (!create_image (&iptr[0], nx, ny))	
        error("Could not create 2D image from scratch");
      wcs_f2i(2,crpix,crval,cdelt,iptr[0]);

      for (iy=0; iy<ny; iy++) {
	fin[1] = iy;
	for (ix=0; ix<nx; ix++) {
	  fin[0] = ix;
	  dofien (fin, 1, &fout, 0.0);   /* do the work --- see: fie.3 */
	  MapValue(iptr[0],ix,iy) = fout;
	  m_min = MIN(m_min,fout);         /* and check for new minmax */
	  m_max = MAX(m_max,fout);
	  total += fout;                   /* add up totals */
	}
      }
    } 
    
    MapMin(iptr[0]) = m_min;
    MapMax(iptr[0]) = m_max;

    dprintf(1,"New min and max in map are: %f %f\n",m_min,m_max);
    dprintf(1,"New total brightness/mass is %f\n",
			total*Dx(iptr[0])*Dy(iptr[0]));
    if (badvalues)
    	warning ("There were %d bad operations in dofie",badvalues);
}

/* 
 *  combine input maps into an output map  -- still 2D only
 */
local void do_combine()
{
    double m_min, m_max, total;
    real  *fin, *fout;
    int    k, ix, iy, iz, nx, ny, nz, offset;
    int    badvalues;
    
    m_min = HUGE; m_max = -HUGE;
    total = 0.0;		/* count total intensity in new map */
    badvalues = 0;		/* count number of bad operations */

    nx = Nx(iptr[0]);
    ny = Ny(iptr[0]);
    nz = Nz(iptr[0]);
    if (nwcs) {
      if (nwcs==3) 
	wcs_f2i(2,crpix,crval,cdelt,iptr[0]);
      else
	warning("Not enough WCS information given (%d/3 keywords) to replace it",nwcs);
    }

    fin = (real *) allocate(nimage*ny*sizeof(real)); 
    fout = (real *) allocate(ny*sizeof(real));
        
    for (iz=0; iz<nz; iz++)
    for (ix=0; ix<nx; ix++) {
        for (k=0; k<nimage; k++) {       /* prepare input column buffer */
            offset = ny*k;
            for (iy=0; iy<ny; iy++)
                fin[iy+offset] = CubeValue(iptr[k],ix,iy,iz);
        }
        dofien (fin, ny, fout, 0.0); /* do the work --- see: nemofie.3 */
        for (iy=0; iy<ny; iy++) {             /* write buffer back to map-0 */
            CubeValue(iptr[0],ix,iy,iz) = (real) fout[iy];
            m_min = MIN(m_min,fout[iy]);         /* and check for new minmax */
            m_max = MAX(m_max,fout[iy]);
            total += fout[iy];
        }
    }
    free(fin);
    free(fout);    

    MapMin(iptr[0]) = m_min;
    MapMax(iptr[0]) = m_max;

    dprintf(1,"New min and max in map are: %f %f\n",m_min,m_max);
    dprintf(1,"New total brightness/mass is %f\n",
			total*Dx(iptr[0])*Dy(iptr[0]));
    if (badvalues)
    	warning("There were %d bad operations in dofie",badvalues);
    
}

