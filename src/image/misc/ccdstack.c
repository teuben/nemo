/* 
 * CCDSTACK: stack images, with simple gridding option
 *
 *   21-may-2021:    derived from ccdmoms, but should not need to allocate MAXIMAGE, just use 2
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <extstring.h>
#include <ctype.h>

string defv[] = {
  "in=???\n       Input image files, first image sets the WCS",
  "out=???\n      Output image file",
  "weight=\n      Scalar weight per image [1 for all]",
  "sigma=f\n      Are the weights still SIGMA (t), or straight weights (f)",
  "bad=0\n        Bad value to ignore",
  "VERSION=0.1\n  21-may-2021 PJT",
  NULL,
};

string usage = "stack images, with simple gridding option if WCS differs";
string cvsid="$Id$";


#ifndef HUGE
# define HUGE 1.0e35
#endif

#define MAXIMAGE 100

imageptr iptr[MAXIMAGE];	/* pointers to (input) images */
real iwt[MAXIMAGE];             /* scalar weight per image */
real     badval;
int      nimage;                /* actual number of input images */
bool     Qsigma;

local void do_combine(void);



void nemo_main ()
{
    string *fnames, *wnames;
    stream  instr;                      /* input files */
    stream  outstr;                     /* output file */
    int     size[3], nx, ny, nz;        /* size of scratch map */
    int     noper;                      /* number of images needed from oper */
    int     i, j, k, l, n;

    Qsigma = getbparam("sigma");
    badval = getrparam("bad");

    fnames = burststring(getparam("in"), ", ");  /* input file names */
    nimage = xstrlen(fnames, sizeof(string)) - 1;
    if (nimage > MAXIMAGE) error("Too many images %d > %d", nimage, MAXIMAGE);
    dprintf(0,"Using %d images\n",nimage);

    n = nemoinpr(getparam("weight"), iwt, nimage);
    if (n<0)
      error("Parsing %s", getparam("weight"));
    else if (n==0)
      for (l=0; l<nimage; l++)  iwt[l] = 1.0;
    else if (n!=nimage)
      error("Cannot handle %d values for weight=",n);
    if (Qsigma)
      for (l=0; l<nimage; l++)  iwt[l] = 1/(iwt[l]*iwt[l]);
    
    dprintf(0,"weight=%g %g ...\n", iwt[0], iwt[1]);
    

    outstr = stropen (getparam("out"),"w");  /* open output file first ... */

    for (l=0; l<nimage; l++) {
        instr   = stropen(fnames[l],"r");    /* open file */
        iptr[l] = NULL;        /* make sure to init it right */
        read_image (instr, &iptr[l]);
        dprintf (2,"Image %d read in, minmax %g %g\n",l,MapMin(iptr[l]),MapMax(iptr[l]));
        strclose(instr);        /* close input file */
    }
    do_combine();

    write_image(outstr,iptr[0]);         /* write image to file */
    strclose(outstr);
}


/* 
 *  combine input maps into an output map  --
 *
 */
local void do_combine()
{
    double m_min, m_max;
    real  *fin, *win, fout;
    int    k, l, ix, iy, iz, nx, ny, nz, offset;
    int    badvalues;
    imageptr wptr;
    
    m_min = HUGE; m_max = -HUGE;
    badvalues = 0;		/* count number of bad operations */

    nx = Nx(iptr[0]);
    ny = Ny(iptr[0]);
    nz = Nz(iptr[0]);

    create_cube(&wptr, nx, ny, nx);
        
    for (l=1; l<nimage; l++) {
      dprintf(0,"Adding image %d\n",l);
      for (iz=0; iz<Nz(iptr[l]); iz++)
	for (iy=0; iy<Ny(iptr[l]); iy++)
	  for (ix=0; ix<Nx(iptr[l]); ix++) {
	    if (l==1) {  // initialize first time around
	      if (CubeValue(iptr[0],ix,iy,iz) == badval)
		CubeValue(wptr,ix,iy,iz) = 0.0;
	      else
		CubeValue(wptr,ix,iy,iz) = 1.0;		
	    } // accumulate
	    if (CubeValue(iptr[l],ix,iy,iz) != badval) {
	      CubeValue(iptr[0],ix,iy,iz) += CubeValue(iptr[l],ix,iy,iz);
	      CubeValue(wptr,ix,iy,iz) += 1.0;
	    } 
	  }
      
    }
    for (iz=0; iz<nz; iz++)
      for (iy=0; iy<ny; iy++)
	for (ix=0; ix<nx; ix++) {
	  if (CubeValue(wptr,ix,iy,iz) > 0)
	    CubeValue(iptr[0],ix,iy,iz) /= CubeValue(wptr,ix,iy,iz);
	  // CubeValue(iptr[0],ix,iy,iz) = CubeValue(wptr,ix,iy,iz);	  
	  if (CubeValue(iptr[0],ix,iy,iz) > m_max) m_max = CubeValue(iptr[0],ix,iy,iz);
	  if (CubeValue(iptr[0],ix,iy,iz) < m_min) m_min = CubeValue(iptr[0],ix,iy,iz);
	}
    
	
    MapMin(iptr[0]) = m_min;
    MapMax(iptr[0]) = m_max;

    dprintf(0,"New min and max in map are: %f %f\n",m_min,m_max);
    if (badvalues)
    	warning("There were %d bad operations in dofie",badvalues);
    
}



