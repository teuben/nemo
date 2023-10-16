/* 
 * CCDMOMS: take a moment,or combine, images/cubes from different files
 *          like ccdmom with axis=4
 *
 *	quick and dirty:  15-dec-2020		pjt
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <moment.h>
#include <extstring.h>
#include <ctype.h>

string defv[] = {
  "in=???\n       Input image files",
  "out=???\n      Output image file",
  "mom=1\n	  Moment to take [0=sum,1=mean,2=disp,-1=median]",
  "weight=\n      Scalar weight per image [1 for all]",
  "axis=\n        WCS value per image, if needed",
  "win=\n         Input image weights, one per input (not implemented)",
  "sigma=f\n      Are the weights still SIGMA (t), or straight weights (f)",
  "VERSION=0.3\n  18-dec-2020 PJT",
  NULL,
};

string usage = "(weighted) moment along a set of images";
string cvsid="$Id$";


#ifndef HUGE
# define HUGE 1.0e35
#endif

#define MAXIMAGE 20

imageptr iptr[MAXIMAGE];	/* pointers to (input) images */
imageptr wptr[MAXIMAGE];	/* pointers to (input) weight images */
real iwt[MAXIMAGE];             /* scalar weight per image */
real iax[MAXIMAGE];             /* WCS axis value per image */
int      nimage;                /* actual number of input images */
int      mom;                   /* special moment to take */
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
    bool    Qwin = hasvalue("win");

    mom = getiparam("mom");
    Qsigma = getbparam("sigma");    

    fnames = burststring(getparam("in"), ", ");  /* input file names */
    nimage = xstrlen(fnames, sizeof(string)) - 1;
    if (nimage > MAXIMAGE) error("Too many images %d > %d", nimage, MAXIMAGE);
    dprintf(0,"Using %d images\n",nimage);

    if (Qwin) {
	wnames = burststring(getparam("win"), ", ");	
	n = xstrlen(wnames, sizeof(string)) - 1;
	if (n != nimage) error("Not the correct number of win= images (%d != %d)",n,nimage);
	if (n>0)
	  warning("The win= keyword is currently not implemented");
    }

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
    

    n = nemoinpr(getparam("axis"), iax, nimage);
    if (n<0)
      error("Parsing %s", getparam("axis"));
    else if (n==0)
      for (l=0; l<nimage; l++)  iax[l] = l;
    else if (n!=nimage)
      error("Cannot handle %d values for axis=",n);
    dprintf(0,"axis=%g %g ...\n", iax[0], iax[1]);
    if (n>0)
      warning("The axis= keyword is currently not implemented");


    outstr = stropen (getparam("out"),"w");  /* open output file first ... */

    for (l=0; l<nimage; l++) {
        instr   = stropen(fnames[l],"r");    /* open file */
        iptr[l] = NULL;        /* make sure to init it right */
        wptr[l] = NULL;        /* make sure to init it right */	
        read_image (instr, &iptr[l]);
        dprintf (2,"Image %d read in\n",l);
        if (l>0) {                   /* check size consistency */
            if (Nx(iptr[l]) != Nx(iptr[0]))
                error ("Input map %d does have different Nx\n",l);
            if (Ny(iptr[l]) != Ny(iptr[0]))
                error ("Input map %d does have different Ny\n",l);
            if (Nz(iptr[l]) != Nz(iptr[0]))
                error ("Input map %d does have different Nz\n",l);
        }
        strclose(instr);        /* close input file */
	if (Qwin) {
	  dprintf(1,"win=%s weight map\n",wnames[l]);
	  instr = stropen(wnames[l],"r"); 
	  read_image (instr, &wptr[l]);
	  if (Nx(wptr[l]) != Nx(iptr[0]))
	    error ("Input map %d does have different Nx\n",l);
	  if (Ny(wptr[l]) != Ny(iptr[0]))
	    error ("Input map %d does have different Ny\n",l);
	  if (Nz(wptr[l]) != Nz(iptr[0]))
	    error ("Input map %d does have different Nz\n",l);
	  strclose(instr);
	}
	
    }
    do_combine();

    write_image (outstr,iptr[0]);         /* write image to file */
    strclose(outstr);
}


/* 
 *  combine input maps into an output map  --
 *
 *  the axis= has not been implemented here
 */
local void do_combine()
{
    double m_min, m_max;
    real  *fin, *win, fout;
    int    k, l, ix, iy, iz, nx, ny, nz, offset;
    int    badvalues;
    Moment m;
    
    m_min = HUGE; m_max = -HUGE;
    badvalues = 0;		/* count number of bad operations */

    nx = Nx(iptr[0]);
    ny = Ny(iptr[0]);
    nz = Nz(iptr[0]);

    dprintf(0,"mom=%d\n",mom);
    ini_moment(&m, 2, nimage);

    fin = (real *) allocate(nimage*sizeof(real)); 
        
    for (iz=0; iz<nz; iz++)
      for (iy=0; iy<ny; iy++)
	for (ix=0; ix<nx; ix++) {
	  reset_moment(&m);
	  for (l=0; l<nimage; l++) {
	    fin[l] = CubeValue(iptr[l],ix,iy,iz);
	    if (wptr[l]) {
	      iwt[l] = CubeValue(wptr[l],ix,iy,iz);
	      if (Qsigma) iwt[l] = 1.0/(iwt[l]*iwt[l]);
	    }
	    accum_moment(&m, fin[l], iwt[l]);
	  } // l
	  if (mom==0)
	    fout = sum_moment(&m);
	  else if (mom==1)
	    fout = mean_moment(&m);
	  else if (mom==2)
	    fout = sigma_moment(&m);
	  else if (mom==-1)
	    fout = median_moment(&m);
	  else
	    error("Invalid mom=%d",mom);
	    
	  CubeValue(iptr[0],ix,iy,iz) = fout;
	  m_min = MIN(m_min,fout);
	  m_max = MAX(m_max,fout);
        } // ix
    // iz,iy
    
    free(fin);

    MapMin(iptr[0]) = m_min;
    MapMax(iptr[0]) = m_max;

    dprintf(0,"New min and max in map are: %f %f\n",m_min,m_max);
    if (badvalues)
    	warning("There were %d bad operations in dofie",badvalues);
    
}



