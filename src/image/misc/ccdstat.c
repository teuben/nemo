/* 
 * CCDSTAT
 *
 *	 4-oct-88  V1.0 created
 *	 5-nov-93  V1.1 nemo V2, and using moment.h
 *      16-sep-95  V1.2 added min/max
 *       9-may-03  V1.3 added bad=, and added #points,
 *       5-jun-03  V1.4 added win=, a weight map
 *      14-nov-04   1.5 provide nppb= correction for chi2 computation   PJT
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <image.h>
#include <moment.h> 

string defv[] = {
    "in=???\n       Input image filename",
    "min=\n         Minimum overrride",
    "max=\n         Maximum overrride",
    "bad=\n         Use this as masking value to ignore data",
    "win=\n         Optional input map for weights",
    "npar=0\n       Number of fitting parameters assumed for chi2 calc",
    "nppb=1\n       Optional correction 'number of points per beam' for chi2 calc",
    "VERSION=1.5\n  15-nov-04 PJT",
    NULL,
};

string usage="basic statistics of an image, optional chi2 calculation";

string	infile;	        		/* file names */
stream  instr;				/* file streams */

imageptr iptr=NULL;			/* will be allocated dynamically */
imageptr wptr=NULL;                     /* optional weight map */


int    nx,ny,nz,nsize;			/* actual size of map */
double xmin,ymin,zmin,dx,dy,dz;
double size;				/* size of frame (square) */
double cell;				/* cell or pixel size (square) */



nemo_main()
{
    int  i, j, k;
    real x, xmin, xmax, mean, sigma, skew, kurt, bad, w;
    Moment m;
    bool Qmin, Qmax, Qbad, Qw;
    real nu, nppb = getdparam("nppb");
    int npar = getiparam("npar");
    
    instr = stropen (getparam("in"), "r");
    read_image (instr,&iptr);
    strclose(instr);

    if (hasvalue("win")) {
      instr = stropen (getparam("win"), "r");
      read_image (instr,&wptr);
      strclose(instr);
      if (Nx(iptr) != Nx(wptr)) error("X sizes of in/win don't match");
      if (Ny(iptr) != Ny(wptr)) error("X sizes of in/win don't match");
      if (Nz(iptr) != Nz(wptr)) error("X sizes of in/win don't match");
      Qw = TRUE;
    } else
      Qw = FALSE;

    nx = Nx(iptr);	
    ny = Ny(iptr);
    nz = Nz(iptr);
    Qmin = hasvalue("min");
    if (Qmin) xmin = getdparam("min");
    Qmax = hasvalue("max");
    if (Qmax) xmax = getdparam("max");
    Qbad = hasvalue("bad");
    if (Qbad) bad = getdparam("bad");

    ini_moment(&m,4);
    for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
        for (k=0; k<nz; k++) {
            x =  CubeValue(iptr,i,j,k);
            if (Qmin && x<xmin) continue;
            if (Qmax && x>xmax) continue;
            if (Qbad && x==bad) continue;
	    w = Qw ? CubeValue(wptr,i,j,k) : 1.0;
            accum_moment(&m,x,w);
        }
      }
    }
    if (npar > 0) {
      nu = n_moment(&m)/nppb - npar;
      if (nu < 1) error("%g: No degrees of freedom",nu);
      printf("chi2= %g\n", show_moment(&m,2)/nu/nppb);
      printf("df= %g\n", nu);
    } else {
      nsize = nx * ny * nz;
    
      mean = mean_moment(&m);
      sigma = sigma_moment(&m);
      skew = skewness_moment(&m);
      kurt = kurtosis_moment(&m);
      
      printf ("Min=%f  Max=%f\n",min_moment(&m), max_moment(&m));
      printf ("Number of points     : %d\n",n_moment(&m));
      printf ("Mean and dispersion  : %f %f\n",mean,sigma);
      printf ("Skewness and kurtosis: %f %f\n",skew,kurt);
      printf ("%d/%d out-of-range points discarded\n",nsize-n_moment(&m), nsize);
    }
}
