/* 
 * CCDSTAT
 *
 *	 4-oct-88  V1.0 created
 *	 5-nov-93  V1.1 nemo V2, and using moment.h
 *      16-sep-95  V1.2 added min/max
 *       9-may-03  V1.3 added bad=, and added #points,
 *       5-jun-03  V1.4 added win=, a weight map
 *      14-nov-04   1.5 provide nppb= correction for chi2 computation   PJT
 *       5-jan-05   1.6 added a total
 *      30-jan-05   1.7 added an optional median
 *      24-may-06   1.8 add mmcount=
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
    "median=f\n     Optional display of the median value",
    "mmcount=f\n    Count occurances of min and max",
    "sort=qsort\n   Sorting routine (not activated yet)",
    "VERSION=1.8\n  24-may-06 PJT",
    NULL,
};

string usage="basic statistics of an image, optional chi2 calculation";

string cvsid="$Id$";

string	infile;	        		/* file names */
stream  instr;				/* file streams */

imageptr iptr=NULL;			/* will be allocated dynamically */
imageptr wptr=NULL;                     /* optional weight map */


int    nx,ny,nz,nsize;			/* actual size of map */
double xmin,ymin,zmin,dx,dy,dz;
double size;				/* size of frame (square) */
double cell;				/* cell or pixel size (square) */

real get_median(int n, real *x);


nemo_main()
{
    int  i, j, k;
    real x, xmin, xmax, mean, sigma, skew, kurt, median, bad, w, *data;
    Moment m;
    bool Qmin, Qmax, Qbad, Qw, Qmedian, Qmmcount = getbparam("mmcount");
    real nu, nppb = getdparam("nppb");
    int npar = getiparam("npar");
    int ngood = 0;
    int min_count, max_count;

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
    Qmedian = getbparam("median");
    if (Qmedian)
      data = (real *) allocate(nx*ny*nz*sizeof(real));
    
    ini_moment(&m,4,0);
    for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
        for (k=0; k<nz; k++) {
            x =  CubeValue(iptr,i,j,k);
            if (Qmin && x<xmin) continue;
            if (Qmax && x>xmax) continue;
            if (Qbad && x==bad) continue;
	    w = Qw ? CubeValue(wptr,i,j,k) : 1.0;
            accum_moment(&m,x,w);
	    if (Qmedian) data[ngood++] = x;
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
      printf ("Sum                  : %f\n",show_moment(&m,1));
      if (Qmedian)
	printf ("Median               : %f\n",get_median(ngood,data));

      if (Qmmcount) {
	min_count = max_count = 0;
	xmin = min_moment(&m);
	xmax = max_moment(&m);
	for (i=0; i<nx; i++) {
	  for (j=0; j<ny; j++) {
	    for (k=0; k<nz; k++) {
	      x =  CubeValue(iptr,i,j,k);
	      if (x==xmin) min_count++;
	      if (x==xmax) max_count++;
	    }
	  }
	} /* i */
	printf("Min_Max_count        : %d %d\n",min_count,max_count);
      }
      printf ("%d/%d out-of-range points discarded\n",nsize-n_moment(&m), nsize);
    }
}

int compar_real(real *a, real *b)
{
  return *a < *b ? -1 : *a > *b ? 1 : 0;
}

real get_median(int n, real *x)
{
  qsort(x,n,sizeof(real),compar_real);
  if (n % 2)
    return  x[(n-1)/2];
  else
    return 0.5 * (x[n/2] + x[n/2-1]);
}
