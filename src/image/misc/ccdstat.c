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
 *       6-nov-06   1.9 output the total mass/luminosity as well as sum of densities
 *      13-feb-13   2.0 ignore axis if N=1 
 *      17-apr-13   3.0 allow each Z plane to be treated seperately
 *      23-apr-13       half= option, use compute_robust_moment()
 *       6-jun-13   3.2 added torben, cleanup a bit
 *      28-feb-14   3.3 added tab=
 *      25-aug-14   3.5 added maxpos=, fixed duplicate header, so format='ascii.commented_header'
 *                      works in astropy
 *    26-dec-2019   3.6 (not finished yet) enable some openmp sections of code
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
    "npar=0\n       Number of fitting parameters assumed for chi2 calc (full data only)",
    "nppb=1\n       Optional correction 'number of points per beam' for chi2 calc",
    "median=f\n     Optional display of the median value",
    "torben=f\n     Use torben method for median instead",
    "robust=f\n     Compute robust median",
    "mmcount=f\n    Count occurances of min and max",
    "maxpos=f\n     Add location of where the max occured",
    "half=f\n       Only use half (negative) values and symmetrize them",
    "maxmom=4\n     Control how many moments are computed",
    "ignore=t\n     (for summing) Ignore cell width when N=1 (assumed infinity)",
    "sort=qsort\n   Sorting routine (not activated yet)",
    "planes=-1\n    -1: whole cube in one      0=all planes   start:end:step = selected planes",
    "tab=\n         If given, print out data values",
    "VERSION=3.5b\n 26-dec-2019 PJT",
    NULL,
};

string usage="basic statistics of an image, optional chi2 calculation";

string cvsid="$Id$";

string	infile;	        		/* file names */
stream  instr;				/* file streams */
stream  tabstr = NULL;                  /* output table? */

imageptr iptr=NULL;			/* will be allocated dynamically */
imageptr wptr=NULL;                     /* optional weight map */



int    nx,ny,nz,nsize;			/* actual size of map */
double xmin,ymin,zmin,dx,dy,dz;
double size;				/* size of frame (square) */
double cell;				/* cell or pixel size (square) */
int *planes = NULL;                     /* selected planes , if used */

real get_median(int n, real *x);

extern real median(int,real*);
extern real median_q1(int,real*);
extern real median_q3(int,real*);


nemo_main()
{
    int  i, j, k, ki;
    real x, y, z, xmin, xmax, mean, sigma, skew, kurt,  bad, w, *data;
    real dmin, dmax;
    real sum, sov, q1, q2, q3;
    Moment m;
    bool Qmin, Qmax, Qbad, Qw, Qmedian, Qrobust, Qtorben, Qmmcount = getbparam("mmcount");
    bool Qx, Qy, Qz, Qone, Qall, Qign = getbparam("ignore");
    bool Qhalf = getbparam("half");
    bool Qmaxpos = getbparam("maxpos");
    real nu, nppb = getdparam("nppb");
    int npar = getiparam("npar");
    int ngood;
    int ndat = 0;
    int nplanes;
    int min_count, max_count;
    int maxmom = getiparam("maxmom");
    int maxpos[2];
    char slabel[32];

    instr = stropen (getparam("in"), "r");
    read_image (instr,&iptr);
    strclose(instr);
    nx = Nx(iptr);	
    ny = Ny(iptr);
    nz = Nz(iptr);
    dprintf(1,"# data order debug:  %f %f\n",Frame(iptr)[0], Frame(iptr)[1]);
    if (hasvalue("tab")) tabstr = stropen(getparam("tab"),"w");

    planes = (int *) allocate((nz+1)*sizeof(int));
    nplanes = nemoinpi(getparam("planes"),planes,nz+1);
    Qall = (planes[0]==-1);
    if (planes[0]==0) {
      Qone = FALSE;
      nplanes = nz;
      for (k=0; k<nz; k++)
	planes[k] = k;
    } else if (!Qall) {
      Qone = (!Qall && nplanes==1);
      for (k=0; k<nplanes; k++) {
	if (planes[k]<1 || planes[k]>nz) error("%d is an illegal plane [1..%d]",planes[k],nz);
	planes[k]--;
      }
    }
 
    if (hasvalue("win")) {
      instr = stropen (getparam("win"), "r");
      read_image (instr,&wptr);
      strclose(instr);
      if (Nx(iptr) != Nx(wptr)) error("X sizes of in=/win= don't match");
      if (Ny(iptr) != Ny(wptr)) error("Y sizes of in=/win= don't match");
      if (Nz(iptr) != Nz(wptr)) error("Z sizes of in=/win= don't match");
      Qw = TRUE;
    } else
      Qw = FALSE;

    Qmin = hasvalue("min");
    if (Qmin) xmin = getdparam("min");
    Qmax = hasvalue("max");
    if (Qmax) xmax = getdparam("max");
    Qbad = hasvalue("bad");
    if (Qbad) bad = getdparam("bad");
    Qmedian = getbparam("median");
    Qrobust = getbparam("robust");
    Qtorben = getbparam("torben");
    if (Qtorben) Qmedian = TRUE;
    if (Qmedian || Qrobust || Qtorben) {
      ndat = nx*ny*nz;
      dprintf(1,"Need spare array size %d\n",ndat);
      data = (real *) allocate(ndat*sizeof(real));
    }

    sov = 1.0;       /* volume of a pixel/voxel */
    Qx = Qign && Nx(iptr)==1;
    Qy = Qign && Ny(iptr)==1;
    Qz = Qign && Nz(iptr)==1;
    sov *= Qx ? 1.0 : Dx(iptr);
    sov *= Qy ? 1.0 : Dy(iptr);
    sov *= Qz ? 1.0 : Dz(iptr);
    strcpy(slabel,"*");
    if (!Qx) strcat(slabel,"Dx*");
    if (!Qy) strcat(slabel,"Dy*");
    if (!Qz) strcat(slabel,"Dz*");
    /* make it always 10 in lenght */
    for (i=strlen(slabel); i<10; i++)
      strcat(slabel," ");
    
    if (maxmom < 0) {
      warning("No work done, maxmom<0");
      stop(0);
    }

    if (Qall) {                 /* treat cube as one data block */

      ini_moment(&m,maxmom,ndat);
      ngood = 0;
#if 0
      // not working yet
      #pragma omp parallel \
	shared(ngood, nx,ny,nz, iptr,wptr,m, xmin,xmax,bad,Qhalf,Qmin,Qmax,Qbad,Qw,data,tabstr,Qmedian) \
        private(i,j,k,x,w)
      #pragma omp for
#endif      
      for (k=0; k<nz; k++) {
	for (j=0; j<ny; j++) {
	  for (i=0; i<nx; i++) {
            x =  CubeValue(iptr,i,j,k);    // iptr->cube[k,j,k]
	    if (Qhalf && x>=0.0) continue;
            if (Qmin  && x<xmin) continue;
            if (Qmax  && x>xmax) continue;
            if (Qbad  && x==bad) continue;
	    w = Qw ? CubeValue(wptr,i,j,k) : 1.0;
            accum_moment(&m,x,w);
	    if (Qhalf && x<0) accum_moment(&m,-x,w);
	    if (Qmedian) data[ngood++] = x;
	    if (tabstr) fprintf(tabstr,"%g\n",x);
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
	sum = mean = sigma = skew = kurt = 0;
	if (maxmom > 0) {
	  mean = mean_moment(&m);
	  sum = show_moment(&m,1);
	}
	if (maxmom > 1)
	  sigma = sigma_moment(&m);
	if (maxmom > 2)
	  skew = skewness_moment(&m);
	if (maxmom > 3)
	  kurt = kurtosis_moment(&m);
	
	printf ("Number of points       : %d\n",n_moment(&m));
	printf ("Min and Max            : %f %f\n",min_moment(&m), max_moment(&m));
	printf ("Mean and dispersion    : %f %f\n",mean,sigma);
	printf ("Skewness and kurtosis  : %f %f\n",skew,kurt);
	printf ("Sum and Sum%s  : %f %f\n",slabel,sum,sum*sov);   /* align trick */
	if (Qmedian) {
	  if (Qtorben) {
	    printf ("Median Torben         : %f (%d)\n",median_torben(ngood,data,min_moment(&m),max_moment(&m)),ngood);
	  } else {
	    printf ("Median                : %f\n",get_median(ngood,data));
	    q2 = median(ngood,data);
	    q1 = median_q1(ngood,data);
	    q3 = median_q3(ngood,data);
	    printf ("Q1,Q2,Q3              : %f %f %f\n",q1,q2,q3);
	  }
#if 1
	  if (ndat>0)
	    printf ("MedianL               : %f\n",median_moment(&m));
#endif
	}
	if (Qrobust) {
	  compute_robust_moment(&m);
	  printf ("Mean Robust           : %f\n",mean_robust_moment(&m));
	  printf ("Sigma Robust          : %f\n",sigma_robust_moment(&m));
	  printf ("Median Robust         : %f\n",median_robust_moment(&m));
	}

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
	  printf("Min_Max_count         : %d %d\n",min_count,max_count);
	}
	printf ("%d/%d out-of-range points discarded\n",nsize-n_moment(&m), nsize);
      }

    } else {             /* treat each plane seperately */

      /* tabular output, one line per (selected) plane */

      printf("# iz z min  max  N  mean sigma skew kurt sum sumsov ");
      if (Qmedian) printf(" med1 med2");
      if (Qrobust) printf(" rN rmean rsigma rmed]");
      if (Qmaxpos) printf(" maxposx maxposy");
      printf("\n");

      ini_moment(&m,maxmom,ndat);
      for (ki=0; ki<nplanes; ki++) {
	reset_moment(&m);
	k = planes[ki];
	z = Zmin(iptr) + k*Dz(iptr);
	ngood = 0;
	for (j=0; j<ny; j++) {
	  for (i=0; i<nx; i++) {
            x =  CubeValue(iptr,i,j,k);
            if (Qmin && x<xmin) continue;
            if (Qmax && x>xmax) continue;
            if (Qbad && x==bad) continue;
	    if (Qmaxpos) {
	      if (i==0 && j==0) { dmax = x; maxpos[0] = 0; maxpos[1] = 0;}
	      else if (x>dmax) {  dmax = x; maxpos[0] = i; maxpos[1] = j;}
	    }
	    w = Qw ? CubeValue(wptr,i,j,k) : 1.0;
            accum_moment(&m,x,w);
	    if (Qmedian) data[ngood++] = x;
	  }
	}

	nsize = nx * ny * nz;
	sum = mean = sigma = skew = kurt = 0;
	if (maxmom > 0) {
	  mean = mean_moment(&m);
	  sum = show_moment(&m,1);
	}
	if (maxmom > 1)
	  sigma = sigma_moment(&m);
	if (maxmom > 2)
	  skew = skewness_moment(&m);
	if (maxmom > 3)
	  kurt = kurtosis_moment(&m);
	if (n_moment(&m) == 0) {
	  printf("# %d no data\n",k+1);
	  continue;
	}
	printf("%d %f  %f %f %d  %f %f %f %f  %f %f",
	       k+1, z, min_moment(&m), max_moment(&m), n_moment(&m),
	       mean,sigma,skew,kurt,sum,sum*sov);
	if (Qmedian) {
	  printf ("   %f",get_median(ngood,data));
	  if (ndat>0) printf (" %f",median_moment(&m));
	}
	if (Qrobust) {
	  compute_robust_moment(&m);
	  printf ("   %d %f %f %f",n_robust_moment(&m), mean_robust_moment(&m),
		  sigma_robust_moment(&m), median_robust_moment(&m));
	}
	if (Qmaxpos) {
	  printf("   %d %d",maxpos[0],maxpos[1]);
	}
#if 0	  
	if (Qmmcount) {
	  min_count = max_count = 0;
	  xmin = min_moment(&m);
	  xmax = max_moment(&m);
	  for (i=0; i<nx; i++) {
	    for (j=0; j<ny; j++) {
	      x =  CubeValue(iptr,i,j,k);
	      if (x==xmin) min_count++;
	      if (x==xmax) max_count++;
	    }
	  } /* i,j */
	  printf(" %d %d",min_count,max_count);
	}
	printf ("%d/%d out-of-range points discarded\n",nsize-n_moment(&m), nsize);
#endif
	printf("\n");
      } /* ki */
    }
}


/* 
 *   local version of a median. The one in moment.c and median.c are sorting
 *   pointers , this one sorts the data 
 */

int compar_real(real *a, real *b)
{
  return *a < *b ? -1 : *a > *b ? 1 : 0;
}

real get_median(int n, real *x)
{
  dprintf(1,"get_median: n=%d\n",n);

  qsort(x,n,sizeof(real),compar_real);
  if (n % 2)
    return  x[(n-1)/2];
  else
    return 0.5 * (x[n/2] + x[n/2-1]);
}
