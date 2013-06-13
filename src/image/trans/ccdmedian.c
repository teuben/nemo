/* 
 * CCDMEDIAN: median filtering of an image
 *
 *	29-jan-01	PJT	q&d for Mark Wolfire
 *      30-jan-01       PJT     experimenting with NR's sort routines
 *      30-jul-04       PJT     optional X,Y range selection
 *       1-nov-04       PJT     experiment to cf. to python kpno_soft/imsubtract.py
 *       2-nov-04       PJT     added nstep= cheat mode for ShiPing Lai
 *      14-jul-11       PJT     0.6 fixed edge problem
 *       7-aug-12       PJT     0.7 optional median method
 *      12-jun-13       PJT     0.8 mean option  (average)
 *                      
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
        "in=???\n       Input image file",
	"out=???\n      Output image file",
	"n=5\n		(odd) size of filter",
	"x=\n           Optional subselection of the X range (min,max)",
	"y=\n           Optional subselection of the Y range (min,max)",
	"nstep=1\n      Cheat mode: replicate each nstep pixels",
	"fraction=0.5\n Fraction of positive image values in subtract mode",
	"mode=median\n  Mode: median, average, subtract",
	"torben=f\n     Median method",
	"VERSION=0.8\n  12-jun-2013 PJT",
	NULL,
};

string usage = "median filter of an image";

string cvsid = "$Id$";



#if 0
#define CVI(x,y,z)  CubeValue(iptr,x,y,z)
#define CVO(x,y,z)  CubeValue(optr,x,y,z)
#else
#define CVI(x,y)    MapValue(iptr,x,y)
#define CVO(x,y)    MapValue(optr,x,y)
#endif

real median(int, real *, real);
real mean(int, real *, real);
real subtract(int, real *, real);
void sort0(int, real *);
void sort1(int, real *);
void sort2(int, real *);
void sort3(int, real *);
void sort4(int, real *);

extern real median_torben(int n, real *x, real xmin, real xmax);

#define sort sort0

/* @todo
 *
 *   miriad 2048 ran map  size=2 -> 3.1"   size=4 -> 12.3"
 *   nemo                      5 -> 1-"         9    60"
 *   this is terribly much slower than miriad,why?
 *   idl is even faster, size=5->1.9"  9->5.3"
 */


void get_range(string axis, int *ia)
{
  int n = nemoinpi(getparam(axis),ia,2);
  if (n==2) return;
  if (n<0) 
    error("%s=%s parsing error: needs two values (min,max)",axis,getparam(axis));
  error("%s=%s needs two values (min,max)",axis,getparam(axis));
}

void nemo_main()
{
    stream  instr, outstr;
    int     nx, ny, nz;
    int     nstep,nstep1;
    int     i,j,k, n, n1, i1, j1, m;
    int     ix[2], iy[2];
    imageptr iptr=NULL, optr;      /* pointer to images */
    real    *vals, fraction;
    string  mode = getparam("mode");
    bool Qmedian = (*mode == 'm');
    bool Qmean = (*mode == 'a');

    nstep = getiparam("nstep");
    if (nstep%2 != 1) error("step size %d needs to be odd",nstep);
    nstep1 = (nstep-1)/2;

    n = getiparam("n");
    if (Qmedian)
      dprintf(1,"Median filter size %d\n",n);
    else if (Qmean) 
      dprintf(1,"Mean filter size %d\n",n);
    else
      dprintf(1,"Subtraction filter size %d\n",n);
    if (n%2 != 1) error("filter size %d needs to be odd",n);
    n1 = (n-1)/2;
    vals = (real *) allocate (sizeof(real) * (n*n + 1));

    instr = stropen(getparam("in"), "r");
    read_image( instr, &iptr);
    nx = Nx(iptr);	
    ny = Ny(iptr);
    nz = Nz(iptr);
    if (nz > 1) error("Cannot do 3D cubes properly; use 2D");

    if (hasvalue("x") && hasvalue("y")) {
      get_range("x",ix);
      get_range("y",iy);
    } else {
      ix[0] = 0;
      ix[1] = nx-1;
      iy[0] = 0;
      iy[1] = ny-1;
    }
    dprintf(1,"Xrange: %d - %d   Yrange: %d - %d\n",ix[0],ix[1],iy[0],iy[1]);
      
    outstr = stropen(getparam("out"), "w");
    create_cube(&optr,nx,ny,nz);
    Dx(optr) = Dx(iptr);
    Dy(optr) = Dy(iptr);
    Dz(optr) = Dz(iptr);
    Xmin(optr) = Xmin(iptr);
    Ymin(optr) = Ymin(iptr);
    Zmin(optr) = Zmin(iptr);

    if (nstep > 1) {
      warning("Cheat mode nstep=%d",nstep);

      for (j=nstep1; j<ny-nstep1; j+=nstep) {
	for (i=nstep1; i<nx-nstep1; i+=nstep) {
	  if (j<n1 || j >= ny-n1 || j < iy[0] || j > iy[1]) {
	    CVO(i,j) = CVI(i,j);
	    continue;
	  }
	  if (i<n1 || i >= nx-n1 || i < ix[0] || i > ix[1]) {
	    CVO(i,j) = CVI(i,j);
	    continue;
	  }
	  m = 0;
	  for (j1=j-n1; j1<=j+n1; j1++)
	    for (i1=i-n1; i1<=i+n1; i1++)
	      vals[m++] = CVI(i1,j1);
	  CVO(i,j) = median(m,vals,fraction);
	  for (j1=j-nstep1; j1<=j+nstep1; j1++)
	    for (i1=i-nstep1; i1<=i+nstep1; i1++)
	      CVO(i1,j1) = CVO(i,j);
	}
      }
    } else {

      for (j=0; j<ny; j++) {
	for (i=0; i<nx; i++) {
	  if (j<n1 || j >= ny-n1 || j < iy[0] || j > iy[1]) {
	    CVO(i,j) = CVI(i,j);
	    continue;
	  }
	  if (i<n1 || i >= nx-n1 || i < ix[0] || i > ix[1]) {
	    CVO(i,j) = CVI(i,j);
	    continue;
	  }
	  m = 0;
	  for (j1=j-n1; j1<=j+n1; j1++)
	    for (i1=i-n1; i1<=i+n1; i1++)
	      vals[m++] = CVI(i1,j1);

	  if (Qmedian)
	    CVO(i,j) = median(m,vals,fraction);
	  else if (Qmean)
	    CVO(i,j) = mean(m,vals,fraction);
	  else
	    CVO(i,j) = subtract(m,vals,fraction);
	}
      }

    }
    write_image(outstr, optr);
}

real median(int n, real *x, real fraction)
{
  sort(n, x);

  /* TODO: do the correct thing for odd and even n 
   * also look at Knuth's O(N) idea and ShiPing's idea of 
   * histogram reduction 
   * Could also consider fraction != 0.5 
   */
  return x[(n-1)/2];
}

real mean(int n, real *x, real fraction)
{
  real sum = 0.0;
  int i;

  sum = 0.0;
  for (i=0; i<n; i++)
    sum += x[i];
  return sum/n;
}

real subtract(int n, real *x, real fraction)
{
  int i,k;
  sort(n, x);
  for (i=0; i<n; i++)
    if (x[i] > 0) break;
  k = (int) ((n-i)*fraction);
  if (k<0) k=0;
  if (k>=n) k=n-1;
  return x[k];
}

/* very simply shell sort ; k&r pp108 */

void sort0(int n, real *x)
{
    int   gap, i, j;
    real  temp;
    for (gap=n/2; gap>0; gap /= 2)
        for (i=gap; i<n; i++)
            for (j=i-gap; j>=0; j -= gap) {
                if (x[j] <= x[j+gap])
                    break;
                temp = x[j];
                x[j] = x[j+gap];
                x[j+gap] = temp;
            }
}


/* 
 * Numerical Recipes quick-sort routine ;; still in NR notation index-1 based
 */
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

#define NSTACK 50
#define M      7

static int *istack = NULL;

void sort1(int n, real *arr)
{
  int i,ir=n,j,k,l=1;
  int jstack=0, *istack;
  real a,temp;

  if (istack == NULL)
    istack = (int *) allocate(sizeof(int) * NSTACK);
  
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	a=arr[j];
	for (i=j-1;i>=1;i--) {
	  if (arr[i] <= a) break;
	  arr[i+1]=arr[i];
	}
	arr[i+1]=a;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1])
	if (arr[l+1] > arr[ir]) {
	  SWAP(arr[l+1],arr[ir])
	    }
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir])
	  }
      if (arr[l+1] > arr[l]) {
	SWAP(arr[l+1],arr[l])
	  }
      i=l+1;
      j=ir;
      a=arr[l];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j]);
      }
      arr[l]=arr[j];
      arr[j]=a;
      jstack += 2;
      if (jstack > NSTACK) error("NSTACK too small in sort.");
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
}


/*
 * NR heap sort  
 */
void sort3(int n, real *ra)
{
  int i,ir,j,l;
  real rra;

  ra--;    /* NR is 1 based, nemo uses 0-based */
  
  if (n < 2) return;
  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra=ra[--l];
    } else {
      rra=ra[ir];
      ra[ir]=ra[1];
      if (--ir == 1) {
	ra[1]=rra;
	break;
      }
    }
    i=l;
    j=l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
	ra[i]=ra[j];
	i=j;
	j <<= 1;
      } else j=ir+1;
    }
    ra[i]=rra;
  }
}

/* 
 * another heap sort ?
 */

void sort4(int n, real *ra)
{
  *ra = 0.0;
}
