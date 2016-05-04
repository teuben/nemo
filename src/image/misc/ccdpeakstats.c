/* 
 * CCDPEAKSTATS
 *
 *       1-nov-2015  V0.1   Cloned off ccdstat to play with segments and peaks, the fast way
 *       3-may-2015  V0.2   print FWHM, not sigma
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <image.h>
#include <moment.h> 

string defv[] = {
    "in=???\n       Input image filename",
    "sigma=???\n    Noise level",
    "nsigma=5\n     Number of sigma over which segments are detected",
    "minchan=5\n    Min number of channels over cutoff to count as a segment",
    "maxgap=2\n     Max number of channels below cutoff allowed in a segments",
    "x=\n           Pick an X pixel (0....)",
    "y=\n           Pick an Y pixel (0....)",
    "z=\n           Pick a segment (two numbers) if to override sigma/nsigma",
    "tab=\n         If given, print out data values",
    "VERSION=0.3\n  3-may-2016 PJT",
    NULL,
};

string usage="peak statistics of an image";

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



nemo_main()
{
    int  i, j, k, nseg, ndet, *w, *s0, *s1;
    int x,y;
    bool Qpix;
    real sigma, nsigma, cutoff;
    real *sp;
    int minchan = getiparam("minchan");
    int maxgap  = getiparam("maxgap");
    int ns, seg[2];

    instr = stropen (getparam("in"), "r");
    read_image (instr,&iptr);
    strclose(instr);
    nx = Nx(iptr);	
    ny = Ny(iptr);
    nz = Nz(iptr);
    if (hasvalue("tab")) tabstr = stropen(getparam("tab"),"w");
    sigma = getrparam("sigma");
    nsigma = getrparam("nsigma");
    cutoff = nsigma * sigma;
    ns = nemoinpi(getparam("z"),seg,2);

    Qpix = (hasvalue("x") && hasvalue("y"));
    if (Qpix) {
      x = getiparam("x");
      y = getiparam("y");
    }

    sp = (double *) allocate(nz*sizeof(double));
    w  = (int *) allocate(nz*sizeof(int));
    s0 = (int *) allocate(nz*sizeof(int));
    s1 = (int *) allocate(nz*sizeof(int));

    if (Qpix) {
      printf("# Searching @ %d,%d above %g\n",x,y,cutoff);
      for (k=0; k<nz; k++) {
	sp[k] = CubeValue(iptr,x,y,k);
      }
      nseg = line_segments(nz, sp, w, s0, s1, cutoff,minchan,maxgap);
      printf("# Found %d segments\n",nseg);
      for (k=0; k<nseg; k++)
	printf("# seg %d %d\n", s0[k], s1[k]);
      for (k=0; k<nz; k++) 
	printf("%d %d %g\n",k,w[k],sp[k]);
    } else {
      ndet = 0;
      for (i=0; i<nx; i++) {
	for (j=0; j<ny; j++) {
	  for (k=0; k<nz; k++)
	    sp[k] = CubeValue(iptr,i,j,k);
	  nseg = line_segments(nz, sp, w, s0, s1, cutoff,minchan,maxgap);
	  if (nseg>0) {
	    ndet++;
	    dprintf(1,"%d %d  %d\n",i,j,nseg);
	  }
	  print_segments(i,j,nz,sp,nseg,s0,s1);
	}
      }
      dprintf(0,"Found %d/%d points with segments\n",ndet,nx*ny);
    }
}

int widx(int n, int *w, int start, int value)
{
  int i = start;
  while (i<n) {
    if (w[i] == value) return i;
    i++;
  }
  return -1;
}

void wset(int n, int *w, int start, int end, int value)
{
  int i;
  for(i=start; i<end; i++)
    w[i] = value;
}

line_segments(int n, real *spec, int *w, int *s0, int *s1, real cutoff, int minchan,  int maxgap)
{
  int i, i0, i1, i2, i3, il, ig, nseg;
  for (i=0; i<n; i++)
    w[i] = (spec[i] < cutoff ?   0   :  1);

  nseg = 0;

  i0 = 0;
  while (i0 >= 0) {
    i1 = widx(n,w,i0,1);
    if (i1<0) break;
    i2 = widx(n,w,i1,0);
    if (i2<0) break;
    i3 = widx(n,w,i2,1);
    if (i3<0) {
      il = i2-i1;
      if (il>=minchan) {
	s0[nseg] = i1;
	s1[nseg] = i2;
	nseg++;
      }
      break;
    }
    ig = i3-i2;
    if (ig<=maxgap) {
      wset(n,w,i2,i3,1);
      i0 = i1;
      continue;
    } else {
      il = i2-i1;
      if (il>=minchan) {
	s0[nseg] = i1;
	s1[nseg] = i2;
	nseg++;
      }
      i0 = i2;
      continue;
    }
  }
  return nseg;
}

print_segments(int ix, int iy, int nz, real *sp, int nseg, int *s0, int *s1)
{
  int i, k, klen;
  real sum0, sum1, sum2, peak;
  
  for (i=0; i<nseg; i++) {
    sum0 = sum1 = sum2 = peak = 0.0;
    klen = s1[i]-s0[i]+1;
    //for (k=s0[i]; k<=s1[i]; k++) {
    for (k=s0[i]; k<s1[i]; k++) {
      printf("# %d %f\n",k,sp[k]);
      sum0 += sp[k];
      sum1 += k*sp[k];
      sum2 += k*k*sp[k];
      if (peak < sp[k]) peak = sp[k];
    }
    sum1 = sum1/sum0;     /* <k>  */
    sum2 = sum2/sum0 - sum1*sum1;
    if (sum2<0) sum2 = 0.0;
    sum2 = sqrt(sum2);    /* <kw>  */
    sum2 = 2.354820 * sum2;  /* convert to FWHM */
    printf("%g %g %g %d %d %d\n",peak,sum1,sum2,klen,ix,iy);
  }
}
