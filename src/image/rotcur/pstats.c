/*
 * simulated the noise in a profile to determine the
 * expected uncertainty in the first moment
 *
 * sigma_v = k * sqrt(N) * sigma_I * Delta V
 *
 * where 'k' is the constant to be determined.
 *
 * For (the default) box between -1 and 1  and N=2 (so Delta V = 2)
 * k = 0.25, which can be derived from first principles.
 * For N large k = 1/sqrt(12), which can also be derived analytically.
 *
 * This C program is about 10 times faster than the equivalent python script
 */

#include <nemo.h>
#include <moment.h>

string defv[] = {
    "n=10000\n        number of experiments",
    "k=2\n            number of samples between -1 and 1",
    "eps=0.01\n       normalized noise (1/SNR)",
    "gauss=-1\n       if positive, this is the gaussian sigma (=FWHM/2.355)",
    "hanning=f\n      Hanning smooth?",
    "seed=0\n         random seed",
    "VERSION=0.3\n    10-Jun-2020 XYZ",
    NULL,
};

string usage="profile stats helper";

string cvsid="$Id:$";

void hanning(int n, real *x, real *tmp);

void nemo_main()
{
  int n = getiparam("n");
  int k = getiparam("k");
  bool Qhanning = getbparam("hanning");
  real gauss = getrparam("gauss");
  real eps = getrparam("eps");
  real dx, *x, *y, xysum, ysum, *z;
  real vm, vs, k1;
  //int seed = init_xrandom(getparam("seed"));
  int seed = set_xrandom(getiparam("seed"));
  Moment m1, m2;
  int i, j;

  // gaussian?
  if (gauss > 0) gauss = 2*gauss*gauss;
  dprintf(1,"gauss %g\n",gauss);
  
  // get X array
  dx = 2.0 / (k-1);
  x = (real *) allocate(k*sizeof(real));
  y = (real *) allocate(k*sizeof(real));
  z = (real *) allocate(k*sizeof(real));
  for (j=0; j<k; j++)  x[j] = -1.0 + j*dx;

  // stats
  ini_moment(&m1, 2, 0);

  // loop
  for (i=0; i<n; i++) {
    xysum = ysum = 0.0;
    if (gauss < 0)
      for (j=0; j<k; j++)      
	y[j] = grandom(1.0,eps);
    else
      for (j=0; j<k; j++)
	y[j] = exp(-x[j]*x[j]/gauss) + grandom(0.0,eps);
    if (Qhanning)
      hanning(k, y, z);
    for (j=0; j<k; j++) {
      xysum += x[j]*y[j];
      ysum  += y[j];
      //printf("%g %g\n",x[j],y);
    }
    xysum /= ysum;
    //printf("%g\n",xysum);
    accum_moment(&m1,xysum,1.0);
  }
  // report
  vm = mean_moment(&m1);
  vs = sigma_moment(&m1);
  k1 = vs / sqrt(k) / eps / (2.0/(k-1));
  printf("%d %d %g %g  %g %g   %g %d\n",n,k,eps, getrparam("gauss"), vm, vs, k1, seed);
}



void hanning(int n, real *x, real *tmp)
{
  int i;
  for (i=0; i<n; i++)
    tmp[i] = x[i];
  for (i=1; i<n-1; i++)
    x[i] = 0.25*tmp[i-1] + 0.5*tmp[i] + 0.25*tmp[i+1];
}
