
// See scaling.c where the original program is kept

//  @todo   negotiate between using OMP_NUM_THREADS and np=

#include <nemo.h>
#include <stdint.h>

string defv[] = {
    "umin=1\n             Starting value",
    "umax=10000\n         sqrt of Ending",   
    "umax2=0\n            Non-parallel loop",
    "iter=1\n             How many times to iterate and report timing",
    "VERSION=1.2\n        29-aug-2022 PJT",
     NULL,
};

string usage="NEMO version of the OpenMP scaling program";
		 

void nemo_main(void)
{
  int umin4 = getiparam("umin");
  int umax4 = getiparam("umax");
  int umax2 = getiparam("umax2");
  int niter = getiparam("iter");
  uint64_t umin = (uint64_t) umin4 * (uint64_t) umin4;
  uint64_t umax = (uint64_t) umax4 * (uint64_t) umax4;
  double t0, t1, t2 = 0.0;
  bool Qshow = niter == 1;
  extern int np_openmp;  // this is a cheat; see getparam.c

  dprintf(0,"omp_get_num_procs() -> %d\n",np_openmp);
  dprintf(0,"scaling2: umin=%ld umax=%ld\n",umin,umax);
  

  while (niter--) {
    double sum=0.0;
#pragma omp parallel for reduction(+:sum)
    for(uint64_t u=umin; u<umax; u++)
      sum+=1./(u*u);   //  1./u/u takes about 2x longer on intel
    if (Qshow) printf("sum=%g\n", sum);
    if (umax2 > 0) {
      sum=0.0;
      umax = (uint64_t) umax2 * (uint64_t) umax2;
      for(uint64_t u=umin; u<umax; u++)
	sum+=1./(u*u); 
      if (Qshow) printf("sum2=%g\n", sum);	
    }
    t0 = 60*cputime2(0);
    t1 = 60*cputime2(2);
    dprintf(0,"cputime: %d %g %g %g sec\n",iter, t0, t1, t1-t2);
    t2 = t1;
  }
  

}
