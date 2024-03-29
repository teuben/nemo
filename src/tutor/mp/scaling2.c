
// See scaling.c where the original program is kept

//  @todo   negotiate between using OMP_NUM_THREADS and np=
//
//   Processor affinity can affect performance?
//   
//   export GOMP_CPU_AFFINITY="0 1 2 3 4 5 6 7"
//   /usr/bin/time scaling2 1 200000 np=4
//   75.01user 0.00system 0:18.81elapsed 398%CPU
//   75.03user 0.00system 0:18.81elapsed 398%CPU
//   75.03user 0.00system 0:18.81elapsed 398%CPU 
//   unset GOMP_CPU_AFFINITY
//   75.69user 0.02system 0:19.04elapsed 397%CPU
//   75.51user 0.00system 0:18.95elapsed 398%CPU
//   75.65user 0.03system 0:19.01elapsed 397%CPU 


#include <nemo.h>
#include <stdint.h>

string defv[] = {
    "umin=1\n             Starting value",
    "umax=10000\n         sqrt of Ending value",   
    "umax2=0\n            Non-parallel loop, again sqrt of Ending value",
    "iter=1\n             How many times to iterate and report timing",
    "VERSION=1.4\n        17-sep-2023 PJT",
     NULL,
};

string usage="NEMO version of the well scaled OpenMP scaling program";
		 

void nemo_main(void)
{
  int umin4 = getiparam("umin");  // NEMO doesn't have a reliable longlong
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
      //sum+=1./u/u;       //  1./u/u takes about 2x longer on intel but doesn't underflow
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
    dprintf(0,"cputime: %d %g %g %g sec\n", niter+1, t0, t1, t1-t2);
    t2 = t1;
  } // niter
}
