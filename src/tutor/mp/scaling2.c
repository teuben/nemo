
// See scaling.c where the original program is kept

//  @todo   negotiate between using OMP_NUM_THREADS and np=

#include <nemo.h>
#include <stdint.h>

string defv[] = {
    "umin=1\n             Starting value",
    "umax=10000\n          sqrt of Ending",   
    "umax2=0\n            Non-parallel loop",
    "VERSION=1.1\n        13-mar-2021 PJT",
     NULL,
};

string usage="NEMO version of the scaling program";
		 

void nemo_main(void)
{
  int umin4 = getiparam("umin");
  int umax4 = getiparam("umax");
  int umax2 = getiparam("umax2");
  uint64_t umin = (uint64_t) umin4 * (uint64_t) umin4;
  uint64_t umax = (uint64_t) umax4 * (uint64_t) umax4;

  double sum=0.0;
  dprintf(0,"scaling2: umin=%ld umax=%ld\n",umin,umax);
#pragma omp parallel for reduction(+:sum)
  for(uint64_t u=umin; u<umax; u++)
    sum+=1./(u*u);   //  1./u/u takes about 2x longer on intel
  printf("sum=%g\n", sum);
  if (umax2 > 0) {
        sum=0.0;
	umax = (uint64_t) umax2 * (uint64_t) umax2;
	for(uint64_t u=umin; u<umax; u++)
	    sum+=1./(u*u); 
        printf("sum2=%g\n", sum);	
   }

}
