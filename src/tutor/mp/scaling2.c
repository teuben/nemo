
// See scaling.c where the original program is kept

//  @todo   negotiate between using OMP_NUM_THREADS and np=

#include <nemo.h>
#include <stdint.h>

string defv[] = {
    "umin=1\n    Starting value",
    "umax=10000000000\n   Ending",    // not working
    "VERSION=1\n  24-jan-2021 PJT",
     NULL,
};

string usage="NEMO version of the scaling program";
		 

void nemo_main(void)
{
  const uint64_t umin=getlparam("umin");
  const uint64_t umax=10000000000LL;
  //const uint64_t umax=getlparam("umax");
  double sum=0.0;
  dprintf(0,"scaling2: umin=%ld umax=%ld\n",umin,umax);
#pragma omp parallel for reduction(+:sum)
  for(uint64_t u=umin; u<umax; u++)
    sum+=1./(u*u);   //  1./u/u takes about 2x longer on intel
  printf("%g\n", sum);

}
