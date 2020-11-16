// Compile with: gcc scaling.c -std=c99 -fopenmp -O3                                                                                               
// See also discussion on:
//    https://stackoverflow.com/questions/19780554/what-limits-scaling-in-this-simple-openmp-program

#include <stdio.h>
#include <stdint.h>

int main(){

  const uint64_t umin=1;
  const uint64_t umax=10000000000LL;
  double sum=0.;
#pragma omp parallel for reduction(+:sum)
  for(uint64_t u=umin; u<umax; u++)
    sum+=1./(u*u);   //  1./u/u takes about 2x longer on intel
  printf("%e\n", sum);

}
