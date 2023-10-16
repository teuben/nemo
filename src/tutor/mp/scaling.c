// Compile with: gcc scaling.c -std=c99 -fopenmp -O3                                                                                               
// See also discussion on:
//    https://stackoverflow.com/questions/19780554/what-limits-scaling-in-this-simple-openmp-program

#include <stdio.h>
#include <stdint.h>
#include <math.h>

int main(){

  const uint64_t umin=1;
  const uint64_t umax=4000000000LL;    //   4->5 already causes overflow
  double sum=0.;
#pragma omp parallel for reduction(+:sum)
  for(uint64_t u=umin; u<umax; u++) {
    sum+=1./(u*u);       //    10.0"  1./u/u takes about 2x longer on intel
    //sum+=1./u/u;         //  15.6"  
    //sum+=1/pow(u,2.0);   // 8.6"  compiler optimizes this out!!!
    //sum+=1/pow(u,1.33);  // 15.3"
    // sum+=(1/u)*(1/u);   //  18.5"
  }
  printf("%e\n", sum);

}
