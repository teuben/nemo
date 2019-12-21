/*
 *   Simple Example OMP usage
 */

#include <omp.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main (void)
{
  int nthreads, my_id;
  int i, imax = 1e8,j,jmax;
  double x=0;

#ifdef _OPENMP
  printf("Compiled by an OpenMP-compliant implementation.\n");

  #pragma omp parallel  shared(x)
  { 
    nthreads = omp_get_num_threads();
    my_id = omp_get_thread_num();
    for (i=0;i<imax/nthreads;i++)
      x += cos(i*x);
  }
  printf("omp: making %d cosines to calculate x= %3.3e took.. \n",imax,x);

#else
  printf("Normal compilation\n");
    
  for (i=0;i<imax;i++)
    x += cos(i*x);
  printf("cpu1: making %d cosines to calculate x= %3.3e took.. \n",imax,x);
#endif
    
 
}
