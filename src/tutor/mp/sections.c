/*
 * sections:    possibly the easiest to OMP ?   careful about
                the work call. they might need to mutex
 */

#include <nemo.h>

string defv[] = {
  "nmax=1000\n         array size will be square of this",
  "iter=10\n           how many times to iterate each block",
  "seed=123\n          seed for xrandom",
  "VERSION=0.1\n       7-jan-2020 PJT",
  NULL,
};

string usage="benchmark openmp overhead of starting a tast via sections";

string cvsid="$Id:$";

void work1(int n, double *x, int iter);
void work2(int n, double *x, int iter);

void nemo_main()
{
  int n = getiparam("nmax");
  int seed = init_xrandom(getparam("seed"));
  int iter1 = getiparam("iter");
  int iter2 = iter1;
  int i;
  int n1 = n*n;
  int n2 = n*n;
  real *x1 = (real *) allocate(n1*sizeof(double));
  real *x2 = (real *) allocate(n2*sizeof(double));
  // cannot omp this for loop: xrandom has no mutex
  for (i=0; i<n1; ++i) {
    x1[i] = xrandom(0.0,1.0);
    x2[i] = xrandom(0.0,1.0);
  }
  
  #pragma omp parallel sections
  {
    work1(n1,x1,iter1);
    #pragma omp section
    work2(n2,x2,iter2);
  }
}

void work1(int n, real *x, int iter)
{
  printf("work1 %d\n",iter);
  int i;
  real sum = 0.0;

  while (iter--)
    for (i=0; i<n; ++i)
      sum += x[i];
  printf("sum1=%g\n",sum);
}

void work2(int n, real *x, int iter)
{
  printf("work2 %d\n",iter);
  int i;
  real sum = 0.0;

  while (iter--)
    for (i=0; i<n; ++i)
      sum += x[i];      
  printf("sum2=%g\n",sum);
}
