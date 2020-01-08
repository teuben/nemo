/*
 * sections:    possibly the easiest to OMP ?   careful about
                the work call. they might need to mutex
 */

#include <nemo.h>

string defv[] = {
  "n=1000000\n      array size will be square of this",
  "iter=10\n        how many times to iterate each block",
  "seed=123\n       seed for xrandom",
  "mode=0\n         Different openmp experiments",
  "VERSION=0.2\n    7-jan-2020 PJT",
  NULL,
};

string usage="benchmark openmp overhead of starting a tast via sections";

string cvsid="$Id:$";

void work1(int n, double *x, int iter);
void work2(int n, double *x, int iter);
void work1r(int n, double *x, int iter);
void work2r(int n, double *x, int iter);

void nemo_main()
{
  int n = getiparam("n");
  int seed = init_xrandom(getparam("seed"));
  int iter1 = getiparam("iter");
  int iter2 = iter1;
  int i;
  int n1 = n;
  int n2 = n;
  int mode = getiparam("mode");
  real *x1 = (real *) allocate(n1*sizeof(double));
  real *x2 = (real *) allocate(n2*sizeof(double));

  dprintf(0,"n=%d iter=%d\n",n,iter1);
  // cannot omp this for loop: xrandom has no mutex
  for (i=0; i<n1; ++i) {
    x1[i] = xrandom(0.0,1.0);
    x2[i] = xrandom(0.0,1.0);
  }

  if (mode==0) {
    dprintf(0,"mode=0: vanilla sections\n");    
    #pragma omp parallel sections
    {
      work1(n1,x1,iter1);
      #pragma omp section
      work2(n2,x2,iter2);
    }
  } else if (mode == 1) {
    dprintf(0,"mode=1: just work1r\n");
    work1r(n1,x1,iter1);
  } else {
    dprintf(0,"mode=2: reduc on work1r and work2r\n");
    #pragma omp parallel sections
    {
      work1r(n1,x1,iter1);
      #pragma omp section
      work2r(n2,x2,iter2);
    }

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

void work1r(int n, real *x, int iter)
{
  printf("work1 %d\n",iter);
  int i;
  real sum = 0.0;

  while (iter--)
    #pragma omp parallel shared(n, x, sum) private(i)
    #pragma omp for reduction(+:sum)   
    for (i=0; i<n; ++i)
      sum += x[i];
  printf("sum1=%g\n",sum);
}

void work2r(int n, real *x, int iter)
{
  printf("work2 %d\n",iter);
  int i;
  real sum = 0.0;

  while (iter--)
    #pragma omp parallel shared(n, x, sum) private(i)
    #pragma omp for reduction(+:sum)   
    for (i=0; i<n; ++i)
      sum += x[i];      
  printf("sum2=%g\n",sum);
}
