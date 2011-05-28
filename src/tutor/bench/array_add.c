/*
 * benchmark openmp vs. pthreads
 * 
 * example taken from
 * http://www.futurechips.org/tips-for-power-coders/open-mp-pthreads.html
 *
 *     28-mar-2011    Created  // PJT
 */
 

#include <nemo.h>



string defv[] = {
  "n=1000\n       Number of times to execute sum_st",
  "VERSION=1\n    28-may-2011",
  NULL,
};

string usage = "bench openmp";

#define M 10000


void sum_st(int *A, int *B, int *C)
{
  int i;
#pragma omp parallel for
  for(i = 0; i < M; i++)
    A[i] = B[i] + C[i];
}



/*  pthreads version */
#if 0
struct params {
  int *A;
  int *B;
  int *C;
  int tid;
  int size;
  int nthreads;
};

void *compute_parallel(void *_p){
  params *p      = (params*) _p;
  int tid        = p->tid;
  int chunk_size = (p->size / p->nthreads);
  int start      = tid * chunk_size;
  int end        = start + chunk_size;
  for(int i = start; i < end; i++)     p->A[i] = p->B[i] + p->C[i];
  return 0;
}

void sum_mt(int *A, int *B, int *C){
  int nthreads = 4;
  int size = 10000;
  pthread_t threads[nthreads]; //array to hold thread information
  params *thread_params = (params*) malloc(nthreads * sizeof(params));

  for(int i = 0; i < nthreads; i++){
    thread_params[i].A        = A;
    thread_params[i].B        = B;
    thread_params[i].C        = C;
    thread_params[i].tid      = i;
    thread_params[i].size     = size;
    thread_params[i].nthreads = nthreads;
    pthread_create(&threads[i], NULL, compute_parallel, (void*) &thread_params[i]);
  }

  for(int i = 0; i < nthreads; i++){
    pthread_join(threads[i], NULL);
  }
  free(thread_params);

}
#endif

nemo_main() 
{
  int i, n = getiparam("n");
  int *A = (int *) allocate(sizeof(int)*M);
  int *B = (int *) allocate(sizeof(int)*M);
  int *C = (int *) allocate(sizeof(int)*M);

  for (i=0; i<M; i++) {
    A[i] = xrandom(0.0,M);
    B[i] = xrandom(0.0,M);
  }

  while(n--)
    sum_st(A,B,C);
}

