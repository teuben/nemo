/*
 * randomizer.c:   select M out of N

 */

#include <stdinc.h>


int *randomizer(int n, int m)
{
  int *sel, *idx;
  int i,j,k;

  if (m > n) error("Randomizer: m=%d >= n=%d",m,n);

  idx = (int *) allocate(n*sizeof(int));
  sel = (int *) allocate(m*sizeof(int));

  for (i=0; i<n; i++)
    idx[i] = i;

  for (i=0, k=0;  i<m;  i++, k++) {
    j = (int) xrandom(0.0,n-i);
    sel[k] = idx[j];
    idx[j] = idx[n-1-i];
  }
  free(idx);
  return sel;
}

#ifdef TESTBED

#include <getparam.h>

string defv[] = {
  "n=10\n	      How many in sample",
  "m=5\n              How many to pick",
  "seed=0\n           Seed",
  "VERSION=1.0\n      28-feb-08 PJT",
  NULL,
};

string usage="test randomizer()";

nemo_main()
{
  int i;
  int n=getiparam("n");
  int m=getiparam("m");
  int *idx;
  int seed;

  seed = init_xrandom(getparam("seed"));
  dprintf(1,"seed=%d\n",seed);

  idx = randomizer(n,m);
  
  for (i=0; i<m; i++)
    printf("%d\n",idx[i]);

  //  free(idx);
}

#endif
