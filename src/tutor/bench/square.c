/*
 *  benchmark large matrix operations
 *
 *  I found that:
 *    - no difference if n big or small, arguing cache not important
 *    - oddly enough if iter was big i saw some degraded speeds
 *   e.q. FUN(a,b)  (1/sqrt(a+b))
 *   square 10000 50  > Mops=4.68242
 *   square 10000 10    Mops=33.2447
 *          1000 100    

 *   10 10 3  -> 56500  local
 *               56512  in library
 */

#include <nemo.h>
#include <timers.h>

string defv[] = {
    "n=1000\n        size of vector",
    "iter=10\n      number of times to repeat the benchmark",
    "m=0\n           blocksize in vector to process",
    "seed=0\n        xrandom seed",
    "VERSION=0.2\n   24-Apr-2004 XYZ",
    NULL,
};

string usage="benchmark cross product/matrix operations";

extern real xrandom(real,real);

#define FUN(a,b)  (1/sqrt(a+b))

nemo_main()
{
  int n=getiparam("n");
  int iter=getiparam("iter");
  int m=getiparam("m");
  int seed = init_xrandom(getparam("seed"));
  int i,j,k,l;
  real *x, sum;
  real t0,t1,t2;

  init_timers(100);
  stamp_timers(0);

  x = (real *) allocate(n*sizeof(real));

  for (i=0; i<n; i++)         /* init the whole array */
    x[i] = xrandom(0.0,1.0);
  for (i=0; i<m; i++)         /* cache it in again ? */
    x[i] = xrandom(0.0,1.0);

  sum = 0.0;
  t0 = cputime();
  stamp_timers(1);
  if (m==0) {                         /* do it in one sweep, the N^2 algorithm */
    stamp_timers(2);
    for (l=0; l<iter;  l++)
      for (j=0; j<n; j++)
	for (i=0; i<n; i++)
	  sum += FUN(x[i],x[j]);
    stamp_timers(3);
  } else {                            /* N/M times a small M*M patch that may be in cache */
    stamp_timers(2);
    for (l=0; l<iter; l++)
      for (k=0; k<n-m; k++)
	for (j=k; j<k+m; j++)
	  for (i=k; i<k+m; i++)
	    sum += FUN(x[i],x[j]);
    stamp_timers(3);
  }
  stamp_timers(4);
  t1 = cputime();
  if (m)
    printf("%d %d %d sum=%lg Mops=%lg\n",
	 n,iter,m,sum,iter*m*m*n/(t1-t0)/60.0/1e6);
  else
    printf("%d %d %d sum=%lg Mops=%lg\n",
	 n,iter,m,sum,iter*n*n/(t1-t0)/60.0/1e6);
  stamp_timers(5);
  printf("%Ld %Ld %Ld %Ld %Ld\n",
	 diff_timers(0,1),
	 diff_timers(1,2),
	 diff_timers(2,3),
	 diff_timers(3,4),
	 diff_timers(4,5));
}
