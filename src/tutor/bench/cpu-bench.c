
/*
 * CPU-BENCH:  benchmark some floating point operations
 *             a.k.a. mytime in starlab
 */

#include <nemo.h>
#include <timers.h>

typedef int (*action)(int);

#define MAXTEST 11
#define NMAX  300000
/* #define VSTATIC  */

/* Work vectors: */

#if defined(VSTATIC)
 real v1[NMAX], v2[NMAX], v3[NMAX], v4[NMAX];
 int  iv1[NMAX], iv2[NMAX], iv3[NMAX];
 int  ind1[NMAX], ind2[NMAX];
 int  v_size = NMAX;
#else
 real *v1, *v2, *v3, *v4;
 int  *iv1, *iv2, *iv3;
 int  *ind1, *ind2;
 int  v_size = 0;
#endif

real s1, s2;
real x = 1.0, y = -2.5, z = 3.14;

#include "routines.c"

/*------------------------------------------------------------------------*/

string defv[] = {
  "test=0\n          test to run (0=all)\n"
  "                  1   vr = vi\n"
  "                  2   vi = vr\n"
  "                  3   vi1 = vi2 + vi3\n"
  "                  4   v1 = v2\n"
  "                  5   s += v\n"
  "                  6   s += v + v\n"
  "                  7   s += v * v\n"
  "                  8   v += s\n"
  "                  9   v *= s\n"
  "                  10  v *= s (partly unrolled)\n"
  "                  11  v1 = v2 + s\n"
  "                  12  v1 = v2 * s\n"
  "                  13  v1 = v2 / s\n"
  "                  14  v1 += v2\n"
  "                  15  v1 = v2 + v3\n"
  "                  16  v1 *= v2\n"
  "                  17  v1 = v2 * v3\n"
  "                  18  v1 /= v2\n"
  "                  19  v1 = v3 / v2\n"
  "                  20  v1 = s*v2 + s\n"
  "                  21  v1 = s*v2 + v3\n"
  "                  22  v1 = v2*v3 + v4\n"
  "                  23  v1 = pseudo_force(v2, v3, v4)\n"
  "                  24  v1 = sqrt(v2)\n"
  "                  25  v1 = abs(v2)\n"
  "                  26  v1 = sin(v2)\n"
  "                  27  v1 = exp(v2)\n"
  "                  28  v1 = pow(v2, s)\n"
  "                  29  v1[iv] = v2\n"
  "                  30  v1[iv] = v2\n"
  "                  31  v1 = v2[iv]\n"
  "                  32  v1 = v2[iv]\n"
  "                  33  if (v2 > 0) v1 = v2",
  "n=30000\n         max size of arrays",
  "vlen=0,3,10,30,100,300,1000,3000,30000\n     Vector lengths to try",
  "count=300000\n    how often shoudl the largest vlen be run",
  "VERSION=3.1\n     24-apr-2004 PJT",
  NULL,
};

string usage = "cpu benchmark";

/*------------------------------------------------------------------------*/

#include <sys/times.h>
#include <unistd.h>

struct tms  buffer;
static long ticks_per_sec = 0;
static real initial_cpu = 0.0;

void cpu_init()
{
  times(&buffer);

  /* Use both system and user time because of ambiguities
     with Linux multiprocessing... */

    initial_cpu = (real) (buffer.tms_utime + buffer.tms_stime);

    ticks_per_sec = sysconf(_SC_CLK_TCK);	/* Clock ticks per second */
}

real cpu_time()
{
    if (!ticks_per_sec) cpu_init();

    times(&buffer);
    return ((real) (buffer.tms_utime + buffer.tms_stime - initial_cpu))
      			/ ticks_per_sec;
}

/*------------------------------------------------------------------------*/

void reset_v1(int l)
{
    int i;
    for (i = 0; i < l; i++)
	v1[i] = i / (i + 1.0);
}

void initialize(int n)
{
    int i;
    real jreal = 58427.0;	/* Random numbers! */

#if !defined(VSTATIC)
    if (v_size == 0) {
      v1 = (real *) allocate(n*sizeof(real));
      v2 = (real *) allocate(n*sizeof(real));
      v3 = (real *) allocate(n*sizeof(real));
      v4 = (real *) allocate(n*sizeof(real));

      iv1 = (int *) allocate(n*sizeof(int));
      iv2 = (int *) allocate(n*sizeof(int));
      iv3 = (int *) allocate(n*sizeof(int));
      ind1 = (int *) allocate(n*sizeof(int));
      ind2 = (int *) allocate(n*sizeof(int));
      v_size = n;
    }
#endif

    reset_v1(n);

    for (i = 0; i < n; i++) {

      v2[i] = 1.0 + 0.00001/(i+1); /* Make sure v2 > 0 for sqrt and / tests,
				      and close to 1 for repeated multiplies. */
      v3[i] = v2[i] - v1[i];

      iv1[i] = i + 42;
      iv2[i] = -2*i*i;
      iv2[i] = i*i*i - 1000;

      jreal *= 48828125.0;
      while (jreal >= 2147483648.0) jreal /= 2147483648.0;

      ind1[i] = i;
      ind2[i]  = (int)(n*jreal/2147483648.0);
      if (ind2[i] < 0 || ind2[i] >= n)
	  error("ind2[%d] = %d\n", i, ind2[i]);

      v4[i] = jreal/2147483648.0 - 0.5;

  }

  s1 = v1[n/2];
  s2 = v2[n/2];

  cpu_init();
  printf("sizeof(real) = %d\n", sizeof(real));
#if defined(VSTATIC)
  printf("Using static arrays\n");
#else
  printf("Using dynamic arrays\n");
#endif
}

int dummy(int n)
{
  return n;
}

void new_time_action(action a, int factor, char *name, int ntmax)
{
  static int ntest=MAXTEST;
  real dtime0, dtime[MAXTEST], t_start, t_end;
  int i, j, maxtst;

  int vlen[MAXTEST]  = {0, 3, 10, 30, 100, 300, 1000, 3000, 30000, 300000, 3000000};
  int count[MAXTEST] = {10000000, 3000000, 1000000, 300000, 100000,
		      30000, 10000, 3000, 300, 30, 3};

  printf("\n%s\n", name);

  if (ntest==0) {
    error("new code");
    ntest = nemoinpi(getparam("vlen"),vlen,MAXTEST);
    if (ntest<=0) error("parsing");
    
  }

  for (i = 0; i < ntest; i++)
    if (vlen[i] <= ntmax) maxtst = i;

  /* Function call timing. */

  i = 0;
  t_start = cpu_time();
  for (j = 0; j < count[i]; j++) dummy(vlen[i]);
  t_end = cpu_time();
  dtime0 = (t_end - t_start) / count[i];

  for (i = 0; i <= maxtst; i++) {
    real sum;

    t_start = cpu_time();
    stamp_timers(i);
    for (j = 0; j < count[i]; j++) {
      s1 = 1./s1;
      a(vlen[i]);
    }
    stamp_timers(i+1);
    t_end = cpu_time();
    dtime[i] = (t_end - t_start) / count[i] - dtime0;
#if 1
    printf("vlen = %6d, time = %.3e, time/loop = %.3e, Mflops = %.4f Ticks = %Ld\n",
	   vlen[i], t_end - t_start, dtime[i],
	   vlen[i]*1.0e-6 / dtime[i] * factor,
	   diff_timers(i,i+1));
#endif
  } /* i */
  
  for (i = 0; i <= maxtst; i++) {
    printf("vlen = %6d,  time/loop = %.3e, Mflops = %.4f Ticks = %Ld\n",
	   vlen[i], dtime[i],
	   vlen[i]*1.0e-6 / dtime[i] * factor,
	   diff_timers(i,i+1));

  }

  printf("function overhead = %.3e/call\n", dtime0);

  fflush(stdout);

}

void time_action(action a, int factor, char *name, int ntmax)
{
  static int ntest=MAXTEST;
  real dtime0, dtime[MAXTEST], t_start, t_end;
  int i, j, maxtst;

  int vlen[MAXTEST]  = {0, 3, 10, 30, 100, 300, 1000, 3000, 30000, 300000, 3000000};
  int count[MAXTEST] = {10000000, 3000000, 1000000, 300000, 100000,
		      30000, 10000, 3000, 300, 30, 3};

  printf("\n%s\n", name);

  if (ntest==0) {
    error("new code");
    ntest = nemoinpi(getparam("vlen"),vlen,MAXTEST);
    if (ntest<=0) error("parsing");
    
  }

  for (i = 0; i < ntest; i++)
    if (vlen[i] <= ntmax) maxtst = i;

  /* Function call timing. */

  i = 0;
  t_start = cpu_time();
  for (j = 0; j < count[i]; j++) dummy(vlen[i]);
  t_end = cpu_time();
  dtime0 = (t_end - t_start) / count[i];

  for (i = 0; i <= maxtst; i++) {
    real sum;

    stamp_timers(i);
    for (j = 0; j < count[i]; j++) {
      s1 = 1./s1;
      a(vlen[i]);
    }
    stamp_timers(i+1);
  } /* i */

  for (i = 0; i <= maxtst; i++) {
    printf("vlen = %6d,  Ticks = %Ld\n",
	   vlen[i], diff_timers(i,i+1));

  }


}  

void old_time_action(action a, int factor, char *name, int ntmax)
{
  static int ntest=MAXTEST;
  real dtime0, dtime[MAXTEST], t_start, t_end;
  int i, j, maxtst;

  int vlen[MAXTEST]  = {0, 3, 10, 30, 100, 300, 1000, 3000, 30000, 300000, 3000000};
  int count[MAXTEST] = {10000000, 3000000, 1000000, 300000, 100000,
		      30000, 10000, 3000, 300, 30, 3};

  printf("\n%s\n", name);

  if (ntest==0) {
    error("new code");
    ntest = nemoinpi(getparam("vlen"),vlen,MAXTEST);
    if (ntest<=0) error("parsing");
    
  }

  for (i = 0; i < ntest; i++)
    if (vlen[i] <= ntmax) maxtst = i;

  /* Function call timing. */

  i = 0;
  t_start = cpu_time();
  for (j = 0; j < count[i]; j++) dummy(vlen[i]);
  t_end = cpu_time();
  dtime0 = (t_end - t_start) / count[i];

  for (i = 0; i <= maxtst; i++) {
    real sum;

    t_start = cpu_time();
    stamp_timers(i);
    for (j = 0; j < count[i]; j++) {
      s1 = 1./s1;
      a(vlen[i]);
    }
    stamp_timers(i+1);
    t_end = cpu_time();
    dtime[i] = (t_end - t_start) / count[i] - dtime0;
#if 1
    printf("vlen = %6d, time = %.3e, time/loop = %.3e, Mflops = %.4f Ticks = %Ld\n",
	   vlen[i], t_end - t_start, dtime[i],
	   vlen[i]*1.0e-6 / dtime[i] * factor,
	   diff_timers(i,i+1));
#endif
  } /* i */
  
  for (i = 0; i <= maxtst; i++) {
    printf("vlen = %6d,  Ticks = %Ld\n",
	   vlen[i], diff_timers(i,i+1));

  }


  fflush(stdout);

}

int do_test(a,b) {
  stamp_timers(b);
  if (a==0 || b==0) return 1;
  if (a==b) return 1;
  return 0;
}

/*------------------------------------------------------------------------*/

/* Routines to time now moved to separate file "routines.c"... */

/*------------------------------------------------------------------------*/

nemo_main()
{
    int i, n = NMAX;
    int test = getiparam("test");

    if (hasvalue("n")) n = getiparam("n");
#if defined(VSTATIC)
    if (n > NMAX) {
      warning("n=%d , NMAX=%d cannot be exceeded",n,NMAX);
      n = NMAX;
    }
#endif

    initialize(n);
    init_timers(200);
    stamp_timers(0);
    if (do_test(test,1)) time_action(itor, 1, "vr = vi", n);
    if (do_test(test,2)) time_action(rtoi, 1, "vi = vr", n);
    if (do_test(test,3)) time_action(iadd, 1, "vi1 = vi2 + vi3", n);

    if (do_test(test,4)) time_action(vmove, 1, "v1 = v2", n);

    s1 = 1.0; if (do_test(test,5)) time_action(ssum1, 1, "s += v", n);
    s1 = 1.0; if (do_test(test,6)) time_action(ssum2, 2, "s += v + v", n);
    s1 = 1.0; if (do_test(test,7)) time_action(ssum3, 2, "s += v * v", n);

    s1 = 1.00000123;
    reset_v1(n);
    if (do_test(test,8)) time_action(vsadd1, 1, "v += s", n);
    reset_v1(n);
    if (do_test(test,9)) time_action(vsmul1, 1, "v *= s", n);

    reset_v1(n);
    if (do_test(test,10)) time_action(vsmul1a, 1, "v *= s (partly unrolled)", n);
    if (do_test(test,11)) time_action(vsadd2, 1, "v1 = v2 + s", n);
    if (do_test(test,12)) time_action(vsmul2, 1, "v1 = v2 * s", n);
    if (do_test(test,13)) time_action(vsdiv2, 1, "v1 = v2 / s", n);

    reset_v1(n);
    if (do_test(test,14)) time_action(vsum1, 1, "v1 += v2", n);
    if (do_test(test,15)) time_action(vsum2, 1, "v1 = v2 + v3", n);

    reset_v1(n);
    if (do_test(test,16)) time_action(vmul1, 1, "v1 *= v2", n);
    if (do_test(test,17)) time_action(vmul2, 1, "v1 = v2 * v3", n);
    reset_v1(n);
    if (do_test(test,18)) time_action(vdiv1, 1, "v1 /= v2", n);
    if (do_test(test,19)) time_action(vdiv2, 1, "v1 = v3 / v2", n);

/*     {int i; for (i = 0; i < n; i+=1000) printf("%d %f\n", i, v1[i]);} */

    if (do_test(test,20)) time_action(saxpy1, 2, "v1 = s*v2 + s", n);
    if (do_test(test,21)) time_action(saxpy2, 2, "v1 = s*v2 + v3", n);
    if (do_test(test,22)) time_action(saxpy3, 2, "v1 = v2*v3 + v4", n);
    if (do_test(test,23)) time_action(vforce, 18, "v1 = pseudo_force(v2, v3, v4)", n);

    if (do_test(test,24)) time_action(vsqrt, 1, "v1 = sqrt(v2)", n);
    if (do_test(test,25)) time_action(vabs, 1, "v1 = abs(v2)", n);
    if (do_test(test,26)) time_action(vsin, 1, "v1 = sin(v2)", n);
    if (do_test(test,27)) time_action(vexp, 1, "v1 = exp(v2)", n);
    if (do_test(test,28)) time_action(vpow, 1, "v1 = pow(v2, s)", n);

    if (do_test(test,29)) time_action(scatter1, 1, "v1[iv] = v2", n);
    if (do_test(test,30)) time_action(scatter2, 1, "v1[iv] = v2", n);
    if (do_test(test,31)) time_action(gather1, 1, "v1 = v2[iv]", n);
    if (do_test(test,32)) time_action(gather2, 1, "v1 = v2[iv]", n);
    if (do_test(test,33)) time_action(vif, 1, "if (v2 > 0) v1 = v2", n);
    if (test)
        printf("Ticks: %Ld\n",diff_timers(test,test+1));
    else
	for (i=0; i<33; i++)
	    printf("Ticks %d %Ld:\n",i+1,diff_timers(i,i+1));
}
