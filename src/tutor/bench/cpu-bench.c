
/*
 * CPU-BENCH:  benchmark some floating point operations
 *             a.k.a. mytime in starlab
 */

#include <nemo.h>


typedef int (*action)(int);

#define NTEST 10
#define NMAX  30000

/* Work vectors: */

real v1[NMAX], v2[NMAX], v3[NMAX], v4[NMAX];
int  iv1[NMAX], iv2[NMAX], iv3[NMAX];
int  ind1[NMAX], ind2[NMAX];

real s1, s2;
real x = 1.0, y = -2.5, z = 3.14;

#include "routines.c"

/*------------------------------------------------------------------------*/

string defv[] = {
  "test=0\n          test to run (0=all)",
  "n=\n              max size of arrays",
  "VERSION=2.0\n     17-feb-2004 PJT",
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

void initialize()
{
    int i;
    real jreal = 58427.0;	/* Random numbers! */

    reset_v1(NMAX);

    for (i = 0; i < NMAX; i++) {

      v2[i] = 1.0 + 0.00001/(i+1); /* Make sure v2 > 0 for sqrt and / tests,
				      and close to 1 for repeated multiplies. */
      v3[i] = v2[i] - v1[i];

      iv1[i] = i + 42;
      iv2[i] = -2*i*i;
      iv2[i] = i*i*i - 1000;

      jreal *= 48828125.0;
      while (jreal >= 2147483648.0) jreal /= 2147483648.0;

      ind1[i] = i;
      ind2[i]  = (int)(NMAX*jreal/2147483648.0);
      if (ind2[i] < 0 || ind2[i] >= NMAX) {
	  printf("error: ind2[%d] = %d\n", i, ind2[i]);
	  exit(1);
      }

      v4[i] = jreal/2147483648.0 - 0.5;

  }

  s1 = v1[NMAX/2];
  s2 = v2[NMAX/2];

  cpu_init();
  printf("sizeof(real) = %d\n", sizeof(real));
}

int dummy(int n)
{
  return n;
}

void time_action(action a, int factor, char *name, int ntmax)
{
  real dtime0, dtime[NTEST], t_start, t_end;
  int i, j, maxtst;

  int vlen[NTEST]  = {0, 3, 10, 30, 100, 300, 1000, 3000, 30000, 300000};
  int count[NTEST] = {10000000, 3000000, 1000000, 300000, 100000,
		      30000, 10000, 3000, 300, 30};

  printf("\n%s\n", name);

  for (i = 0; i < NTEST; i++)
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
    for (j = 0; j < count[i]; j++) {
      s1 = 1./s1;
      a(vlen[i]);
    }
    t_end = cpu_time();
    dtime[i] = (t_end - t_start) / count[i] - dtime0;

    printf("vlen = %6d, time = %.3e, time/loop = %.3e, Mflops = %.4f\n",
	   vlen[i], t_end - t_start, dtime[i],
	   vlen[i]*1.0e-6 / dtime[i] * factor);
  }

  printf("function overhead = %.3e/call\n", dtime0);

  fflush(stdout);

}
/*------------------------------------------------------------------------*/

/* Routines to time now moved to separate file "routines.c"... */

/*------------------------------------------------------------------------*/

nemo_main()
{
    int n = NMAX;

    if (hasvalue("n")) n = getiparam("n");

    initialize();

    time_action(itor, 1, "vr = vi", n);
    time_action(rtoi, 1, "vi = vr", n);
    time_action(iadd, 1, "vi1 = vi2 + vi3", n);

    time_action(vmove, 1, "v1 = v2", n);

    s1 = 1.0, time_action(ssum1, 1, "s += v", n);
    s1 = 1.0, time_action(ssum2, 2, "s += v + v", n);
    s1 = 1.0, time_action(ssum3, 2, "s += v * v", n);

    s1 = 1.00000123;
    reset_v1(n);
    time_action(vsadd1, 1, "v += s", n);
    reset_v1(n);
    time_action(vsmul1, 1, "v *= s", n);

    reset_v1(n);
    time_action(vsmul1a, 1, "v *= s (partly unrolled)", n);
    time_action(vsadd2, 1, "v1 = v2 + s", n);
    time_action(vsmul2, 1, "v1 = v2 * s", n);
    time_action(vsdiv2, 1, "v1 = v2 / s", n);

    reset_v1(n);
    time_action(vsum1, 1, "v1 += v2", n);
    time_action(vsum2, 1, "v1 = v2 + v3", n);

    reset_v1(n);
    time_action(vmul1, 1, "v1 *= v2", n);
    time_action(vmul2, 1, "v1 = v2 * v3", n);
    reset_v1(n);
    time_action(vdiv1, 1, "v1 /= v2", n);
    time_action(vdiv2, 1, "v1 = v3 / v2", n);

/*     {int i; for (i = 0; i < n; i+=1000) printf("%d %f\n", i, v1[i]);} */

    time_action(saxpy1, 2, "v1 = s*v2 + s", n);
    time_action(saxpy2, 2, "v1 = s*v2 + v3", n);
    time_action(saxpy3, 2, "v1 = v2*v3 + v4", n);
    time_action(vforce, 18, "v1 = pseudo_force(v2, v3, v4)", n);

    time_action(vsqrt, 1, "v1 = sqrt(v2)", n);
    time_action(vabs, 1, "v1 = abs(v2)", n);
    time_action(vsin, 1, "v1 = sin(v2)", n);
    time_action(vexp, 1, "v1 = exp(v2)", n);
    time_action(vpow, 1, "v1 = pow(v2, s)", n);

    time_action(scatter1, 1, "v1[iv] = v2", n);
    time_action(scatter2, 1, "v1[iv] = v2", n);
    time_action(gather1, 1, "v1 = v2[iv]", n);
    time_action(gather2, 1, "v1 = v2[iv]", n);

    time_action(vif, 1, "if (v2 > 0) v1 = v2", n);
}
