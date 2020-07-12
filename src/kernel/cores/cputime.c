/*
 * CPUTIME: return CPU time in minutes
 * ETIME, DTIME: fortran support
 *
 *      28-nov-89   non-VMS version uses HZ variable now
 *                  on most machines (in the us) HZ will be 60
 *		    For SUN and BSD it is defined in
 *			<param.h>
 *	26-may-91   re-introduced the correct Cray code + TOOLBOX from starlab
 *		    Whenever we have POSIX compliant code - this
 *		    nasty piece will be better code...
 *	 4-nov-91   Multiflow _trace_ seems to have funny HZ
 *		    better to use getrusage(2) - Cray doesn't know this
 *	17-feb-94   cleanup; return double, not real
 *	19-may-95   added optional (-DETIME) etime/dtime
 *	18-apr-96   debuglevel of cputool at 4 now.
 *	31-mar-01   now using CLK_TCK instead of HZ (POSIX ?)
 *      21-nov-03   adding clock() for testing, using sysconf to get clk_tck;
 *      15-mar-05   fallback to CLOCKS_PER_SEC for g++
 *       7-nov-18   CLK_TCK erroneously at 1000000 now?
 *	
 */

#include <stdinc.h>
#include <limits.h>     /* solaris hides CLK_TCK here */
#include <time.h>	/* linux hides CLK_TCK here */
#include <sys/times.h>
#include <unistd.h>     /* for sysconf() */

#define TRY_CLOCK
// test  LINUX doesn't do CLK_TCK correct anymore?
#define CLK_TCK 100

/*  CLK_TCK is typically 100, in g++ CLK_TCK isn't known, uses CLOCK_PER_SEC, else fail */

#ifndef CLK_TCK
# ifndef CLOCKS_PER_SEC
#  error CLK_TCK as well as CLOCKS_PER_SEC definitions missing
# else
#  define CLK_TCK CLOCKS_PER_SEC
# endif
#endif

static long clk_tck = 0;

double cputime2(int mode)
{
    struct tms buffer;
    clock_t pt = 0;

    if (clk_tck == 0) {
      clk_tck = sysconf(_SC_CLK_TCK);
      dprintf(4,"cputime: initialized clk_tck=%ld [CLK_TCK=%ld]\n",clk_tck,CLK_TCK);
    }

    if (times(&buffer) == -1) error("times(2) call failed");
#ifdef TRY_CLOCK
    pt = clock();
#endif
    dprintf(4,"cputool: times: usr=%ld sys=%ld clock=%ld\n",buffer.tms_utime, buffer.tms_stime,pt);

    if (mode == 0)
        return buffer.tms_utime / ((double)CLK_TCK * 60.0);       /* return minutes */
    else if (mode == 1)
        return buffer.tms_stime / ((double)CLK_TCK * 60.0);       /* return minutes */
    else if (mode == 2)
        return pt/( (double)CLOCKS_PER_SEC*60.0);
}

double cputime()
{
  return cputime2(0);
}


/* Now some routines that some fortran users like:
 * ETIME returns the elapsed runtime in seconds for the calling process.
 * DTIME returns the elapsed time since the last call to dtime,
 *       or the start of execution on the first call.
 * The returned argument tarray contains user time in the first
 * element  and  system  time  in  the second element.  Elapsed
 * time, the returned value, is the  sum  of  user  and  system
 * time.
 */

#if defined(ETIME)
float etime(float *tarr)
{
    struct tms buffer;

    if (times(&buffer) == -1) error("times(2) call failed");
    tarr[0] = buffer.tms_utime / (1.0*CLK_TCK);
    tarr[1] = buffer.tms_stime / (1.0*CLK_TCK);
    return tarr[0] + tarr[1];
}

float dtime(float *tarr)
{
    struct tms buffer;
    float tmp;
    static float tusr=-1;
    static float tsys=-1;

    if (times(&buffer) == -1) error("times(2) call failed");
    tarr[0] = buffer.tms_utime / (1.0*CLK_TCK);
    tarr[1] = buffer.tms_stime / (1.0*CLK_TCK);
    if (tusr<0) {
        tusr = tarr[0];
        tsys = tarr[1];
    } else {
        tmp = tarr[0];  tarr[0] -= tusr;    tusr = tmp;
        tmp = tarr[1];  tarr[1] -= tsys;    tsys = tmp;
    }
    return tarr[0] + tarr[1];
}
#endif

/*===========================================================================*/

#if defined(TESTBED)

main(int ac,char *av[])
{
    int i,j,k, imax=10, jmax=10, kmax=10;
    double x=1,y=2,z=3;
    
    if (ac>1) imax=atoi(av[1]);
    if (ac>2) jmax=atoi(av[2]);
    if (ac>3) kmax=atoi(av[3]);
    
    for(i=0; i<imax;i++) {
       y += sqrt(z);
       for(j=0; j<jmax;j++) {
          z += sqrt(x);
          for(k=0; k<kmax;k++) {
             x += sqrt(y);
          }
       }
    }
    printf("CPU time = %f min for %d*%d*%d -> %g CLK_TCK=%d\n",
	cputime(),imax,jmax,kmax,x+y+z,CLK_TCK);
}

#endif


/*===========================================================================*/

#if defined(TOOLBOX)

/*-----------------------------------------------------------------------------
 *  main  --  driver to test  cputime() . It also gives useful information
 *            about the floating point speed of the machine at hand.
 *-----------------------------------------------------------------------------
 */

#include <getparam.h>

string defv[] =  {
   "count=4000\n       KiloLoop count (in 1000)",
   "mode=f\n           Float or Int test",
   "VERSION=1.3\n      31-mar-2001 pjt",
   NULL,
};

string usage="cpu benchmark of sorts";

nemo_main()
{
    int cnt = getiparam("count") * 1000;
    string mode = getparam("mode");

    if (*mode == 'f')
        float_test(cnt);
    else if (*mode == 'i')
        int_test(cnt);
    else
        error("Invalid mode: choose f or i");
}

float_test(int cnt)
{
    int  i;
    real a, b, c, d, e, f, g;
    real  t0, t1, t2, t3, t4, t5;
    real  addspeed, subspeed, mulspeed, divspeed, sqrtspeed;

    
    a = 0.0000001;
    b = 1.0000001;
    c = d = e = f = g = 1.0;

    t0 = cputime();

    i = cnt;
    while (i--)
        {
	c += a;
	d *= b;
	e -= a;
	f /= b;
	}

    t1 = cputime();

    printf("To calculate %d kilo-floating point operations (+,-,*,/)\n",4*cnt);
    printf("  took about %.3f seconds.\n", (t1 - t0)*60.0);
    printf("This implies a rough average floating point speed of ");
    printf("%.0f Mflops.\n", 4*cnt/1000000.0/(60.0 *(t1 - t0)));

    t0 = cputime();

    i = cnt;
    while (i--) c += a;
	
    t1 = cputime();

    i = cnt;
    while (i--) d *= b;

    t2 = cputime();

    i = cnt;
    while (i--) e -= a;
	
    t3 = cputime();

    i = cnt;
    while (i--) f /= b;
	
    t4 = cputime();

    i = cnt;
    while (i--) {
	d *= b;                     /* added to fool the optimizing compiler */
	g = sqrt(d);
    }	

    t5 = cputime();

    dummy(a, b, c, d, e, f, g);     /* we have to do something with these    */
                                    /* variables in order to fool the        */
                                    /* optimizing compilers on some machines */

    addspeed = cnt / (t1 - t0);
    mulspeed = cnt / (t2 - t1);
    subspeed = cnt / (t3 - t2);
    divspeed = cnt / (t4 - t3);
    sqrtspeed = cnt / ((t5 - t4) - (t2 - t1));

    printf("\n");
    printf("The following ratios were measured in floating point speed:\n  ");
    printf("add : sub. : mult. : div. : sqrt = ");
    printf("%.2f : %.2f : %.0f : %.2f : %.2f\n", addspeed/mulspeed,
	   subspeed/mulspeed, 1.0, divspeed/mulspeed, sqrtspeed/mulspeed);
}

int_test(int cnt)
{
    int  i;
    int a, b, c, d, e, f, g;
    real  t0, t1, t2, t3, t4, t5;
    real  addspeed, subspeed, mulspeed, divspeed, sqrtspeed;

    
    a = 1;
    b = 3;
    c = d = e = f = g = 1;

    t0 = cputime();

    i = cnt;
    while (i--)
        {
	c += a;
	d *= b;
	e -= a;
	f /= b;
	}

    t1 = cputime();

    printf("To calculate %d int operations (+,-,*,/)\n",4*cnt);
    printf("  took about %.3f minutes.\n", t1 - t0);
    printf("This implies a rough average int speed of ");
    printf("%.0f kips.\n", 4*cnt/1000.0/(60.0 *(t1 - t0)));

    t0 = cputime();

    i = cnt;
    while (i--)
	c += a;

    t1 = cputime();

    i = cnt;
    while (i--)
	d *= b;

    t2 = cputime();

    i = cnt;
    while (i--)
	e -= a;

    t3 = cputime();

    i = cnt;
    while (i--)
	f /= b;

    t4 = cputime();

    t5 = cputime();

    dummy(a, b, c, d, e, f, g);     /* we have to do something with these    */
                                    /* variables in order to fool the        */
                                    /* optimizing compilers on some machines */

    addspeed = cnt / (t1 - t0);
    mulspeed = cnt / (t2 - t1);
    subspeed = cnt / (t3 - t2);
    divspeed = cnt / (t4 - t3);
    sqrtspeed = cnt / ((t5 - t4) - (t2 - t1));

    printf("\n");
    printf("The following ratios were measured in floating point speed:\n  ");
    printf("add : sub. : mult. : div. : sqrt = ");
    printf("%.2f : %.2f : %.0f : %.2f : %.2f\n", addspeed/mulspeed,
	   subspeed/mulspeed, 1.0, divspeed/mulspeed, sqrtspeed/mulspeed);
}

dummy(real a, real b, real c, real d, real e, real f, real g) {}


#endif

