/*
 * PROGRESS:  write timely progress messages, and keep the user happy
 *
 *
 *      30-jun-04 V1.0  created
 *       1-nov-04  1.1  made it an 'int' function instead
 */

#include <stdinc.h>
#include <stdarg.h>
#include <unistd.h>


static int bypass = -1;
static double cpu0;

extern double cputime(void);

int progress(double dtime, string fmt, ...)
{
    va_list ap;
    char cp;
    double cpu1;

    if (bypass < 0) {                /* initialize a few things */
      dprintf(0,"Progress bar initialized; dtime=%g cpu0=%g\n",dtime,cpu0);
      if (isatty(fileno(stderr)))
	bypass = 0;
      else
	bypass = 1;
      cpu0 = cputime()*60.0;
    } 
    if (bypass) return 0;

    if (dtime > 0) {
      cpu1 = cputime()*60.0;
      dprintf(1,"check: %g %g\n",cpu0,cpu1);
      if (cpu1-cpu0 < dtime) return 0;
      cpu0 = cpu1;
    }
    if (fmt==0 || *fmt==0) return 1;

    va_start(ap, fmt);              /* ap starts with string 'fmt' */

#if defined(NEED_DOPRNT)
    _doprnt(fmt, ap, stderr);       /* print out on stderr */
#else
    vfprintf(stderr, fmt, ap);      /* print out on stderr */
#endif

    cp = fmt[strlen(fmt)-1];
    if (cp != '\n' && cp != '\r') 
      fprintf(stderr,"\r");         
    fflush(stderr);                 /* flush it NOW */
    va_end(ap);                     /* end varargs */
    return 1;
}

#ifdef TESTBED

string defv[]={
  "n=10\n         number of loops",
  "m=1\n          report every m",
  "sleep=0\n	  delay?",
  "cpu=0\n        cpu delay?",
  "compute=1\n    compute something",
  "int=f\n        check int value",
  "VERSION=1.1\n  1-nov-04 PJT",
  NULL,
};

string usage="testing ";

void do_compute(int n)
{
  double a,b,c;
  int i,j;

  a = 2.0;
  b = 3.0;
  c = 0.0;
  for (i=0; i<n; i++)
  for (j=0; j<n; j++)
    c += sqrt(a)*sqrt(a) + sqrt(b)*sqrt(b);
  dprintf(2,"c=%g\n",c);
}

nemo_main()
{
  int k, n0, n = getiparam("n");
  int m = getiparam("m");
  int isleep = getiparam("sleep");
  double cpu = getdparam("cpu");
  int ncomp = getiparam("compute");
  bool Qint = getbparam("int");

  n0 = n;
  dprintf(1,"Going to print:\n");
  while(n-- > 0) {
    if (n%m == 0) {
      if (ncomp) do_compute(ncomp);
      if (isleep) sleep(isleep);
      if (Qint) {
	k=progress(cpu,0);
	if (k) dprintf(0,"Done %d/%d\r",n,n0);
      } else
	progress(cpu,"Done %d/%d",n,n0);
    }
  }
  dprintf(1,"\n All done.\n");

}

#endif

