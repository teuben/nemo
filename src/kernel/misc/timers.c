/*
 * timers: collection of routines to timers
 *
 */

#include <stdinc.h>

/* 
 * readTSC:   reads the Time Stamp Counter of an i386 processor
 *            and returns it in a 64 bit 'long'
 */

#if defined(__GNUC__) && defined(i386)
long long readTSC(void) {
  /* this assumes long long is 64 bits. unsigned in 32 */
  union {
    long long complete;
    unsigned part[2];
  } ticks;
  __asm__ ( "rdtsc; mov %%eax,%0;mov %%edx,%1"
	      : "=mr" (ticks.part[0]),
	        "=mr" (ticks.part[1])
	      : /* No inputs */
	      : "eax", "edx");
  return ticks.complete;
}
#else
long long readTSC(void) {
  return 0;
}
#endif

static long long *timers = NULL;
static int maxtimers = 0;

void init_timers(int n)
{
  timers = (long long *) allocate(n*sizeof(long long));
  maxtimers = n;
  if (n>1) {
    stamp_timers(0);    stamp_timers(1);
    dprintf(1,"init_timers: overhead of called stamp_timers=%d ticks\n",diff_timers(0,1));
    stamp_timers(0);    stamp_timers(1);
    dprintf(1,"init_timers: overhead of called stamp_timers=%d ticks\n",diff_timers(0,1));
    stamp_timers(0);    stamp_timers(1);
    dprintf(1,"init_timers: overhead of called stamp_timers=%d ticks\n",diff_timers(0,1));
  }
}

void stamp_timers(int i)
{
  if (maxtimers==0) error("init_timers not called");
  if (i >= maxtimers) error("i=%d maxtimers=%d",i,maxtimers);
  timers[i] = readTSC();
}

int diff_timers(int i, int j)
{
  if (maxtimers==0) error("init_timers not called");
  if (i >= maxtimers) error("i=%d maxtimers=%d",i,maxtimers);
  if (j >= maxtimers) error("j=%d maxtimers=%d",j,maxtimers);
  
  return (int) (timers[j]-timers[i]);
}
