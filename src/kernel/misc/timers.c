/*
 * timers: collection of routines to timers
 *
 *         jul-2005     gleaned from some magazine?
 *      18-nov-2012     certified it works on x86_64 as well
 *
 */

#include <stdinc.h>
#include <timers.h>

/* 
 * readTSC:   reads the Time Stamp Counter of an intel processor
 *            and returns it in a 64 bit 'long'
 */

#if defined(__GNUC__) && (defined(i386) || defined(__x86_64))
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
  if (timers) {
    dprintf(2,"init_timers: free up old maxtimers %d \n",maxtimers);
    free(timers);
  }
  if (sizeof(long long) != 8) warning("timers::long long is not 8 bytes");
  if (sizeof(unsigned)  != 4) warning("timers::unsigned  is not 4 bytes");

  timers = (long long *) allocate(n*sizeof(long long));
  maxtimers = n;
  if (n>1) {
    stamp_timers(0);    stamp_timers(1);
    dprintf(1,"init_timers: overhead of called stamp_timers=%ld ticks\n",diff_timers(0,1));
    stamp_timers(0);    stamp_timers(1);
    dprintf(1,"init_timers: overhead of called stamp_timers=%ld ticks\n",diff_timers(0,1));
    stamp_timers(0);    stamp_timers(1);
    dprintf(1,"init_timers: overhead of called stamp_timers=%ld ticks\n",diff_timers(0,1));
  }
}

void stamp_timers(int i)
{
  if (maxtimers==0) error("init_timers not called");
  if (i >= maxtimers) error("i=%d maxtimers=%d",i,maxtimers);
  timers[i] = readTSC();
}

long long diff_timers(int i, int j)
{
  if (maxtimers==0) error("init_timers not called");
  if (i >= maxtimers) error("i=%d maxtimers=%d",i,maxtimers);
  if (j >= maxtimers) error("j=%d maxtimers=%d",j,maxtimers);
  
  return timers[j]-timers[i];
}



#ifdef TESTBED


nemo_main()
{
  int i, n=10;
    init_timers(n+1);
    for (i=0; i<n; i++)
	stamp_timers(i);
    stamp_timers(n);
    for (i=0; i<n; i++)
	printf("Method-1: %Ld\n",diff_timers(i,i+1));

    stamp_timers(0);
    for (i=0; i<n; i++) {
	stamp_timers(i+1);
	printf("Method-2: %Ld\n",diff_timers(i,i+1));
    }
}
#endif
