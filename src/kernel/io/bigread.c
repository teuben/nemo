/*
 *  BIGREAD:    read file in buffers and measure read speed for each buffer
 *
 *    11-jul-2020    V0.1    - drafted    Peter Teuben
 */

#include <nemo.h>
#include <timers.h>
#include <moment.h>

string defv[] = {
    "in=???\n        Input file",
    "bufsize=1\n     Buffer size, in (1024 sytle) GB",
    "mode=1\n        cputime2(mode)    0:usr 1:sys 2:clock",
    "VERSION=0.3\n   11-jul-2020 PJT",
    NULL,
};

string usage="measured buffered read speed";

string cvsid="$Id:$";

#define MAX_TIMERS 10000

void nemo_main()
{
  size_t bufsize = (size_t) (getrparam("bufsize") * 1024 * 1024 * 1024);
  size_t nread;
  stream instr = stropen(getparam("in"),"r");
  int mode = getiparam("mode");
  double t0, t1, s;
  int i,n;
  char *buffer;
  Moment m;

  printf("BufSize: %ld bytes\n", bufsize);
  buffer = (char *) malloc(bufsize);
  if (buffer == NULL) error("Cannot allocate %ld",bufsize);

  ini_moment(&m,2,MAX_TIMERS+1);

  init_timers(MAX_TIMERS+1);  
  t0 = cputime2(mode) * 60.0;
  dprintf(0,"Start-time %g sec, mode=%d\n",t0,mode);
  
  for(n=0;;n++) {
    nread = fread(buffer, sizeof(char), bufsize, instr);
    stamp_timers(n);
    dprintf(2,"%d %ld\n",n,nread);
    if (nread <= 0) break;
    t1 = cputime2(mode)*60.0;
    s = nread/1024.0/1024.0/(t1-t0);
    printf("%d %g\n", n, s);
    accum_moment(&m, s, 1.0);
    t0 = t1;
  }
  n=n-1;
  for (i=0; i<n; i++)
    dprintf(1,"Timers: %d %lld\n",i,diff_timers(i,i+1));

  t1 = cputime2(mode)*60.0;
  dprintf(0,"End-time %g sec\n",t1);

  printf("Speed: n=%d min=%g max=%g mean=%g median=%g rms=%g (all)\n",
	 n_moment(&m),min_moment(&m),max_moment(&m),mean_moment(&m),
	 median_moment(&m),sigma_moment(&m));
  compute_robust_moment(&m);
  printf("Speed: n=%d min=%g max=%g mean=%g median=%g rms=%g (robust)\n",
	 n_robust_moment(&m),min_robust_moment(&m),max_robust_moment(&m),
	 mean_robust_moment(&m),median_robust_moment(&m),sigma_robust_moment(&m));
}
