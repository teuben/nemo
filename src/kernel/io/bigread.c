/*
 *  BIGREAD:    read file in buffers and measure read speed for each buffer
 *
 *    11-jul-2020    V0.1    - drafted    Peter Teuben
 */

#include <nemo.h>
#include <timers.h>

string defv[] = {
    "in=???\n        Input file",
    "bufsize=1\n     Buffer size, in (1024 sytle) GB",
    "mode=1\n        cputime2(mode)    0:usr 1:sys 2:clock",
    "VERSION=0.2\n   11-jul-2020 PJT",
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
  double t0, t1;
  int i,n;
  char *buffer;

  printf("BufSize: %ld bytes\n", bufsize);
  buffer = (char *) malloc(bufsize);
  if (buffer == NULL) error("Cannot allocate %ld",bufsize);

  init_timers(MAX_TIMERS+1);  
  t0 = cputime2(mode) * 60.0;
  dprintf(0,"Start-time %g sec, mode=%d\n",t0,mode); 
  for(n=0;;n++) {
    nread = fread(buffer, sizeof(char), bufsize, instr);
    stamp_timers(n);
    dprintf(1,"%d %ld\n",n,nread);
    if (nread <= 0) break;
    t1 = cputime2(mode)*60.0;
    dprintf(0,"%d %g\n", n, nread/1024.0/1024.0/(t1-t0));
    t0 = t1;
  }
  n=n-1;
  for (i=0; i<n; i++)
    dprintf(1,"Method-2: %d %lld\n",i,diff_timers(i,i+1));

  t1 = cputime2(mode)*60.0;
  dprintf(0,"End-time %g sec\n",t1);
}
