/*
 * run: tools to run external non-NEMO programs
 *
 *      17-sep-2013     finally written
 *
 */

#include <stdinc.h>
#include <run.h>


int run_mkdir(string name)
{
  if (mkdir(name,0755)) {
    error("Run directory %s already exists",name);
    return 1;
  }
  return 0;
}

int run_cd(string name)
{
  if (chdir(name)) {
    error("Cannot change directory to %s",name);
    return 1;
  }
  return 0;
}

int run_sh(string cmd)
{
  int retval;
  retval = system(cmd);
  dprintf(1,"%s: returning %d\n",cmd,retval);
  return retval;
}



#ifdef TESTBED


nemo_main()
{
  warning("nothing here yet");
}
#endif
