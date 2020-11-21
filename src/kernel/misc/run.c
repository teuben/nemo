/*
 * run: tools to run external non-NEMO programs
 *
 *      17-sep-2013     finally written
 *      26-dec-2017     add include's
 *
 */

#include <stdinc.h>
#include <run.h>

/* mkdir */
#include <sys/stat.h>
#include <sys/types.h>
/* chdir */
#include <unistd.h>

FILE *popen(const char *command, const char *type);  //C99


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

/*
 *   echo file | exe
 */

int run_popen1(string exe, string file)
{
  FILE *fp = popen(exe,"w");
  if (fp==NULL) return 1;
  fprintf(fp,"%s",file);
  fclose(fp);
  return 0;
}

#ifdef TESTBED


void nemo_main(void)
{
  dprintf(0,"wc on hello:");
  run_popen1("wc","hello");
}
#endif
