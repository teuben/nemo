/*
 *  preprocessor to run (a series of) CGS
 *
 *  The first argument must be a (non-existent) directory name, in
 *  which a new 'galaxy.dat' file will be written, and in which
 *  galaxy will be run. The 'galaxy' program is assumed to be in the
 *  PATH
 *
 *	3-nov-2005		written
 */

#include <stdinc.h>
#include <getparam.h>

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

string defv[] = {
    "out=???\n       directory in which GALAXY will be run",
    "nbody=40000\n   nbody",
    "nrad=80\n       nradii",
    "maxstep=1000\n  max",
    "dt=0.0025\n     step",
    "tstart=0.0\n    start",
    "tstop=4.0\n     stop",
    "dtmin=0.001\n   min",
    "dtmax=0.01\n    max",
    "flag=1\n        flag",
    "mass=1.0\n      mass",
    "freqcmss=2000\n freq",
    "freqdiag=50\n   freq",
    "exe=CGS.exe\n   name of CGS executable",
    "in=\n           optional input snapshot (nemo format)",
    "VERSION=0.1\n   3-nov-05 PJT",
    NULL,
};

string usage = "run CGS";

#define MAXCARDS 256

static char *database[MAXCARDS];    /* full lines */
static char *names[MAXCARDS];       /* individual namelist */
static int ncards = 0;

int nemo_main()
{
  int i, nbody;
  real scale, dt, dtout, dtlog, tstop;
  string exefile = getparam("exe");
  string rundir = getparam("out");
  stream datstr;
  char fullname[256];
  char command[256];
  
  nbody = getiparam("nbody");

  make_rundir(rundir);


  sprintf(fullname,"%s/%s",rundir,"PARAMETER.DAT");
  datstr = stropen(fullname,"w");    
  fprintf(datstr,"%d\n",getiparam("nrad"));
  fprintf(datstr,"%d\n",nbody);
  fprintf(datstr,"%d\n",getiparam("maxstep"));
  fprintf(datstr,"%d\n",getiparam("freqcmss"));
  fprintf(datstr,"%d\n",getiparam("freqdiag"));
  fprintf(datstr,"%g\n",getdparam("dt"));
  fprintf(datstr,"%g\n",getdparam("tstart"));
  fprintf(datstr,"%g\n",getdparam("tstop"));
  fprintf(datstr,"%g\n",getdparam("mass"));
  fprintf(datstr,"%d\n",getiparam("flag"));
  fprintf(datstr,"%g\n",getdparam("dtmax"));
  fprintf(datstr,"%g\n",getdparam("dtmin"));

  strclose(datstr);

#if 0  
#if 0
    sprintf(fullname,"%s/%s",rundir,"galaxy.ini");
    datstr = stropen(fullname,"w");
    strclose(datstr);
#else
    sprintf(command,"snap2ini %s > %s/%s",
	    getparam("in"),rundir,"galaxy.ini");
    dprintf(1,"COMMAND: %s\n",command);
    system(command);
#endif
#endif
    goto_rundir(rundir);

    run_program(exefile);
}


goto_rundir(string name)
{
    if (chdir(name))
        error("Cannot change directory to %s",name);
}

make_rundir(string name)
{
    if (mkdir(name, 0755))
        error("Run directory %s already exists",name);
}

write_namelist(char *name)
{
    int i;
    FILE *fp;

    fp = fopen(name,"w");
    if (fp==NULL) {
        fprintf(stderr,"File %s could not be written\n",name);
        exit(1);
    }
    for (i=0; i<ncards; i++)
        fprintf(fp,"%s",database[i]);
    fclose(fp);
}



run_program(string exe)
{
#if 1
    system(exe);
#else
    if (execlp(exe,NULL)) {
        fprintf(stderr,"Problem executing %s\n",exe);
        exit(1);
    }
#endif
}
