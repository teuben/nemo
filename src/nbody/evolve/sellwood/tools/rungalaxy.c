/*
 *  preprocessor to run (a series of) GALAXY:
 *
 *  The first argument must be a (non-existent) directory name, in
 *  which a new 'galaxy.dat' file will be written, and in which
 *  galaxy will be run. The 'galaxy' program is allowed to be in the
 *  unix PATH

 *	24-jun-1997		written
 *	18-jul-2001		default is galaxy.exe
 *
 * Todo:  incoorporate some kind of snap2ini script in the code.
 */

#include <stdinc.h>
#include <getparam.h>

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

string defv[] = {
    "in=???\n        input snapshot (nemo format)",
    "out=???\n       directory in which GALAXY will be run",
    "scale=15.0\n    number of grid cells per unit length",
    "dt=0.05\n       integration time step",
    "dtout=0.5\n     time between particle outputs",
    "dtlog=0.1\n     time between integral checks",
    "tstop=1.0\n     end time",
    "grid=33,33,33\n  number of grid cells in (x,y,z,)",
    "exe=galaxy.exe\n name of GALAXY executable",
    "VERSION=1.2\n    18-jul-01 PJT",
    NULL,
};

string usage = "run galaxy";

#define MAXCARDS 256

static char *database[MAXCARDS];    /* full lines */
static char *names[MAXCARDS];       /* individual namelist */
static int ncards = 0;

static char *haskey(char *, char *);

int nemo_main()
{
    int i, n, ngrid[3];
    real scale, dt, dtout, dtlog, tstop;
    string exefile = getparam("exe");
    string rundir = getparam("out");
    stream datstr;
    char fullname[256];

    n = nemoinpi(getparam("grid"),ngrid,3);
    if (n>0 && n<=3) {
        for (i=n; i<3; i++)
            ngrid[i] = ngrid[i-1];
    } else
        error("%d Syntax error: %s (need 1,2 or 3 integers)",
                n,getparam("grid"));
    scale = getdparam("scale");
    dt = getdparam("dt");
    dtout = getdparam("dtout");
    dtlog = getdparam("dtlog");
    tstop = getdparam("tstop");

    make_rundir(rundir);


    sprintf(fullname,"%s/%s",rundir,"galaxy.dat");
    datstr = stropen(fullname,"w");    
    fprintf(datstr,"%d %d %d\n",ngrid[0],ngrid[1],ngrid[2]);
    fprintf(datstr,"%g\n",scale);
    fprintf(datstr,"%g\n",dt);
    fprintf(datstr,"%g\n",dtout);
    fprintf(datstr,"%g\n",dtlog);
    fprintf(datstr,"%g\n",tstop);
    strclose(datstr);


    sprintf(fullname,"%s/%s",rundir,"galaxy.ini");
    datstr = stropen(fullname,"w");
    strclose(datstr);

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
