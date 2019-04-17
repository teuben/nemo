/*
 *  preprocessor to run (a series of) CGS
 *
 *
 *	3-nov-2005     written
 *     12-dec-2005     added writing pos/vel files with freqout=
 *     16-mar-2006     0.4 update version 
 *     22-mar-2006     0.5 also read forces and potentials
 *     22-may-2006     1.0 readied for more formal release within NEMO; no more freq* but ndt*
 *     17-sep-2013     1.1 new run interface
 */

#include <stdinc.h>
#include <getparam.h>
#include <run.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

string defv[] = {
    "outdir=???\n    directory in which GALAXY will be run",
    "nbody=40000\n   nbody, if no in= given",
    "nrad=80\n       nradii in code",
    "maxstep=1000\n  max number of integration steps",
    "dt=0.0025\n     initial integration step",
    "tstart=0.0\n    start time",
    "tstop=4.0\n     stop time",
    "dtmin=0.001\n   min integration time",
    "dtmax=0.01\n    max integration time",
    "flag=1\n        data creation flag",
    "mass=1.0\n      total mass of system",
    "virial=1.0\n    initial virial ratio if flag=2",
    "ndtcmss=2000\n  number of steps between center of mass adjustments",
    "ndtdiag=50\n    number of steps between diagnostics",
    "ndtout=10\n     number of steps between snapshot output",
    "in=\n           optional input snapshot (nemo format)",
    "nemo=t\n        convert data to NEMO and cleanup ASCII",
    "options=\n      Optional output:  phi(potential), acc (forces) [not implemented yet]",
    "exe=CGS.exe\n   name of CGS executable",
    "VERSION=2.1a\n  9-apr-2019 PJT",
    NULL,
};

string usage = "NEMO frontend to run CGS";

string cvsid="$Id$";



void nemo_main()
{
  int i, nbody, flag;
  real scale, dt, dtout, dtlog, tstop; 
  real virial = getrparam("virial");
  bool Qnemo = getbparam("nemo");
  bool Qpot = TRUE;
  bool Qacc = TRUE;
  string exefile = getparam("exe");
  string rundir = getparam("outdir");
  string infile;
  stream datstr;
  char fullname[256], command[256], line[256];
  string rmcmd = "rm -f";


  run_mkdir(rundir);

  if (hasvalue("in")) {
    infile = getparam("in");
    flag = 3;       /* discard the flag= parameter, fix it here */
    sprintf(command,"snapprint %s x,y,z    > %s/%s",infile,rundir,"initPOS.dat");
    dprintf(0,">>> %s\n",command);
    system(command);
    sprintf(command,"snapprint %s vx,vy,vz > %s/%s",infile,rundir,"initVEL.dat");
    dprintf(0,">>> %s\n",command);
    system(command);
    sprintf(command,"cat %s/%s | wc -l",rundir,"initVEL.dat");
    datstr = popen(command,"r");
    if (fgets(line,256,datstr) == 0) {
      warning("bad read finding nbody, trying keyword");
      nbody = getiparam("nbody");
    } else {
      line[strlen(line)-1] = 0;
      nbody = atoi(line);
    }
    dprintf(0,"line=%s\n",line);
    dprintf(0,"Using snapshot in=%s with nbody=%d\n",infile,nbody);
  } else {
    infile = 0;
    flag = getiparam("flag");
    nbody = getiparam("nbody");
  }

  if (flag == 2) {
    sprintf(fullname,"%s/%s",rundir,"initial_virial.dat");    
    datstr = stropen(fullname,"w");    
    fprintf(datstr,"%g\n",virial);
    strclose(datstr);
  }

  sprintf(fullname,"%s/%s",rundir,"PARAMETER.DAT");
  datstr = stropen(fullname,"w");    
  fprintf(datstr,"%d\n",getiparam("nrad"));
  fprintf(datstr,"%d\n",nbody);
  fprintf(datstr,"%d\n",getiparam("maxstep"));
  fprintf(datstr,"%d\n",getiparam("ndtcmss"));
  fprintf(datstr,"%d\n",getiparam("ndtdiag"));
  fprintf(datstr,"%d\n",getiparam("ndtout"));
  fprintf(datstr,"%g\n",getdparam("dt"));
  fprintf(datstr,"%g\n",getdparam("tstart"));
  fprintf(datstr,"%g\n",getdparam("tstop"));
  fprintf(datstr,"%g\n",getdparam("mass"));
  fprintf(datstr,"%d\n",flag);
  fprintf(datstr,"%g\n",getdparam("dtmax"));
  fprintf(datstr,"%g\n",getdparam("dtmin"));
  strclose(datstr);

  run_cd(rundir);

  sprintf(command,"%s > CGS.log",exefile);
  dprintf(0,">>> %s\n",command);
  run_sh(command);

  if (Qnemo) {
    /* only supporting Qpot and Qacc, pos and vel always written  */
    if (Qpot && Qacc)
      sprintf(command,"tabtos fort.90 snap.out nbody,time skip,pos,vel,acc,phi; %s fort.90",rmcmd);
    else if (Qpot)
      sprintf(command,"tabtos fort.90 snap.out nbody,time skip,pos,vel,skip,skip,skip,phi; %s fort.90",rmcmd);
    else if (Qacc)
      sprintf(command,"tabtos fort.90 snap.out nbody,time skip,pos,vel,acc; %s fort.90",rmcmd);
    else
      sprintf(command,"tabtos fort.90 snap.out nbody,time skip,pos,vel; %s fort.90",rmcmd);
    dprintf(0,">>> %s\n",command);
    system(command);
  }
  
}

