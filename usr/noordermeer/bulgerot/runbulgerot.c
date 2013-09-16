/*
 *
c     This program calculates the rotation curve for an oblate spheroidal 
c     bulge. It assumes a projected surface density that follows a Sersic 
c     profile. The central surface magnitude M_0, characteristic radius R_0 and
c     Sersic index n_sers must be given. Also, the inclination and intrinsic 
c     axis ratio of the bulge must be given. Finally, the distance and the 
c     Galactic foreground extinction of the galaxy is needed.
c
c     Copyright (C) 2007, Edo Noordermeer
c     E-mail: edo.noordermeer@gmail.com
c
c     If you use this program for a scientific publication, I would appreciate 
c     a reference to the paper 'The rotation curves of flattened Sersic bulges
c     (Noordermeer 2008)'.
 */


#include <nemo.h>

string defv[] = {
  "outdir=???\n               output run directory",
  "mu0=15\n                   apparent central R-band surface magnitude M_0",
  "r0=1\n                     characteristic radius R_0 in arcseconds",
  "n=1\n                      Sersic index",
  "inc=0\n                    Inclination (0=face on)",
  "q=1\n                      Axis ratio b/a (1=round)",
  "dist=10\n                  Distance in Mpc",
  "ar=0\n                     galactic foreground extinction",
  "r3=0,0.01,250\n            rstart,rstep,nradii for rotation curve, in kpc",
  "galaxy=sersic\n            Identification string",
  "exe=bulgerot\n             Name of executable (needs to be in $PATH)",
  "VERSION=1.2\n              15-sep-2013 PJT",
  NULL,
};


string usage="rotation curve of an oblate spheroidal sersic bulge";

string cvsid="$Id$";
 
int run_program(string);


nemo_main()
{
  string exefile = getparam("exe");
  string rundir  = getparam("outdir");
  string infile  = "bulgerot.in";
  string datfile = "bulgerot.dat";
  string logfile = "bulgerot.log";
  real r[3];
  int n;
  char dname[256], command[256];
  stream datstr;
  
  n = nemoinpr(getparam("r3"),r,3);
  if (n!=3) error("parsing error %s: r3 needs rmin,rstep,nradii",getparam("r3"));
  n = r[2];

  make_rundir(rundir);
  sprintf(dname,"%s/%s",rundir,infile);
  datstr = stropen(dname,"w");
  fprintf(datstr,"%s\n",datfile);
  fprintf(datstr,"%s\n",getparam("galaxy"));
  fprintf(datstr,"%g %g %g\n",
	  getrparam("mu0"),getrparam("r0"),getrparam("n"));
  fprintf(datstr,"%g %g %g %g\n",
	  getrparam("inc"),getrparam("q"),1e6*getrparam("dist"),getrparam("ar"));
  fprintf(datstr,"%g %g %d\n",r[0],r[1],n);
  strclose(datstr);

  sprintf(dname,"%s/%s",rundir,"runme");
  datstr = stropen(dname,"w");
  sprintf(command,"#! /bin/sh\n%s < %s >%s 2>&1\n",exefile,infile,logfile);
  dprintf(1,"\n%s",command);
  fprintf(datstr,"%s",command);
  strclose(datstr);

  goto_rundir(rundir);
  if (run_program("sh runme"))
    error("Problem executing runme in %s",rundir);
}

void goto_rundir(string name)
{
    if (chdir(name))
        error("Cannot change directory to %s",name);
    return 0;
}

void make_rundir(string name)
{
    if (mkdir(name, 0755))
        error("Run directory %s already exists",name);
}

int run_program(string cmd)
{
  int retval;
  retval = system(cmd);
  dprintf(1,"%s: returning %d\n",cmd,retval);
  return retval;
}
