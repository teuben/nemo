/*
 *  preprocessor to run Spekkens and Sellwood 2007 'velfit' program
 *
 *	27-oct-2008    initial version
 *      17-sep-2013    use new run interface..
 *      11-aug-2014    V1.2 : allow inp= to mirror how rundiskfit works 
 * 
 *  @TODO: 
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <run.h>

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

string defv[] = {
  "in=???\n              input velocity field",
  "out=???\n             output model field",
  "par=???\n             output parameter file",
  "center=0,0\n          center of vel field",
  "vsys=0\n              systemic velocity",
  "pa=0\n                position angle of disk (E of N)",
  "inc=0\n               inclination angle of disk ",
  "rmsism=0\n            ISM turbulence",
  "flags=1,1,1,0,0,0\n   Flags which to be fitted/done (geom, center, vsys, radial, bar, error)",
  "m=2\n                 harmonic order of bar perturbation (1 or 2)",
  "seed=-50\n            Random Seed",
  "bootstrap=200\n       Number of bootstrap samples if error to be determined",
  "j=1.0\n               correllation length",
  "rcirc=50\n            Radius beyond which no non-circular rmotions fitted",
  "radii=\n              Ring Radii at which velocity field components are extracted",

  /* inp= is optional, it's the original way to run velfitss07 */

  "inp=\n                Input parameter file (optional way to run)",
  "exe=velfitss07\n      name of VELFITSS07 executable in your $PATH",
  "par=velfit.inc\n      parameter file passed for VELFITSS07",
  "VERSION=1.2\n         11-aug-2014 PJT",
  NULL,
};

string usage = "frontend to run Spekkens & Sellwood's 2007 velfit code";

string cvsid="$Id$";

#define MAXRAD   100
#define MAXFLAGS   6

int nemo_main()
{
  int i, j, n, m, nrad, seed, flags[6], nflags, ncenter, bootstrap, nradii;
  real scale, dt, dtout, dtlog, tstop, tsnap, mass;
  real center[2], pa, inc, vsys, rmsism, rcirc, radii[MAXRAD];
  string infile = getparam("in");
  string outfile = getparam("out");
  string parfile = getparam("par");
  string inpfile = getparam("inp");
  stream datstr;
  
  string exefile = getparam("exe");
  char dname[256];
  char command[256];
  FILE  *fp;

  if (hasvalue("inp"))  {
    warning("sloppy version that only accepts inp=");
    fp = popen(exefile,"w");
    sprintf(dname,"'%s'\n",inpfile);
    fprintf(fp,dname);
    fclose(fp);
    stop(0);
  } 

#if 0
  seed = init_xrandom(getparam("seed"));
#else
  seed = getiparam("seed");
#endif
  
  nflags  = nemoinpi(getparam("flags"),flags,MAXFLAGS);
  ncenter = nemoinpr(getparam("center"),center,2);
  nradii  = nemoinpr(getparam("radii"),radii,MAXRAD);
  j = getiparam("j");


  /* sanity checks before writing out input file */
  if (nflags != MAXFLAGS) error("flags= needs %d values",MAXFLAGS);
  if (ncenter != 2) error("center= needs 2 values");

  if (seed >=0 ) error("seed must be negative for SS07");
  if (strlen(infile)  > 100) warning("infile too long");
  if (strlen(outfile) > 100) warning("outfile too long");
  if (strlen(parfile) > 100) warning("parfile too long");


  sprintf(dname,"%s","velfit.inp");
  datstr = stropen(dname,"w!");    
  fprintf(datstr,"# input file for velfitss07\n");
  fprintf(datstr,"'%s'\n",infile);                 /* invfile */
  fprintf(datstr,"'%s'\n",parfile);                /* outpfile */
  fprintf(datstr,"'%s'\n",outfile);                /* outmfile */
  fprintf(datstr,"%g\n",center[0]);                /* xcen*/
  fprintf(datstr,"%g\n",center[1]);                /* ycen*/
  fprintf(datstr,"%g\n",getrparam("vsys"));        /* vsys */
  fprintf(datstr,"%g\n",getrparam("pa"));          /* pa */
  fprintf(datstr,"%g\n",getrparam("inc"));         /* incl */
  fprintf(datstr,"%g\n",getrparam("rmsism"));      /* eISM */
  for (i=0; i<MAXFLAGS; i++)
    fprintf(datstr,"%s ", flags[i] ? "T" : "F");   /* disk,centre,systemic,radial,bisymm,uncert */
  fprintf(datstr,"\n");
  fprintf(datstr,"%d\n",getiparam("m"));           /* order */
  fprintf(datstr,"%d\n",seed);
  fprintf(datstr,"%d\n",getiparam("bootstrap"));   /* nunc */
  fprintf(datstr,"%g\n",getrparam("j"));           /* junc */
  fprintf(datstr,"%g\n",getrparam("rcirc"));       /* maxr */
  for (i=0; i<nradii; i++)
    fprintf(datstr,"%g\n", radii[i]);              /* sma(i) */
  

  strclose(datstr);

  if (run_popen1(exefile,"'velfit.inp'"))
    error("Problem executing %s",exefile);
}

/*

  Input Directive File (IDF)     qs=quoted string  b=bool  r=real

  #  comment line, doesn't count
  :comment
  qs:vinfile
  qs:poutfile
  qs:moutfile
  r:xcen 
  r:ycen
  r:vsys
  r:pa 
  r:inc
  r:rmsism
  b:fitgeo b:fitcen b:fitvsys b:fitm0 b:fitm12 b:fitunc
  i:m
  i:seed
  i:bootstrap
  r:j
  r:rcirc
  # now the array starts with K sample radii  r1...rK
  r:r[]

 */
