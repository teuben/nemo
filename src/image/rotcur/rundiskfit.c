/*
 *  preprocessor to run 'diskfit' 2012 version  :  rundiskfit, vs. rundiskfits12
 *
 *	13-sep-2012    2nd initial version, but very sloppy still. Only inp= is allowed
 * 
 *  @TODO: 
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

string defv[] = {
  "inp=???\n             Input parameter file",

  /* anything below this not used yet */

  "in=\n              input velocity field",
  "out=\n             output model field",
  "par=\n             output parameter file",
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
  
  /* anything above this not used yet */

  "regrid=f\n            Rasterize ascii tabular data into an image?",
  "exe=diskfit\n         name of VELFITSS07 executable in your $PATH",
  "VERSION=0.1\n         13-sep-2012 PJT",
  NULL,
};

string usage = "frontend to run Spekkens & Sellwood's 2012 diskfit code";

string cvsid="$Id$";

int run_program(string);

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
  stream datstr;
  string inpfile = getparam("inp");
  
  string exefile = getparam("exe");
  char dname[256];
  char command[256];
  bool  Qraster = getbparam("regrid");
  FILE  *fp;


  if (hasvalue("inp"))  {
    warning("sloppy version that only accepts inp= and regrid=");
    fp = popen(exefile,"w");

    sprintf(dname,"'%s'\n",inpfile);
    fprintf(fp,dname);
    if (Qraster)
      fprintf(fp,"y");
    else
      fprintf(fp,"n");
    fclose(fp);

    stop(0);
  } else
    error("sloppy version that only accepts inp= and regrid=");



  /* anything below this is still the old SS07 V1.0 setup */

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

  if (run_program(exefile))
    error("Problem executing %s",exefile);
}


int run_program(string exe)
{
  FILE *fp;
  fp = popen(exe,"w");
  fprintf(fp,"'velfit.inp'");

  fclose(fp);
  return 0;
}
