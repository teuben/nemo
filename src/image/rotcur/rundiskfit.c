/*
 *  preprocessor to run 'diskfit' 2012 version  :  rundiskfit, vs. rundiskfits12
 *
 *	13-sep-2012    2nd initial version, but very sloppy still. Only inp= is allowed
 *      18-sep-2012    
 *      17-sep-2013    new run interface
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

  "mode=phot\n        Fitting mode: phot or vel",
  "in=\n              input data (intensities or velocities)",

  "out=\n             output model field",
  "par=\n             output parameter file",

  "sky=800,14.14,4.0\n   PHOT: image params:   sky, sky sig, gain",
  "region=1,1,256,256\n  FITS Region to fit:   xlow,ylow,xrange,yrange",
  "sample=10.0,35.0,0.4,2,0.5\n  FITS sampling: regrad,regpa,regeps,istepout,pixscale",

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
  
  /* simple method if inp= is given */

  "inp=\n                Input parameter file",
  "regrid=f\n            Rasterize ascii tabular data into an image?",
  "exe=diskfit\n         name of VELFITSS07 executable in your $PATH",
  "VERSION=0.2\n         11-aug-2014 PJT",
  NULL,
};

string usage = "frontend to run Spekkens & Sellwood's 2012 diskfit code";

string cvsid="$Id$";


#define MAXRAD   100
#define MAXFLAGS   6

int nemo_main()
{
  int i, j, n, m, nrad, seed, flags[6], nflags, ncenter, bootstrap, nradii;
  real scale, dt, dtout, dtlog, tstop, tsnap, mass;
  real center[2], pa, inc, vsys, rmsism, rcirc, radii[MAXRAD];
  real sky[3], sample[5];
  int region[4];
  string mode   = getparam("mode");
  string infile = getparam("in");
  string outfile = getparam("out");
  string parfile = getparam("par");
  stream datstr;
  string inpfile = getparam("inp");
  
  string exefile = getparam("exe");
  char dname[256];
  char command[256];
  bool  Qraster = getbparam("regrid");
  bool  Qphot;
  FILE  *fp;

  warning("NEMO interface to DISKFIT under development");

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
  } 

  if (*mode == 'p')        /* photometry */
    Qphot = TRUE;
  else if (*mode == 'v')   /* velocities / kinematics */
    Qphot = FALSE;
  else
    error("mode=phot or mode=vel are the only valid options");



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

  if (Qphot) {
    if (nemoinpr(getparam("sky"),sky,3) != 3) error("sky= needs 3 values");
    if (nemoinpi(getparam("region"),region,4) != 4) error("region= needs 4 values");
    if (nemoinpr(getparam("sample"),sample,5) != 5) error("sample= needs 5 values");
  } else {
    warning("mode=vel not yet supported");
  }


  sprintf(dname,"%s","diskfit.inp");
  datstr = stropen(dname,"w!");    
  fprintf(datstr,"# input file for diskfits\n");
  if (Qphot) {
    fprintf(datstr,"phot\n");
    fprintf(datstr,"N/A\n");
    fprintf(datstr,"'%s'\n",infile);               /* invfile */   
    fprintf(datstr,"%g %g %g\n",sky[0],sky[1],sky[2]);
    fprintf(datstr,"%d %d %d %d\n",region[0],region[1],region[2],region[3]);
    fprintf(datstr,"%g %g %g %d %g\n",sample[0],sample[1],sample[2],(int)sample[3],sample[4]);
  } else {
    fprintf(datstr,"vels\n");
    fprintf(datstr,"F F\n");    // only for text file ala velfitss07 now
    fprintf(datstr,"'%s'\n",infile);               /* invfile */   
    fprintf(datstr,"N/A\n");
    fprintf(datstr,"N/A\n");
    fprintf(datstr,"N/A\n");
  }
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

  if (run_popen1(exefile,"'diskfit.inp'"))
    error("Problem executing %s",exefile);
}

/*

  Input Directive File (IDF)

  #  line with comment
  s:mode              # this is for mode=phot
  b:fits b:velms      # only used if mode=='vel', else skipped
  qs:vinfile
  r:sky r:sky_sig r:gain
  i:xlow i:ylow i:xrange i:yrange
  r:regrad r:regpa r:regeps i:istepout r:pixscale
  qs:outfile
  b:fit_pa b:fit_eps b:fit_cen
  r:pa r:eps
  r:xcen r:ycen
  b:fit_bar b:fit_bar_pa b:fit_bar_eps r:bar_pa r:bar_eps
  # ignored in mode=phot 
  b:fit_bul b:fit_r_e r:r_e
  b:fit_sersic_n b:fit_bul_eps  r:sersic_n r:bul_eps
  r:seeing
  r:lambda_1  r:lambda_2
  b:toggle i:seed  i:nunc r_junc
  b:verbose
  r:bar_rmin r:bar_rmax
  r:ring[]
  
  #  line with comment
  s:mode              # this is for mode=vel
  b:fits b:velms      # only used if mode=='vel', else skipped
  qs:vinfile
  r:sky r:sky_sig r:gain      # ignored for mode=phot
  i:xlow i:ylow i:xrange i:yrange
  r:regrad r:regpa r:regeps i:istepout r:pixscale
  qs:outfile
  b:fit_pa b:fit_eps b:fit_cen
  r:pa r:eps
  r:xcen r:ycen

now it all switches

  b:fit_bar b:fit_bar_pa b:fit_bar_eps r:bar_pa r:bar_eps
  # ignored in mode=phot 
  b:fit_bul b:fit_r_e r:r_e
  b:fit_sersic_n b:fit_bul_eps  r:sersic_n r:bul_eps
  r:seeing
  r:lambda_1  r:lambda_2
  b:toggle i:seed  i:nunc r_junc
  b:verbose
  r:bar_rmin r:bar_rmax
  r:ring[]
  
  

 

 */
