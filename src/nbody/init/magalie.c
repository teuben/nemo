/*  magalie:   make a galaxy (MaGalie; Boily et al. prescription)
 *
 *  <dark-ages> Original fortran code               Boily et al.
 *  21-mar-04   interface created in rainy Strasbourg           PJT/CB
 *
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <history.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

string defv[] = {
  "out=???\n      output snapshot (a rundirectory $out.tmpdir is also created)",
 
  "ndisk=20000\n  Number of particles in disk  (use 0 to skip this component)",
  "nbulge=5000\n  Number of particles in bulge (use 0 to skip this component)",
  "nhalo=40000\n  Number of particles in halo  (use 0 to skip this component)",
  "seed=0\n       Random seed",
  "cleanup=f\n    cleanup run directory after use (not used yet)",
  "VERSION=1.1\n  23-mar-04 PJT",
  NULL,
};

string usage="Boily et al. composite bulge-disk-halo model";

void goto_rundir(string name);
void make_rundir(string name);
void run_program(string cmd);

void nemo_main(void)
{
  char rundir[256];
  stream datstr, histr;
  string out=getparam("out");
  int seed, nbulge, ndisk, nhalo;
  bool Qcleanup = getbparam("cleanup");

  warning("** program is in development **");
  
  seed =   init_xrandom(getparam("seed")); 
  
  nbulge = getiparam("nbulge");
  ndisk  = getiparam("ndisk");
  nhalo  = getiparam("nhalo");
  
  datstr = stropen(out,"w");           /* a dummy write ; should not fail */
  strclose(datstr);
  
  sprintf(rundir,"%s.tmpdir",out);     /* create and change to (tmp) rundir */
  make_rundir(rundir);
  goto_rundir(rundir);
  
  datstr = stropen("magalie.in","w");      /* create input file */
  
  fprintf(datstr,"%d             !! Random generator seed (ran1)\n",seed);
  fprintf(datstr,"n               ! Output particles smooth length?\n");
  fprintf(datstr,"%d             !! Number of disc particles (M_d=1)\n",ndisk);
  fprintf(datstr,"3.              ! Mass of disk (reset = 1 for default)\n");
  fprintf(datstr,"1.              ! Disc scalelength\n");
  fprintf(datstr,".1              ! Disc scale height ( = 1/5 length )\n");
  fprintf(datstr,"2.43            ! Solar radius\n");
  fprintf(datstr,"1.5             ! Toomre Q parameter @ Sun ( was 3/2 -> close to 1)\n");
  fprintf(datstr,".1              ! disc particle smoothing length (set to numerical resolution)\n");
  fprintf(datstr,"2.              ! z-max (= 10.*scale height)\n");
  fprintf(datstr,"10.             ! Max disc part radius\n");
  
  fprintf(datstr,"n               ! Gas in disc?\n");
  fprintf(datstr,"0               ! Number of disc gas particles\n");
  fprintf(datstr,"0.              ! mass of disc gas particles\n");
  fprintf(datstr,"0.              ! temperature of gas in disc\n");
  fprintf(datstr,"0.              ! scale height of gas in disc z0gas\n");
  fprintf(datstr,"0.              ! max z value of gas in disc (= 10 z0gas)\n");
  fprintf(datstr,"0.              ! max cylindrical radius of gas in disc\n");
  fprintf(datstr,"0.              ! min    "          "    "   "   "  "\n");
  fprintf(datstr,"n               ! include gas self-gravity (in disc)? [y/n]\n");
  
  fprintf(datstr,"y               ! Add bulge?\n");
  fprintf(datstr,"0.75            ! Mass of bulge\n");
  fprintf(datstr,"2.              ! Bulge scale length (hernquist 'a')\n");
  fprintf(datstr,"y               ! bulge self-gravity?\n");
  fprintf(datstr,"%d             !! N part in  bulge\n",nbulge);
  fprintf(datstr,"10.             ! max radius for bulge particles\n");
  fprintf(datstr,".01             ! softening length for particles\n");
  fprintf(datstr,"n               ! Non-spherical bulge?\n");
  fprintf(datstr,".89             ! Value of minor axis ratio ( c/a < 1 )\n");
  fprintf(datstr,"1.              ! Z-max value for bulge particles\n");
  fprintf(datstr,"5               ! Number of Simpson integration steps\n");
  fprintf(datstr,"n               ! Include rotation ? (flip z-momentum component)\n");
  fprintf(datstr,".0              ! Fraction of stars with aligned momentum (0<f<1)\n");
  
  fprintf(datstr,"y               ! Include a halo?\n");
  fprintf(datstr,"y               ! halo self-gravity ?\n");
  fprintf(datstr,"20.             ! max radius of halo\n");
  fprintf(datstr,"%d             !! Halo particle number\n",nhalo);
  fprintf(datstr,"0.01            ! halo particle softening length\n");
  fprintf(datstr,"LH              ! Halo type (LH or IS)\n");
  fprintf(datstr,"10.             ! Halo mass\n");
  fprintf(datstr,"2.              ! Halo core radius  (need two lengths if IS halo)\n");
  fprintf(datstr,"2.              ! Second length : dummy (LH) or truncation radius (IS)\n");
  fprintf(datstr,"n               ! Non-spherical halo?\n");
  fprintf(datstr,".5              ! aspect ratio (spheroid only)\n");
  fprintf(datstr,"n               ! Galactic satellite\n");
  fprintf(datstr,".1              ! Satellite mass\n");
  fprintf(datstr,"1.              ! Satellite scale length (LH model only)\n");
  fprintf(datstr,"5.,0.,0.        ! xyz position of the satellite in Galaxy frame of reference\n");
  fprintf(datstr,"5.              ! maximum satellite size (compare with sat length scale)\n");
  fprintf(datstr,"y               ! Include satellite self-gravity ? (yes = live model only)\n");
  fprintf(datstr,"0               ! number of satellite particles\n");
  fprintf(datstr,"0.1             ! smoothing length for satelllite particles\n");
  fprintf(datstr,"b               ! data format for output (a)scii or (b)inary\n");
  fprintf(datstr,"nbody           ! formatted as LH, nbody\n");
  fprintf(datstr,"m.dat           ! name of dataset (note - truncated at cr character)\n");
  strclose(datstr);
   
  histr = stropen("history","w");   /* maintain NEMO style history */
  put_history(histr);
  strclose(histr);
  
  datstr = stropen("make-it","w!"); /* create shell script to be run */
  fprintf(datstr,"#! /bin/sh\n");
  fprintf(datstr,"# created by NEMO's magalie wrapper program\n");
  fprintf(datstr,"magalie.exe < magalie.in\n");
  fprintf(datstr,"rm -f ../%s\n",out);
  fprintf(datstr,"unfio in=m.dat block=0 type=f | tabtos - ../%s block1=m,x,y,z,vx,vy,vz options=wrap nbody=%d\n",
	  out, ndisk+nbulge+nhalo);
  strclose(datstr);
  run_program("chmod +x make-it; ./make-it");   /* run it ! */
  if (Qcleanup) {
    dprintf(0,"Removing the run directory %s.tmpdir\n",out);  
    sprintf(rundir,"cd ..; rm -rf %s.tmpdir",out);  
    run_program(rundir);   
  }
}

void goto_rundir(string name)
{
  if (chdir(name))
    error("Cannot change directory to %s",name);
}

void make_rundir(string name)
{
  if (mkdir(name, 0755))
    error("Run directory %s already exists",name);
}

void run_program(string cmd)
{
  system(cmd);
}
