/*  magalie:   make a galaxy (MaGalie; Boily et al. prescription)
 *
 *  <dark-ages> Original fortran code               Boily et al.
 *  21-mar-04   1.1 interface created in rainy Strasbourg           PJT/CB
 *  24-mar-04   1.2 added most primary keywords, at 37,000ft        PJT
 *  23-jan-05   1.2a   fixed bulge mass encoding error     - courtesy J.J.Fleck
 *  24-mar-06   1.2b   fixed bulge radius encoding error   PJT
 *  12-jul-06   1.2c   merged two versions - PJT
 *  19-jul-06   1.3    use header=   PJT

 * TODO: merge in history of how magalie was called...
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

  "hdisk=0.2\n    Scaleheight of the disk (in units of scale length)",
  "rsolar=2.43\n  Solaris radius for calibration of Vsolar=220",
  "Qtoomre=1.5\n  Toomre stability parameter",
  "cdisk=10\n     Cutoff radius of disk",
  "zdisk=2\n      Height cutoff of disk",

  "mbulge=0.25\n  Mass of bulge (in units of disk mass)",
  "rbulge=2\n     Scalelength of bulge",
  "cbulge=10\n    Cutoff radius of bulge",
  "zbulge=1\n     Cutoff height of bulge",
  "qbulge=\n      Axis ratio C/A of bulge, in case flattened",
  "fbulge=\n      Streaming fraction of bulge, in case streaming",

  "mhalo=5\n      Mass of halo (in units of disk mass)",
  "rhalo=2\n      Scalelength of halo",
  "chalo=20\n     Cutoff radius of halo",
  "c2halo=20\n    Cutoff radius of halo for ISO (for LH it is c2halo)",
  "qhalo=\n      Axis ratio C/A of halo, in case flattened",
  "type=LH\n      Type of halo: LH or ISO",

  "seed=0\n       Random seed",
  "cleanup=t\n    cleanup run directory after use",
  "header=\n      use an explicit unfio header size of 4 or 8",
  "VERSION=1.3a\n 22-jul-06 PJT",
  NULL,
};

string usage="Boily et al. composite bulge-disk-halo model";

string cvsid="$Id$";

void goto_rundir(string name);
void make_rundir(string name);
void run_program(string cmd);

double getdparamq(string,double);
string getsparamq(string);

void nemo_main(void)
{
  char rundir[256];
  stream datstr, histr;
  string out=getparam("out");
  char hdrkey[32];
  int seed, nbulge, ndisk, nhalo;
  bool Qcleanup = getbparam("cleanup");

  seed =   init_xrandom(getparam("seed")); 
  if (hasvalue("header")) 
    sprintf(hdrkey,"header=%d", getiparam("header"));
  else
    sprintf(hdrkey," ");
  
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
  fprintf(datstr,"1.              ! Mass of disk (reset = 1 for default)\n");
  fprintf(datstr,"1.              ! Disc scalelength\n");
  fprintf(datstr,"%g             !! Disc scale height ( = 1/5 length )\n",getdparam("hdisk"));
  fprintf(datstr,"%g             !! Solar radius\n",getdparam("rsolar"));
  fprintf(datstr,"%g             !! Toomre Q parameter @ Sun ( was 3/2 -> close to 1)\n",getdparam("Qtoomre"));
  fprintf(datstr,".1              ! disc particle smoothing length (set to numerical resolution)\n");
  fprintf(datstr,"%g             !! z-max (= 10.*scale height)\n",getdparam("zdisk"));
  fprintf(datstr,"%g             !! Max disc part radius\n",getdparam("cdisk"));
  
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
  fprintf(datstr,"%g             !! Mass of bulge\n",getdparam("mbulge"));
  fprintf(datstr,"%g             !! Bulge scale length (hernquist 'a')\n",getdparam("rbulge"));
  fprintf(datstr,"y               ! bulge self-gravity?\n");
  fprintf(datstr,"%d             !! N part in  bulge\n",nbulge);
  fprintf(datstr,"%g             ! max radius for bulge particles\n",getdparam("cbulge"));
  fprintf(datstr,".01             ! softening length for particles\n");
  fprintf(datstr,"%s             *! Non-spherical bulge?\n", getsparamq("qbulge"));
  fprintf(datstr,"%g             *! Value of minor axis ratio ( c/a < 1 )\n",getdparamq("qbulge",1.0));
  fprintf(datstr,"%g             !! Z-max value for bulge particles\n",getdparam("zbulge"));
  fprintf(datstr,"5               ! Number of Simpson integration steps\n");
  fprintf(datstr,"%s             *! Include rotation ? (flip z-momentum component)\n",getsparamq("fbulge"));
  fprintf(datstr,"%g             *! Fraction of stars with aligned momentum (0<f<1)\n",getdparamq("fbulge",0.0));
  
  fprintf(datstr,"y               ! Include a halo?\n");
  fprintf(datstr,"y               ! halo self-gravity ?\n");
  fprintf(datstr,"%g             !! max radius of halo\n",getdparam("chalo"));
  fprintf(datstr,"%d             !! Halo particle number\n",nhalo);
  fprintf(datstr,"0.01            ! halo particle softening length\n");
  fprintf(datstr,"%s             !! Halo type (LH or IS)\n",getparam("type"));
  fprintf(datstr,"%g             !! Halo mass\n",getdparam("mhalo"));
  fprintf(datstr,"%g             !! Halo core radius  (need two lengths if IS halo)\n",getdparam("c2halo"));
  fprintf(datstr,"2.              ! Second length : dummy (LH) or truncation radius (IS)\n");
  fprintf(datstr,"%s             *! Non-spherical halo?\n",getsparamq("qhalo"));
  fprintf(datstr,"%g             *! aspect ratio (spheroid only)\n",getdparamq("qhalo",1.0));

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
  fprintf(datstr,"# created by NEMO's magalie wrapper program VERSION=%s\n",getparam("VERSION"));
  fprintf(datstr,"magalie.exe < magalie.in\n");
  fprintf(datstr,"rm -f ../%s\n",out);
  fprintf(datstr,"unfio in=m.dat block=0 type=f %s | tabtos - ../%s block1=m,x,y,z,vx,vy,vz options=wrap nbody=%d\n",
	  hdrkey, out, ndisk+nbulge+nhalo);
  strclose(datstr);
  run_program("chmod +x make-it; ./make-it > magalie.log 2>&1");   /* run it ! */
  if (Qcleanup) {
    dprintf(1,"Removing the run directory %s.tmpdir\n",out);  
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
  if (system(cmd))
    warning("Some error occurred running: %s",cmd);
}

double getdparamq(string key, double def)
{
  return hasvalue(key) ? getdparam(key) : def ;
}

string getsparamq(string key)
{
  return hasvalue(key) ? "y" : "n";
}
