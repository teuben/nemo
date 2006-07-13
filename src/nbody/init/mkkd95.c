/*
 *  mkkd95:   make a Kuijken-Dubinski (1995) compostite B-D-H model
 *       This is simply a driver that runs the GalactICS routines
 *	
 *  6-mar-04    Created              PJT
 *  11-mar-04   NEMO style seed=  instead of iseed=    pjt
 *              added nmodel=
 *
 *  22-mar-04   V1.4  stdout/err from kd95 routines now to logfile     PJT
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
  "out=???\n        output snapshot (a rundirectory $out.tmpdir is also created)",
 
  "ndisk=8000\n     Number of particles in disk  (use 0 to skip this component)",
  "nbulge=4000\n    Number of particles in bulge (use 0 to skip this component)",
  "nhalo=6000\n     Number of particles in halo  (use 0 to skip this component)",

  "psi0=-4.6\n       in.dbh: Psi0    (HALO)",
  "v0=1.42\n         in.dbh: v0",
  "q=1\n             in.dbh: q",
  "rck2=0.1\n        in.dbh: (rc/rk)^2",
  "ra=0.8\n          in.dbh: ra",

  "md=0.867\n        in.dbh: M_d      (DISK)",
  "rd=1\n            in.dbh: R_d",
  "router=5\n        in.dbh: R_outer",
  "zd=0.1\n          in.dbh: z_d",
  "drtrunc=0.5\n     in.dbh: dR_trunc",

  "rhob=14.45\n      in.dbh: rho_b    (BULGE)",
  "psicut=-2.3\n     in.dbh: psi_cut",
  "sigb=0.714\n      in.dbh: sib_b",

  "dr=0.01\n         in.dbh: delta_r",
  "nr=2400\n         in.dbh: nr",
  "lmax=10\n         in.dbh: number of harmonics (even)",

  "sigvr0=0.47\n     in.diskdf: central radial velocity dispersion",
  "sigr0=1.0\n       in.diskdf: scalelength of sig_r^2",
  
  "ncorr=50\n        in.diskdf: number of intervals for correction function",
  "niter=10\n        in.diskdf: number of iterations",
  
  "fstreamb=0.75\n   in.bulge:",
  "fstreamh=0.5\n    in.halo:",
  
  "seed=0\n          in.- Random Seed (NEMO style) ",
  "zerocm=t\n        in.- Center the snapshot?",

  "bin=\n            directory in which KD95 binaries live (otherwise assume $PATH)",
  "model=\n          Select base model A, B, C or D ** not properly implemented yet **",
  "nmodel=1\n        Number of models to make",
  "cleanup=t\n       Cleanup the temporary rundir?",

  "VERSION=1.5\n     12-jul-06 PJT",
  NULL,
};

string usage="Kuijken-Dubinski-95 composite bulge-disk-halo model";

void goto_rundir(string name);
void make_rundir(string name);
void run_program(string cmd);
void model(char *m);


void nemo_main(void)
{
    char rundir[256], comment;
    stream datstr, histr;
    string out=getparam("out");
    string kd_bindir = getparam("bin");
    int seed, zerocm, nbulge, ndisk, nhalo;
    int imodel, nmodel = getiparam("nmodel");
    bool Qcleanup = getbparam("cleanup");

    if (hasvalue("model"))
      model(getparam("model"));           /*  patch up (putparam) if needed */

    seed =   -init_xrandom(getparam("seed"));  /* make sure it's negative */
    zerocm = getbparam("zerocm") ? 1 : 0;
 
    nbulge = getiparam("nbulge");
    ndisk  = getiparam("ndisk");
    nhalo  = getiparam("nhalo");

    dprintf(0,"mkkd95: Ndisk=%d Nbulge=%d Nhalo=%d\n",ndisk,nbulge,nhalo);

    datstr = stropen(out,"w");           /* a dummy write ; should not fail */
    strclose(datstr);

    sprintf(rundir,"%s.tmpdir",out);     /* create and change to (tmp) rundir */
    make_rundir(rundir);
    goto_rundir(rundir);

    datstr = stropen("in.dbh","w");      /* create input file */
    fprintf(datstr,"%s\n",  
	    nhalo > 0 ? "y" : "n");
    fprintf(datstr,"%g %g %g %g %g\n",
	    getdparam("psi0"),
	    getdparam("v0"),
	    getdparam("q"),
	    getdparam("rck2"),
	    getdparam("ra"));

    fprintf(datstr,"%s\n",  
	    ndisk > 0 ? "y" : "n");
    fprintf(datstr,"%g %g %g %g %g\n",
	    getdparam("md"),
	    getdparam("rd"),
	    getdparam("router"),
	    getdparam("zd"),
	    getdparam("drtrunc"));

    fprintf(datstr,"%s\n",  
	    nbulge > 0 ? "y" : "n");
    fprintf(datstr,"%g %g %g\n",
	    getdparam("rhob"),
	    getdparam("psicut"),
	    getdparam("sigb"));
    fprintf(datstr,"%g %d\n",
	    getdparam("dr"),
	    getiparam("nr"));
    fprintf(datstr,"%d\n",
	    getiparam("lmax"));
    fprintf(datstr,"%s\n",
	    "dbh.ps/ps");
    strclose(datstr);


    datstr = stropen("in.diskdf","w"); /* create input file */
    fprintf(datstr,"%g %g\n",
	    getdparam("sigvr0"),
	    getdparam("sigr0"));
    fprintf(datstr,"%d\n",
	    getiparam("ncorr"));
    fprintf(datstr,"%d\n",
	    getiparam("niter"));
    fprintf(datstr,"%s\n",
	    "diskdf.ps/ps");
    strclose(datstr);

    histr = stropen("history","w");/* maintain history */
    put_history(histr);
    strclose(histr);

    for (imodel=0; imodel<nmodel; imodel++) {   /* loop over making nmodel of them */
      dprintf(1,"Creating model %d\n",imodel+1);
      comment = (imodel==0) ? ' ' : '#';        /* comment out this out beyond 1st model */

      datstr = stropen("in.bulge","w!"); /* create input file for genbulge */
      fprintf(datstr,"%g\n",
	      getdparam("fstreamb"));
      fprintf(datstr,"%d\n",
	      nbulge);
      fprintf(datstr,"%d\n",
	      seed);
      fprintf(datstr,"%d\n",
	      zerocm);
      fprintf(datstr,"%s\n",
	      "dbh.dat");
      strclose(datstr);
      
      datstr = stropen("in.disk","w!"); /* create input file for gendisk*/
      fprintf(datstr,"%d\n",
	      ndisk);
      fprintf(datstr,"%d\n",
	      seed);
      fprintf(datstr,"%d\n",
	      zerocm);
      fprintf(datstr,"%s\n",
	      "dbh.dat");
      strclose(datstr);
      
      datstr = stropen("in.halo","w!"); /* create input file for genhalo */
      fprintf(datstr,"%g\n",
	      getdparam("fstreamh"));
      fprintf(datstr,"%d\n",
	      nhalo);
      fprintf(datstr,"%d\n",
	      seed);
      fprintf(datstr,"%d\n",
	      zerocm);
      fprintf(datstr,"%s\n",
	      "dbh.dat");
      strclose(datstr);
      
      
      datstr = stropen("make-it","w!"); /* create shell script to be run */
      fprintf(datstr,"#! /bin/csh -f\n");
      fprintf(datstr,"# created by NEMO's mkkd95 program\n");
      fprintf(datstr,"setenv _POSIX2_VERSION 1\n");
      if (*kd_bindir) {
	fprintf(datstr,"%c %s/dbh < in.dbh\n",comment,kd_bindir);
	fprintf(datstr,"%c %s/getfreqs\n",comment,kd_bindir);
	fprintf(datstr,"%c %s/diskdf < in.diskdf\n",comment,kd_bindir);
	fprintf(datstr,"%s/genbulge < in.bulge  > bulge\n",kd_bindir);
	fprintf(datstr,"%s/gendisk  < in.disk   > disk\n",kd_bindir);
	fprintf(datstr,"%s/genhalo  < in.halo   > halo\n",kd_bindir);
	fprintf(datstr,"%s/mergerv disk bulge halo > galaxy\n",kd_bindir);
      } else {
	fprintf(datstr,"%c dbh < in.dbh\n",comment);
	fprintf(datstr,"%c getfreqs\n",comment);
	fprintf(datstr,"%c diskdf < in.diskdf\n",comment);
	fprintf(datstr,"genbulge < in.bulge  > bulge\n");
	fprintf(datstr,"gendisk  < in.disk   > disk\n");
	fprintf(datstr,"genhalo  < in.halo   > halo\n");
	fprintf(datstr,"mergerv disk bulge halo > galaxy\n");
      }
      if (imodel==0) {
	fprintf(datstr,"rm -f ../%s\n",out);
	fprintf(datstr,"tabtos galaxy ../%s nbody,time mass,pos,vel headline=%d\n",out,seed);
      } else {
	fprintf(datstr,"tabtos galaxy - nbody,time mass,pos,vel headline=%d>> ../%s\n",seed,out);
      }
      fprintf(datstr,"echo DEBUG; cat in.bulge\n");
      strclose(datstr);
      
      run_program("chmod +x make-it; ./make-it > kd95.log 2>&1");   /* run it ! */
      
      seed = -init_xrandom(getparam("seed"));  /* make sure seed is negative again for KD95 */
    } /* imodel */
    if (Qcleanup) {     /* remove the run directory */
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

typedef struct _mpar {
  string name, par[4];
} mpar;


mpar ModelPars[] = {
  { "nbulge",  { "4000", "1000",  "2000", "1000" }},
  { "ndisk",   { "8000", "1000",  "4000", "1000" }},
  { "nhalo",   { "6000", "1000", "10000", "1000" }},

  { "fstreamb",{ "0.75", "0.5",   "0.5",  "0.5" }},
  { "fstreamh",{ "0.5",  "0.5",   "0.5",  "0.5" }},

  {  "psi0",   { "-4.6", "-5.202", "-6.0", "-7.0"}},
  {  "v0",     { "1.42", "1.36",   "1.32", "1.30"}},

  /* q     always 1   */
  /* rck2  always 0.1 */
  /* ra    always 0.8 */
  /* md    always 0.867 */
  /* rd    always 1 */
  /* router always 5 */
  /* zd     always 0.1 */
  /* drtrunc always 0.5 */


  /* rhob    always 14.45 */
  { "psicut",   { "-2.3",  "-2.89",  "-3.7",  "-4.7"}},
  /* sigb    always 0.714 */

  /* deltar  always 0.01 */
  { "nr",       {  "2400", "3200", "5000", "7500"}},


  /* lmax    always 10 */

  /* sigvr0   always 0.47 */
  /* sigr0    always 1.0 */
  
  /* ncorr    always 50 */
  /* niter    always 10 */

};

/*
 *  model: this sets new default values for a number
 *
 *    there is a bug in getparam.c - putparam() should
 *    only be done if the commandline version wasn't 
 *    used!!! 
 *         
 */

void model(char *m)
{
  int idx, i, n = sizeof(ModelPars)/sizeof(mpar);

  idx = -1;
  if (*m == 'A') idx=0;
  if (*m == 'B') idx=1;
  if (*m == 'C') idx=2;
  if (*m == 'D') idx=3;
  if (idx < 0) error("invalid model=%s\n",m);

  printf("We found %d; idx=%d\n",n,idx);
  for (i=0; i<n; i++)
    putparam(ModelPars[i].name, ModelPars[i].par[idx]);

}
