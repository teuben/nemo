/*
 *  mkkd95:   make a Kuijken-Dubinski (1995) compostite B-D-H model
 *       This is simply a driver that runs the GalactICS routines
 *	
 *  6-mar-04    Created              PJT
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
  "out=???\n        output snapshot (a rundirectory $out.tmpdir is also created)",
 
  "nbulge=4000\n    Number of particles in bulge (use 0 to skip this component)",
  "ndisk=8000\n     Number of particles in disk  (use 0 to skip this component)",
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
  
  "iseed=-1\n        in.- Random Seed (kd95 style) ",
  "zerocm=t\n        in.- Center the snapshot?",

  "bin=.\n           directory in which KD95 binaries live",
  "model=A\n         Select base model A, B, C or D",

  "VERSION=1.1\n     6-mar-04 PJT",
  NULL,
};

string usage="Kuijken-Dubinsky-95 composite bulge-disk-halo model";

void goto_rundir(string name);
void make_rundir(string name);
void run_program(string cmd);
void model(char *m);


void nemo_main(void)
{
    char rundir[256];
    stream datstr, histr;
    string out=getparam("out");
    string kd_bindir = getparam("bin");
    int iseed, zerocm, nbulge, ndisk, nhalo;

    model(getparam("model"));           /*  patch up (putparam) if needed */

    iseed = getiparam("iseed");
    zerocm = getbparam("zerocm") ? 1 : 0;
 
    nbulge = getiparam("nbulge");
    ndisk  = getiparam("ndisk");
    nhalo  = getiparam("nhalo");

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

    datstr = stropen("in.bulge","w"); /* create input file */
    fprintf(datstr,"%g\n",
	    getdparam("fstreamb"));
    fprintf(datstr,"%d\n",
	    getiparam("nbulge"));
    fprintf(datstr,"%d\n",
	    iseed);
    fprintf(datstr,"%d\n",
	    zerocm);
    fprintf(datstr,"%s\n",
	    "dbh.dat");
    strclose(datstr);

    datstr = stropen("in.disk","w"); /* create input file */
    fprintf(datstr,"%d\n",
	    getiparam("ndisk"));
    fprintf(datstr,"%d\n",
	    iseed);
    fprintf(datstr,"%d\n",
	    zerocm);
    fprintf(datstr,"%s\n",
	    "dbh.dat");
    strclose(datstr);

    datstr = stropen("in.halo","w"); /* create input file */
    fprintf(datstr,"%g\n",
	    getdparam("fstreamh"));
    fprintf(datstr,"%d\n",
	    getiparam("nhalo"));
    fprintf(datstr,"%d\n",
	    iseed);
    fprintf(datstr,"%d\n",
	    zerocm);
    fprintf(datstr,"%s\n",
	    "dbh.dat");
    strclose(datstr);


    histr = stropen("history","w");/* maintain history */
    put_history(histr);
    strclose(histr);

    datstr = stropen("make-it","w"); /* create shell script to be run */
    fprintf(datstr,"#! /bin/sh\n");
    fprintf(datstr,"# created by NEMO's mkkd95 program\n");
    fprintf(datstr,"%s/dbh < in.dbh\n",kd_bindir);
    fprintf(datstr,"%s/getfreqs\n",kd_bindir);
    fprintf(datstr,"%s/diskdf < in.diskdf\n",kd_bindir);
    fprintf(datstr,"%s/genbulge < in.bulge  > bulge\n",kd_bindir);
    fprintf(datstr,"%s/gendisk  < in.disk   > disk\n",kd_bindir);
    fprintf(datstr,"%s/genhalo  < in.halo   > halo\n",kd_bindir);
    fprintf(datstr,"%s/mergerv disk bulge halo > galaxy\n",kd_bindir);
    fprintf(datstr,"rm ../%s\n",out);
    fprintf(datstr,"tabtos galaxy ../%s nbody,time mass,pos,vel\n",out);
    strclose(datstr);

    run_program("chmod +x make-it; ./make-it");   /* run it ! */
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
  { "nbulge",  { "4000", "1000",  "2000", "2000" }},
  { "ndisk",   { "8000", "1000",  "4000", "1000" }},
  { "nhalo",   { "6000", "1000",  "1000", "1000" }},
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
