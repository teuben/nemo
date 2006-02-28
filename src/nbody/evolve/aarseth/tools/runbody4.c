/*
 *  MAIN.C: NEMO driver for nbody4
 *	
 *  28-feb-06    V0.1     Created, in Cambridge           PJT
 */

#include <stdinc.h>
#include <getparam.h>

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

string defv[] = {
    "in=\n            input file (optional - see nbody= ) ",
    "outdir=???\n     output run directory (required)",

    "nbody=1000\n     Total particle number (<= NMAX).",
    "nfix=1\n         Output frequency of data save or binaries; KZ(3 & 6)",
    "ncrit=5\n        Final particle number (alternative termination criterion)",
    "nrand=50000\n    Random number sequence skip",
    "nrun=1\n         Run identification index",

    "eta=0.02\n       Time-step parameter for irregular force polynomial",
    "dtadj=2.0\n      check",
    "deltat=10.0\n    Output time interval in units of the crossing time",
    "tcrit=100.0\n    Termination time in units of the crossing time",
    "qe=2e-5\n        Energy tolerance (restart if DE/E > 5*QE & KZ(2) > 1)",
    "rbar=1.0\n       ?",
    "zmbar=0.5\n      ?",


    "kz=0 1 0 0 1 0 1 0 0 0  0 0 0 1 1 1 0 0 3 4  1 0 2 0 0 1 0 0 0 1  0 0 0 0 0 0 1 0 0 0\n"
      "    Non-zero options for alternative paths (see below)\n"
      " 1  COMMON save on unit 1 at end of run (=2: every 100*NMAX steps)\n"
      " 2  COMMON save on unit 2 at output (=1); restart if DE/E > 5*QE (=2)\n"
      " 3  Basic data written to unit 3 at output time (frequency NFIX)\n"
      " 4  Initial conditions on unit 4 (=1: output; =2: input)\n"
      " 5  Initial conditions (=0: uniform & isotropic; =1: Plummer model)\n"
      " 6  Output of significant binaries (=2: frequency NFIX)\n"
      " 7  Lagrangian radii (=1: unit 6; =2: unit 7; =3: both types)\n"
      " 8  Core radius & density centre (N > 20 only)\n"
      " 9  Individual bodies printed at output time (MIN(5**KZ9,N))\n"
     " 10  No unique density centre (skips velocity modification of RS(I))\n"
     " 11  Modification of ETAI & ETAR by tolerance QE\n"
     " 12  Inelastic mergers (>1: diagnostic output)\n"
     " 13  Escaper removal (R > 2*RTIDE; RTIDE = 10*RSCALE)\n"
     " 14  Skip full predictor loop if NNB > KZ(14) = <NNB> & KZ(14) > 0\n"
     " 15  External potential (=1: Plummer model; =2: logarithmic potential)\n"
     " 16  No scaling of initial conditions\n"
     " 17  Generation of two subsystems (merger experiment)\n"
     " 18  Adjustment of coordinates & velocities to c.m. condition\n"
     " 19  Not used at present (same for # 20)\n"
     " 20  Not used at present (same for # 20)\n"
     " 21  Not used at present (same for # 20)\n"
     " 22  Not used at present (same for # 20)\n"
     " 23  Not used at present (same for # 20)\n"
     " 24  Not used at present (same for # 20)\n"
     " 25  Not used at present (same for # 20)\n"
     " 26  Not used at present (same for # 20)\n"
     " 27  Not used at present (same for # 20)\n"
     " 28  Not used at present (same for # 20)\n"
     " 29  Not used at present (same for # 20)\n"
     " 30  Not used at present (same for # 20)\n"
     " 31  Not used at present (same for # 20)\n"
     " 32  Not used at present (same for # 20)\n"
     " 33  Not used at present (same for # 20)\n"
     " 34  Not used at present (same for # 20)\n"
     " 35  Not used at present (same for # 20)\n"
     " 36  Not used at present (same for # 20)\n"
     " 37  Not used at present (same for # 20)\n"
     " 38  Not used at present (same for # 20)\n"
     " 39  Not used at present (same for # 20)\n"
     " 40  Not used at present (same for # 20)",
 
    "dtmin=1e-4\n     ",
    "rmin=1e-3\n      ",
    "etau=0.2\n       ",
    "eclose=1.0\n     ",
    "gmin=1e-6\n      ",
    "gmax=0.001\n     ",

    "alpha=2.3\n      Power-law index for initial mass function",
    "body1=10.0\n     Maximum particle mass before scaling",
    "bodyn=0.2\n      Minimum particle mass before scaling",
    "nbin0=0\n        ",
    "zmet=0.02\n      ",
    "epoch0=0\n       ",
    "dtplot=10.0\n    ",

    "q=0.0\n          Virial ratio (q=0.5 for virial equilibrium)",

    "apo=6.0\n        kz(5) = 2 or 3",
    "ecc=0.5\n        Eccentricity of relative motion for subsystem (ECC =< 1)",
    "n2=500\n         kz(5) = 2",
    "scale=0.5\n      kz(5) = 2 or 3",
    "dmin=3\n         kz(5) = 3",
    "semi=1e-3\n      ",
    "m1=25.0\n        ",
    "m2=25.0\n        ",
    "ratio=0.0\n      ",
    "range=100.0\n    ",
    
    "kstart=1\n       Running mode (1=new 2=restart 3,4,5=restart w/ new par",
    "tcomp=10.0\n     Maximum allowed running time (minutes) **ignored **",
    "gpid=0\n         Grape ID",

    "VERSION=0.1\n    28-feb-06 PJT",
    NULL,
};

string usage="Hermite N-body code";

#define KZ_MAX  40

#define KZ_OUT   3
#define KZ_INI   4
#define KZ_EXT  15
#define KZ_MER  17
#define KZ_COM  18

nemo_main()
{
  int nbody, nfix, ncrit, nrand, nnbmax, nrun, kstart, gpid, nbin0, n2;
  real eta, dtadj, deltat, tcrit, qe, eps, tcomp;
  int k, nkz, kz[KZ_MAX];
  real dtmin,rmin,etau,eclose,gmin,gmax;
  real alpha, body1, bodyn, zmet, epoch0, dtplot;
  real q, rbar, zmbar;
  real apo, ecc, dmin, semi, scale, m1, m2, ratio, range;

  string exefile = "nbody4";
  string parfile = "nbody4.in";
  string rundir = getparam("outdir");
  string infile;
  char fullname[256], runcmd[256];
  stream datstr, histr;

  kstart = getiparam("kstart");
  tcomp =  getdparam("tcomp");
  gpid = getdparam("gpid");
  
  nbody = getiparam("nbody");
  nfix = getiparam("nfix");
  ncrit = getiparam("ncrit");
  nrand = getiparam("nrand");
  nrun = getiparam("nrun");
  
  eta = getdparam("eta");
  dtadj = getdparam("dtadj");
  deltat = getdparam("deltat");
  tcrit = getdparam("tcrit");
  qe = getdparam("qe");
  rbar =  getdparam("rbar");
  zmbar =  getdparam("zmbar");
  
  nkz = nemoinpi(getparam("kz"),kz,KZ_MAX);
  if (nkz < 0) error("%d Syntax error kz=%s",nkz,getparam("kz"));
  for (k=nkz; k<KZ_MAX; k++) kz[k]=0;
  for (k=0; k<KZ_MAX; k++) dprintf(1,"%d ",kz[k]);
  dprintf(1,"\n");

  dtmin = getdparam("dtmin");
  rmin = getdparam("rmin");
  etau = getdparam("etau");
  eclose= getdparam("eclose");
  gmin = getdparam("gmin");
  gmax = getdparam("gmax");

  
  alpha = getdparam("alpha");
  body1 = getdparam("body1");
  bodyn = getdparam("bodyn");
  nbin0 = getiparam("nbin0");
  zmet = getdparam("zmet");
  epoch0 = getiparam("epoch0");
  dtplot = getiparam("dtplot");

  /* optional things, depending on kz5 */
  apo = getdparam("apo");
  ecc = getdparam("ecc");
  n2 = getiparam("n2");
  dmin = getdparam("dmin");
  scale = getdparam("scale");
  semi = getdparam("semi");
  m1 = getdparam("m1");
  m2 = getdparam("m2");
  ratio = getdparam("ratio");
  range = getdparam("range");



  make_rundir(rundir);

  sprintf(fullname,"%s/%s",rundir,parfile);
  datstr = stropen(fullname,"w");    
  
    /*  New Run */
  
  if (kstart == 1) {
    
    fprintf(datstr,"%d %g %d\n",kstart,tcomp,gpid);
    fprintf(datstr,"%d %d %d %d %d\n",nbody,nfix,ncrit,nrand,nrun);
    fprintf(datstr,"%g %g %g %g %g %g %g\n",eta,dtadj,deltat,tcrit,qe,rbar,zmbar);

    for (k=0; k<KZ_MAX; k++) {
      fprintf(datstr,"%d ",kz[k]);
      if (k%10 == 9)     fprintf(datstr,"\n");
    }

    fprintf(datstr,"%g %g %g %g %g %g\n",dtmin,rmin,etau,eclose,gmin,gmax);
    fprintf(datstr,"%g %g %g %d %g %g %g\n",alpha,body1,bodyn,nbin0,zmet,epoch0,dtplot);

    if (kz[4] == 2)
      fprintf(datstr,"%g %g %g %g %g %g\n",apo,ecc,n2,scale);
    else if (kz[4] == 3) 
      fprintf(datstr,"%g %g %g %g %g %g\n",apo,ecc,dmin,scale);
    else if (kz[4] == 4) 
      fprintf(datstr,"%g %g %g %g %g %g\n",semi,ecc,m1,m2);

    fprintf(datstr,"%g 0  0 0 \n",q);

    if (kz[7] == 1 || kz[7] == 3)
      fprintf(datstr,"%g %g %g %g 0 0 0 \n",semi,ecc,ratio,range);




    strclose(datstr);

    goto_rundir(rundir);
    histr = stropen("history","w");
    put_history(histr);
    strclose(histr);

    if (hasvalue("in")) {
	infile = getparam("in");
	if (*infile == '-') {		/* do something special for pipes */
            sprintf(runcmd,"stou4 %s",infile);
	} else if (*infile == '/') {	/* regular file */
            sprintf(runcmd,"stou4 %s nbody=%d",infile,nbody);
        } else
            error("in=%s must be an absolute pathname or special file",infile);
        dprintf(0,"%s\n",runcmd);
        if (system(runcmd)) error("Error converting input data");
    }

    sprintf(runcmd,"%s < %s",exefile,parfile);
    run_program(runcmd);
  } else {
    error("kstart=%d not yet supported for NBODY4",kstart);
  }
}


goto_rundir(string name)
{
    if (chdir(name))
        error("Cannot change directory to %s",name);
}

make_rundir(string name)
{
    if (mkdir(name, 0755))
        warning("Run directory %s already exists",name);
}

run_program(string cmd)
{
    system(cmd);
}
