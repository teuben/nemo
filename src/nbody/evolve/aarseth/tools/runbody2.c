/*
 *  MAIN.C: NEMO driver for nbody2
 *	
 *	16-jun-97   V1.0    does the basics, nothing fancy      PJT
 *      18-jun-97   V1.1    more logic builtin (see input.f, data.f, scale.f)
 *	21-jun-97   V1.1a   added history
 *	15-jul-97   V1.2    changed out= to outdir=, no def for nbody=     PJT
 *                          added tnext=
 *	 5-mar-98   V1.3    added KZ_ parameters in the code		   PJT
 *      17-mar-06   V1.4    using fullname() for in=                       pjt
 *      17-sep-2013 V1.5    using new run interface                        PJT
 *       9-feb-2017 V2.0    better defaults for snapshot input             PJT
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <history.h>
#include <filefn.h>
#include <run.h>

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>


string defv[] = {
    "in=\n            input file (optional - see nbody= ) ",
    "outdir=???\n     output run directory (required)",

    "nbody=\n         Total particle number (<= NMAX) if in= not given",    
    "nfix=1\n         Output frequency of data save or binaries; KZ(3 & 6)",
    "nrand=0\n        Random number sequence skip",
    "nnbmax=40\n      Maximum number of neighbours (< LMAX)",
    "nrun=1\n         Run identification index",

    "etai=0.02\n      Time-step parameter for irregular force polynomial",
    "etar=0.04\n      Time-step parameter for regular force polynomial",
    "rs0=1.0\n        Initial radius of neighbour sphere",
    "deltat=0.5\n     Output time interval in units of the crossing time",
    "tcrit=5.0\n      Termination time in units of the crossing time",
    "tnext=\n         Next output time for restart runs",
    "qe=0.00002\n     Energy tolerance (restart if DE/E > 5*QE & KZ(2) > 1)",
    "eps=0.05\n       Potential softening parameter",


    "kz=1,2,1,2,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,0\n  Non-zero options for alternative paths (see below)\n"
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
     " 19  Use code units for tcrit/deltat\n"
     " 20  Not used at present",

    "xtpar1=\n        Mass of external Plummer model (KZ(15) = 1; scaled units)",
    "xtpar2=\n        Length scale for Plummer model (KZ(15) = 1)",
    "zmgas=\n         Mass scale for logarithmic potential (KZ(15) = 2)",
    "rgas=\n          Length scale for logarithmic potential (KZ(15) = 2)",

    "alphas=2.3\n     Power-law index for initial mass function",
    "body1=5.0\n      Maximum particle mass before scaling",
    "bodyn=1.0\n      Minimum particle mass before scaling",

    "q=0.0\n          Virial ratio (q=0.5 for virial equilibrium)",
    "vxrot=0.0\n      XY-velocity scaling factor (> 0 for solid-body rotation)",
    "vzrot=0.0\n      Z-velocity scaling factor (not used if VXROT = 0)",
    "rbar=1.0\n       Virial radius in pc (for scaling to physical units)",
    "zmbar=1.0\n      Mean mass in solar units",

    "xcm=\n           Displacement for subsystem (routine SCALE; KZ 17)",
    "ecc=\n           Eccentricity of relative motion for subsystem (ECC =< 1)",

    "kstart=1\n       Running mode (1=new 2=restart 3,4,5=restart w/ new par",
    "tcomp=40.0\n     Maximum allowed running time (minutes)",

    "KZ#=\n           [indexed] Override some kz= keywords",
    "exe=nbody2\n     Name of the executable",
    

    "VERSION=2.2\n    11-feb-2019 PJT",
    NULL,
};

string usage="front end to Ahmad-Cohen N-body code nbody2";

#define KZ_MAX  20

#define KZ_OUT   3
#define KZ_INI   4
#define KZ_EXT  15
#define KZ_MER  17
#define KZ_COM  18

void nemo_main(void)
{
    int nbody, nfix, nrand, nnbmax, nrun, kstart;
    real etai, etar, rs0, deltat, tcrit, qe, eps, tcomp;
    int k, nkz, kz[KZ_MAX];
    real alphas, body1, bodyn;
    real q, vxrot, vzrot, rbar, zmbar;
    string exefile = getparam("exe");
    string parfile = "nbody2.in";
    string rundir = getparam("outdir");
    string infile, fname;
    char dname[256], runcmd[256];
    stream datstr, histr;

    kstart = getiparam("kstart");
    tcomp =  getdparam("tcomp");

    nfix = getiparam("nfix");
    nrand = getiparam("nrand");
    nnbmax = getiparam("nnbmax");
    nrun = getiparam("nrun");

    etai = getdparam("etai");
    etar = getdparam("etar");
    rs0 = getdparam("rs0");
    deltat = getdparam("deltat");
    tcrit = getdparam("tcrit");
    qe = getdparam("qe");
    eps = getdparam("eps");

    nkz = nemoinpi(getparam("kz"),kz,20);
    if (nkz < 0) error("%d Syntax error kz=%s",nkz,getparam("kz"));
    for (k=nkz; k<KZ_MAX; k++) kz[k]=0;
    for (k=0; k<KZ_MAX; k++) dprintf(1,"%d ",kz[k]);
    dprintf(1,"\n");

    for (k=0; k<KZ_MAX; k++) {
      if (indexparam("KZ",k+1)) {
	dprintf(0,"KZ %d=%d\n",k+1,getiparam_idx("KZ",k+1));
	kz[k] = getiparam_idx("KZ",k+1);
      }
    }

    alphas = getdparam("alphas");
    body1 = getdparam("body1");
    bodyn = getdparam("bodyn");

    q = getdparam("q");
    vxrot = getdparam("vxrot");
    vzrot = getdparam("vzrot");
    rbar = getdparam("rbar");
    zmbar = getdparam("zmbar");

    infile = getparam("in");
    fname = fullname(infile);

    run_mkdir(rundir);

    sprintf(dname,"%s/%s",rundir,parfile);
    datstr = stropen(dname,"w");

    if (hasvalue("in")) {
      stream instr = stropen(getparam("in"),"r");
      nbody = get_snap_nbody(instr);
      strclose(instr);
      dprintf(0,"Grabbing nbody=%d\n",nbody);
    } else {
      if (hasvalue("nbody"))
	nbody = getiparam("nbody");
      else
	error("nbody= needs to be specified");
    }
    

    /*  New Run */

  if (kstart == 1) {

    fprintf(datstr,"%d %g\n",kstart,tcomp);
    fprintf(datstr,"%d %d %d %d %d\n",nbody,nfix,nrand,nnbmax,nrun);
    fprintf(datstr,"%g %g %g %g %g %g %g\n",etai,etar,rs0,deltat,tcrit,qe,eps);

    if (hasvalue("in") && kz[KZ_INI-1] != 2) {
      warning("Using in= and setting kz(%d)=2",KZ_INI);
      kz[KZ_INI-1] = 2;
    }

    for (k=0; k<KZ_MAX; k++) fprintf(datstr,"%d ",kz[k]);
    fprintf(datstr,"\n");

    if (kz[KZ_EXT-1] > 0) {           /* external potential parameters */
        if (kz[KZ_EXT-1]==1)
            fprintf(datstr,"%g %g\n",getdparam("xtpar1"),getdparam("xtpar2"));
        else
            fprintf(datstr,"%g %g\n",getdparam("zmgas"),getdparam("rgas"));
    }
            
    fprintf(datstr,"%g %g %g\n",alphas,body1,bodyn);
    fprintf(datstr,"%g %g %g %g %g\n",q,vxrot,vzrot,rbar,zmbar);

    if (kz[KZ_MER-1] > 0)
        fprintf(datstr,"%g %g\n",getdparam("xcm"), getdparam("ecc"));
    strclose(datstr);

    run_cd(rundir);
    histr = stropen("history","w");
    put_history(histr);
    strclose(histr);

    if (hasvalue("in")) {
	if (*infile == '-') {		/* do something special for pipes */
	  sprintf(runcmd,"stou4 %s",infile);
	} else {
	  sprintf(runcmd,"stou4 %s nbody=%d",fname,nbody);
        } 
        dprintf(0,"%s\n",runcmd);
        if (run_sh(runcmd)) error("Error converting input data");
    }

    sprintf(runcmd,"%s < %s",exefile,parfile);
    run_sh(runcmd);
  } else {
    error("kstart=%d not yet supported",kstart);
  }
}


