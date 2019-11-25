/*
 *  MAIN.C: NEMO driver for nbody1
 *	
 *  4-mar-98    cloned off runbody2, added KZ_ parameters     PJT
 *  17-mar-04   using in= didn't enable kz(4) properly        PJT
 *  17-mar-06   using fullname()                              PJT
 *  17-sep-2013  V1.3   using new run interface               PJT
 *   9-feb-2017  V2.0   better defaults for snapshot input    PJT
 *
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
    "nrun=1\n         Run identification index",

    "eta=0.02\n       Time-step parameter for irregular force polynomial",
    "deltat=0.25\n    Output time interval in units of the crossing time",
    "tcrit=2.0\n      Termination time in units of the crossing time",
    "qe=0.00002\n     Energy tolerance (restart if DE/E > 5*QE & KZ(2) > 1)",
    "eps=0.05\n       Potential softening parameter",


    "kz=1,2,1,2,0,1,0,0,0,1,1,0,0,0,1\n  Non-zero options for alternative paths (see below)\n"
      " 1  COMMON save on unit 1 if TCOMP > CPU or if TIME > TCRIT.\n"
      " 2  COMMON save on unit 2 at output (=1); restart if DE/E > 5*QE (=2).\n"
      " 3  Basic data written to unit 3 at output time (frequency NFIX).\n"
      " 4  Initial conditions on unit 4 (=1: output; =2: input).\n"
      " 5  Initial conditions (=0: uniform & isotropic; =1: Plummer).\n"
      " 6  Output of significant binaries.\n"
      " 7  Output of movie frames on unit 7.\n"
      " 8  Generation of two subsystems (merger experiment).\n"
      " 9  Individual bodies printed at output time (MIN(5**KZ9,N)).\n"
      "10  No scaling of initial conditions.\n"
      "11  Modification of ETA by tolerance QE.\n"
      "12  Initial parameters for binary orbit.\n"
      "13  Escaper removal (R > 2*RTIDE; RTIDE = 10*RSCALE).\n"
      "14  Adjustment of coordinates & velocities to c.m. condition.\n"
      "15  Use code units for tcrit/deltat",
    
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
    "exe=nbody1\n     Name of the executable",
    
    "VERSION=2.2\n    6-mar-2019 PJT",
    NULL,
};

string usage="front end for Aarseth nbody1 N-body code";

#define KZ_MAX  15


#define KZ_OUT   3
#define KZ_INI   4
#define KZ_MOV   7
#define KZ_MER   8
#define KZ_BIN  12
#define KZ_COM  14

void nemo_main(void)
{
    int nbody, nfix, nrand, nrun, kstart;
    real eta, deltat, tcrit, qe, eps, tcomp;
    int k, nkz, kz[KZ_MAX];
    real alphas, body1, bodyn;
    real q, vxrot, vzrot, rbar, zmbar;
    string exefile = getparam("exe");
    string parfile = "nbody1.in";
    string rundir = getparam("outdir");
    string infile, fname;
    char dname[256], runcmd[256];
    stream datstr, histr;

    kstart = getiparam("kstart");
    tcomp =  getdparam("tcomp");

    nfix = getiparam("nfix");
    nrand = getiparam("nrand");
    nrun = getiparam("nrun");

    eta = getdparam("eta");
    deltat = getdparam("deltat");
    tcrit = getdparam("tcrit");
    qe = getdparam("qe");
    eps = getdparam("eps");

    nkz = nemoinpi(getparam("kz"),kz,KZ_MAX);
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
    fprintf(datstr,"%d %d %d %d\n",nbody,nfix,nrand,nrun);
    fprintf(datstr,"%g %g %g %g %g\n",eta,deltat,tcrit,qe,eps);

    if (hasvalue("in") && kz[KZ_INI-1] != 2) {
      warning("Using in= and setting kz(%d)=2",KZ_INI);
      kz[KZ_INI-1] = 2;
    }

    for (k=0; k<KZ_MAX; k++) fprintf(datstr,"%d ",kz[k]);
    fprintf(datstr,"\n");

    if (kz[KZ_INI-1] != 2) {
        fprintf(datstr,"%g %g %g\n",alphas,body1,bodyn);
        if (hasvalue("in")) warning("Also gave in=, will be ignored");
    } else
    	if (!hasvalue("in")) warning("Did not supply in=");
    fprintf(datstr,"%g %g %g %g %g\n",q,vxrot,vzrot,rbar,zmbar);

    if (kz[KZ_BIN-1] > 0)
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
        if (system(runcmd)) error("Error converting input data");
	// @todo   grab nbody
    }

    sprintf(runcmd,"%s < %s",exefile,parfile);
    run_sh(runcmd);
    
    sprintf(runcmd,"u3tos in=OUT3 out=OUT3.snap mode=1");
    run_sh(runcmd);    
    
  } else {
    error("kstart=%d not yet supported",kstart);
  }
}

/*
 *	Order of input lines in "nbody1.in" for a new run (KSTART=1)
 *
 *          variables                   condition               where
 *  ----------------------------    -------------------
 *  KSTART TCOMP                                            nbody1.f    MAIN
 *  Nbody NFIX NRAND NRUN                                   input.f     INPUT
 *  ETA DELTAT TCRIT QE EPS                                 input.f     INPUT
 *  KZ(1..15)                                               input.f     INPUT
 *      ALPHAS BODY1 BODYN          KZ(4).NE.2              data.f      DATA
 *      SEMI ECC                    KZ(12).NE.0             data.f      DATA
 *  Q VXROT VZROT RBAR ZMBAR                                scale.f     SCALE
 *      NFRAME DELTAF               KZ(7).GT.0              scale.f     SCALE
 *      XCM ECC                     KZ(8).GT.0              subsys.f    SUBSYS
 *
 */
 
