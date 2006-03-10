/*
 *  MAIN.C: NEMO driver for nbody4
 *	
 *  28-feb-06    V0.1     Created, in Cambridge           PJT
 *   8-mar-06    V0.3     clean up keywords, create Nbody on the fly     PJT
 *  10-mar-06    V0.4     unit10 now binary, parallel to nbody1,2
 *                         
 *  TODO:   only a few special kz() cases were implemented. See nbody4.5 for the full set
 * 
 * VJ bench:       492", code halted tol exceeded at T=40  (3.2 GHz P4)
 * Code bench:    1448"                                                         2642" (laptop, P1.6GHz)
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

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
    "rbar=1.0\n       mean radius of system?",
    "zmbar=0.5\n      mean mass of system, in solar units",

#if 0
    "kz=0 1 0 0 1 0 1 0 0 0  0 0 0 1 1 1 0 0 3 4  1 0 2 0 0 1 0 0 0 1  0 0 0 0 0 0 1 0 0 0\n"
#else
    "kz=0 0 1 0 1 0 5 0 0 0  0 0 0 0 1 0 0 0 3 0  1 0 2 0 1 2 0 0 0 2  0 0 0 0 0 0 1 0 0 0\n",
#endif
      "Non-zero options for alternative paths (see below)\n"
      "       1  COMMON save on unit 1 at end of run (=2: every 100*NMAX steps).\n"
      "       2  COMMON save on unit 2 at output (=1); restart if DE/E > 5*QE (=2).\n"
      "       3  Basic data on unit 3 at output (freq. NFIX; >1: cluster + tail).\n"
      "     # 4  Binary diagnostics on unit 4 (# threshold levels = KZ(4) < 10).\n"
      "       5  Initial conditions (#22 =0; =0: uniform & isotropic sphere;\n"
      "                =1: Plummer; =2: two Plummer models in orbit, extra input;\n"
      "                =3: massive perturber and planetesimal disk, extra input).\n"
      "                =4: massive initial binary, extra input; output on unit 35).\n"
      "       6  Output of significant & regularized binaries (=1, 2, 3 & 4).\n"
      "       7  Lagrangian radii (>0: RSCALE; =2, 3, 4: output units 6, 7, 12;\n"
      "                >=5: density & rms velocity at given radii on unit 26 & 27;\n"
      "                 =6: Lagrangian radii for two mass groups on unit 31 & 32.\n"
      "       8  Primordial binaries (=1 & >=3; >0: BINOUT; >2: BINDAT; >3: HIDAT).\n"
      "       9  Individual bodies printed at output time (MIN(5**KZ9,NTOT)).\n"
      "      10  Diagnostic KS output (>0: begin; >1: end; >=3: each step).\n"
      "    # 11  Synchronization of circular orbits (suppressed).\n"
      "      12  Disk shocks (=1: standard model) or interstellar clouds (< 0).\n"
      "      13  Scaling of time (1: variable by t_cr; 2: variable by t_r;\n"
      "                 -1: constant scaling to t_r; -2: constant scaling to t_c).\n"
      "      14  External force (=1: linearized; -1: cutoff; =2: point-mass galaxy;\n"
      "             =3: point-mass + disk + logarithmic halo in any combination).\n"
      "      15  Triple, quad, chain (#30 > 0) or merger search (>1: full output).\n"
      "      16  Updating of regularization parameters (RMIN, DTMIN & ECLOSE).\n"
      "      17  Modification of ETA (>=1) & ETAU (>1) by tolerance QE.\n"
      "      18  Hierarchical systems (=1: diagnostics; =2: primordial; =3: both).\n"
      "      19  Stellar evolution and mass loss (=1: old supernova scheme;\n"
      "                      =3: Eggleton, Tout & Hurley; >4: Chernoff--Weinberg).\n"
      "      20  Initial mass function (=0,1: Salpeter; >1: various, see IMF).\n"
      "      21  Extra output (>0: model, etc; >1: CENTRE; >2: MTRACE; >3: GLOBAL).\n"
      "      22  Initial conditions on unit 10 (=1: output; =2,3(unscaled): input).\n"
      "      23  Escaper removal (>1: diagnostics in file ESC; =2: angles unit #6;\n"
      "                           >1: tidal tails output if #14 = 3).\n"
      "      24  Initial conditions for subsystem (routine SCALE; KZ(24) = #).\n"
      "    # 25  Partial reflection of KS binary orbit (GAMMA < GMIN; suppressed).\n"
      "      25  HR diagnostics of evolving stars (output of B & S on #82 & 83).\n"
      "      26  Slow-down of two-body motion (=1: KS binary; =2: chain binary).\n"
      "      27  Two-body interactions (-2: RADIUS = 0; -1: collision detection;\n"
      "                                 =1: sequential circ; > 0: collision).\n"
      "      28  (not used).\n"
      "    # 29  Boundary reflection for hot system (suppressed).\n"
      "      30  Chain regularization (=1: basic; >1: main output; >2: each step).\n"
      "      31  Centre of mass correction after energy check.\n"
      "      32  Increase of output intervals (based on single particle energy).\n"
      "      33  Block-step diagnostics at main output (=2: active pipes).\n"
      "    # 34  Roche lobe overflow (suppressed).\n"
      "      35  Time offset (global time from TTOT = TIME + DTOFF; offset = 100).\n"
      "      36  Step reduction for hierarchical systems (not recommended).\n"
      "      37  Step reduction for encounters with high-velocity particles.\n"
      "    # 38  Multiple use of GRAPE-6 (sleep 1 sec after each timer check).\n"
      "      39  Neighbour list (=-1: on host; =0: full list or closest on GRAPE).\n"
      "      40  (not used).\n"
      "    # Currently supressed",

    "dtmin=1e-5\n     time-step criterion for regularization search",
    "rmin=1e-4\n      distance criterion for regularization search",
    "etau=0.2\n       ",
    "eclose=1.0\n     ",
    "gmin=1e-6\n      ",
    "gmax=0.001\n     ",

    "alpha=2.3\n      Power-law index for initial mass function",
    "body1=10.0\n     Maximum particle mass before scaling",
    "bodyn=0.2\n      Minimum particle mass before scaling",
    "nbin0=0\n        number of primordial binaries",
    "zmet=0.02\n      ",
    "epoch0=0\n       epoch of star formation (in Myr)",
    "dtplot=10.0\n    plotting interval (nbody units)",

    "q=0.5\n          Virial ratio (q=0.5 for virial equilibrium)",

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
    "gpid=0\n         Number of grape boards (starting from 0) ** if applicable **",

    "format=%g\n      Format used for fort.10 input conditions if in= used",
    "KZ#=\n           [indexed] Override some kz= keywords",

    "VERSION=0.4\n    10-mar-06 PJT",
    NULL,
};

string usage="Hermite N-body code";

string cvsid="$Id$";

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
  int k, nkz, kz[KZ_MAX], KZ[KZ_MAX];
  int maxidx, bits;
  real dtmin,rmin,etau,eclose,gmin,gmax, tsnap;
  real alpha, body1, bodyn, zmet, epoch0, dtplot;
  real q, rbar, zmbar;
  real apo, ecc, dmin, semi, scale, m1, m2, ratio, range;
  Body *bp, *btab = NULL;

  string exefile = "nbody4";
  string parfile = "nbody4.in";
  string rundir = getparam("outdir");
  string fmt = getparam("format");
  string infile;
  char fullname[256], runcmd[256], fmt7[256];
  stream datstr, histr, instr, outstr;

  kstart = getiparam("kstart");
  tcomp =  getdparam("tcomp");
  gpid = getdparam("gpid");
  
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

  for (k=0; k<KZ_MAX; k++) {
    if (indexparam("KZ",k+1)) {
      dprintf(0,"KZ %d=%d\n",k+1,getiparam_idx("KZ",k+1));
	kz[k] = getiparam_idx("KZ",k+1);
    }
  }

  if (hasvalue("in")) {
    /* read snapshot and grab nbody.... */

    instr = stropen(getparam("in"),"r");
    get_history(instr);
    if (!get_tag_ok(instr, SnapShotTag)) error("no snapshot");
    get_snap(instr, &btab, &nbody, &tsnap, &bits);
    strclose(instr);
    if (kz[21] == 0) {
      warning("patching kz[22]=1 since in= was given");
      kz[21] = 2;
    }
  } else
    nbody = getiparam("nbody");


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

  q = getdparam("q");

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

  /* there are more optional kz-dependant keywords...
     see nbody4.5, they have not been encoded here
     by lack of Cambridge time 
  */



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

    /*
     * two options here: the snapshot can be written in ascii, and read
     * or stou4 can be used, thereby preserving more digits reliably....
     * in some future version this will be aligned with nbody1 and 2
     * or vice versa....
     */

    if (hasvalue("in")) {
#if 0
      dprintf(0,"Using stou4 to convert data for nbody4\n");
      infile = getparam("in");
      if (*infile == '-') {           /* do something special for pipes */
	sprintf(runcmd,"stou4 %s fort.10",infile);
      } else if (*infile == '/') {    /* regular file */
	sprintf(runcmd,"stou4 %s fort.10",infile);
      } else
	error("in=%s must be an absolute pathname or special file",infile);
      dprintf(0,"%s\n",runcmd);
      if (system(runcmd)) error("Error converting input data");
#else
      dprintf(0,"Using ascii printout with format=%s to convert data for nbody4\n",fmt);
      sprintf(fmt7,"%s %s %s %s %s %s %s\n",fmt,fmt,fmt,fmt,fmt,fmt,fmt);
      outstr = stropen("fort.10","w");
      for (bp=btab; bp<btab+nbody; bp++)
	fprintf(outstr,fmt7,
		Mass(bp),
		Pos(bp)[0],Pos(bp)[1],Pos(bp)[2],
		Vel(bp)[0],Vel(bp)[1],Vel(bp)[2]);
      strclose(outstr);
#endif
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
