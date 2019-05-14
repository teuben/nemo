/*
 *  MAIN.C: NEMO driver for nbody6
 *	
 *  19-feb-2019  V0.1     draft for nbody6++                             PJT
 *                         
 *  CAVEAT:    this is optimized to read a NEMO snapshot and process that and
 *             write a NEMO snapshot.
 *             Although it (should) work with nbody6, it was developed with nbody6++
 * 
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <run.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

string defv[] = {
    "in=\n            input file (optional - see nbody= ) ",
    "outdir=???\n     output run directory (required)",

    "nbody=1000\n     Total particle number (<= NMAX) if no in= given",
    "nfix=1\n         Output frequency of data save or binaries; KZ(3 & 6)",
    "ncrit=5\n        Final particle number (alternative termination criterion)",
    "nrand=123\n      Random number seed",
    "nnbopt=20\n      Desired optimal neighbor number (NNBMAX in nbody6)",
    "nrun=1\n         Run identification index",
    "ncomm=10\n       Frequency to store dumping data (nbody6++ only)",

    "etai=0.02\n      Time-step parameter for irregular force polynomial",
    "etar=0.02\n      Time-step parameter for regular force polynomial",
    "rs0=0.1\n        Initial guess for all radii of neighbor spheres",
    "dtadj=$deltat\n  Time interval for parameter adjustment and energy check",
    "deltat=0.25\n    Output time interval",
    "tcrit=10.0\n     Termination time",
    "qe=2e-5\n        Energy tolerance (restart if DE/E > 5*QE & KZ(2) > 1)",
    "rbar=1.0\n       mean radius of system",
    "zmbar=0.5\n      mean mass of system, in solar units",

    "kz=0 0 1 0 0 1 5 0 0 1  0 0 0 0 2 1 1 0 0 0  1 2 0 0 0 2 0 0 0 2  0 0 2 0 0 0 1 0 0 1  0 0 0 0 0 0 0 0 0 0\n",
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
      "      21  Extra output (>0: model, etc; >1: CENTER; >2: MTRACE; >3: GLOBAL).\n"
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
      "      31  Center of mass correction after energy check.\n"
      "      32  Increase of output intervals (based on single particle energy).\n"
      "      33  Block-step diagnostics at main output (=2: active pipes).\n"
      "    # 34  Roche lobe overflow (suppressed).\n"
      "      35  Time offset (global time from TTOT = TIME + DTOFF; offset = 100).\n"
      "      36  Step reduction for hierarchical systems (not recommended).\n"
      "      37  Step reduction for encounters with high-velocity particles.\n"
      "    # 38  Multiple use of GRAPE-6 (sleep 1 sec after each timer check).\n"
      "      39  Neighbour list (=-1: on host; =0: full list or closest on GRAPE).\n"
      "      40  (not used).\n"
      "    # 41..50  TBD",
    

    "dtmin=1e-4\n     time-step criterion for regularization search",
    "rmin=1e-3\n      distance criterion for regularization search",
    "etau=0.1\n       ",
    "eclose=1.0\n     ",
    "gmin=1e-6\n      ",
    "gmax=0.001\n     ",
    "smax=0.25\n      ",

    "alpha=2.3\n      Power-law index for initial mass function",
    "body1=10.0\n     Maximum particle mass before scaling",
    "bodyn=0.2\n      Minimum particle mass before scaling",
    "nbin0=0\n        number of primordial binaries",
    "nhi0=0\n         number of primordial hierarchical binaries",
    "zmet=0.001\n     Metal abundance (0.03..0.0001)",
    "epoch0=0\n       epoch of star formation (in Myr)",
    "dtplot=10.0\n    plotting interval (nbody units)",

    /* some kz(5) > 0 options here: APO,ECC,N2,SCALE SEMI, M1,M2,ZMH,RCUT */

    "q=0.5\n          Virial ratio (q=0.5 for virial equilibrium)",
    "vxrot=0\n        XY velocity scaling",
    "vzrot=0\n        Z velocity scaling",
    "rtide=1\n        Unscaled tidal radius for kz(14)=2 and kz(22) >= 2",

    /* some more:   GMG,RG0 if kz(14)=2 or 3*/

    /* ah, here they are */

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
    "gpid=0\n         Number of grape boards (starting from 0) ** ignored **",

    "format=%g\n      Format used for fort.10 input conditions if in= used",
    "KZ#=\n           [indexed] Override some kz= keywords",
    "nbody6=0\n       run mode : 1=nbody6  0=nbody6++",
    "exe=\n           Name of the executable (default see nbody6=)",

    "VERSION=0.8b\n   13-may-2019 PJT",
    NULL,
};

string usage="NEMO frontend to the NBODY6(++) code";

string cvsid="$Id$";

#define KZ_MAX  50

void nemo_main(void)
{
  int nbody, nfix, ncrit, nrand, nnbopt, nrun, ncomm, kstart, gpid, nbin0, nhi0, n2;
  real etai, etar, rs0, dtadj, deltat, tcrit, qe, eps, tcomp;
  int k, nkz, kz[KZ_MAX], KZ[KZ_MAX];
  int maxidx, bits;
  real dtmin,rmin,etau,eclose,gmin,gmax,smax,tsnap;
  real alpha, body1, bodyn, zmet, epoch0, dtplot;
  real q, rbar, zmbar;
  real apo, ecc, dmin, semi, scale, m1, m2, ratio, range;
  Body *bp, *btab = NULL;

  string exefile = getparam("exe");
  string parfile = "nbody6.in";
  string rundir = getparam("outdir");
  string fmt = getparam("format");
  int nbody6_mode = getiparam("nbody6");
  string infile;
  char dname[256], runcmd[256], fmt7[256];
  stream datstr, histr, instr, outstr;

  warning("This is a not fully tested pre-release version %s",getparam("VERSION"));
  if (!hasvalue("exe")) {
    exefile = (nbody6_mode == 1) ? "nbody6" : "nbody6++";
    dprintf(0,"exefile=%s\n",exefile);
  }

  kstart = getiparam("kstart");
  tcomp =  getdparam("tcomp");
  gpid = getdparam("gpid");
  
  nfix = getiparam("nfix");
  ncrit = getiparam("ncrit");
  nrand = getiparam("nrand");
  nnbopt = getiparam("nnbopt");
  nrun = getiparam("nrun");
  ncomm = getiparam("ncomm");
  
  etai = getdparam("etai");
  etar = getdparam("etar");
  rs0 = getdparam("rs0");
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

  if (nbody6_mode == 1) {
    if (kz[21] == 2) kz[21] = 3;  /* prevent scaling in nbody6 */
  }

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
      warning("patching kz[22]=2 since in= was given");
      if (nbody6_mode == 0)
	kz[21] = 2;          // nbody6++
      else
	kz[21] = 3;          // nbody6
    }
  } else
    nbody = getiparam("nbody");


  dtmin = getdparam("dtmin");
  rmin = getdparam("rmin");
  etau = getdparam("etau");
  eclose= getdparam("eclose");
  gmin = getdparam("gmin");
  gmax = getdparam("gmax");
  smax = getdparam("smax");  

  
  alpha = getdparam("alpha");
  body1 = getdparam("body1");
  bodyn = getdparam("bodyn");
  nbin0 = getiparam("nbin0");
  nhi0 = getiparam("nhi0");
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



  run_mkdir(rundir);

  sprintf(dname,"%s/%s",rundir,parfile);
  datstr = stropen(dname,"w");    
  
    /*  New Run */
  
  if (kstart == 1) {

    if (nbody6_mode == 0)
      fprintf(datstr,"%d 999999.9 1.E6 40 40 640\n",kstart);
    else
      fprintf(datstr,"%d 999999.9\n",kstart);

    if (nbody6_mode == 0)
      fprintf(datstr,"%d %d %d %d %d %d %d\n",nbody,nfix,ncrit,nrand,nnbopt,nrun,ncomm);
    else
      fprintf(datstr,"%d %d %d %d %d %d\n",nbody,nfix,ncrit,nrand,nnbopt,nrun);

    if (nbody6_mode == 0)
      fprintf(datstr,"%g %g %g %g %g %g %g %g %g\n",etai,etar,rs0,dtadj,deltat,tcrit,qe,rbar,zmbar);
    else
      fprintf(datstr,"%g %g %g %g %g %g %g %g %g\n",etai,etar,rs0,dtadj,deltat,tcrit+dtadj,qe,rbar,zmbar);   // work around tcrit bug in nbody6

    for (k=0; k<KZ_MAX; k++) {
      fprintf(datstr,"%d ",kz[k]);
      if (k%10 == 9)     fprintf(datstr,"\n");
    }
    if (nbody6_mode == 0)
      fprintf(datstr,"%g %g %g %g %g %g %g\n",dtmin,rmin,etau,eclose,gmin,gmax,smax);
    else
      fprintf(datstr,"%g %g %g %g %g %g\n",dtmin,rmin,etau,eclose,gmin,gmax);
    
    fprintf(datstr,"%g %g %g %d %d %g %g %g\n",alpha,body1,bodyn,nbin0,nhi0,zmet,epoch0,dtplot);

    if (kz[4] == 2)
      fprintf(datstr,"%g %g %d %g\n",apo,ecc,n2,scale);
    else if (kz[4] == 3) 
      fprintf(datstr,"%g %g %g %g\n",apo,ecc,dmin,scale);
    else if (kz[4] == 4) 
      fprintf(datstr,"%g %g %g %g\n",semi,ecc,m1,m2);

    if (nbody6_mode == 0)
      fprintf(datstr,"%g 0  0 0 \n",q); 
    else
      fprintf(datstr,"%g 0.0 0.0 0.0 %g\n",q,smax);    // READ (5,*)  Q, VXROT, VZROT, RTIDE, SMAX

    if (kz[7] == 1 || kz[7] == 3)
      fprintf(datstr,"%g %g %g %g 0 0 0 \n",semi,ecc,ratio,range);

    /*  this one is not well documented - gleaned by Massimo from source code */
    fprintf(datstr,"20000.0 2 0\n");
    fprintf(datstr,"20000.0 2 0\n");
    fprintf(datstr,"# this last bogus line should not bother nbody6(++) as it's not supposed to be read\n");
    fprintf(datstr,"# other perhaps useful info: nbody6=%d exefile=%s\n",nbody6_mode,exefile);

    strclose(datstr);

    run_cd(rundir);

    histr = stropen("history","w");
    put_history(histr);
    strclose(histr);

    if (hasvalue("in")) {
      dprintf(1,"Using ascii printout with format=%s to convert data for nbody6\n",fmt);
      sprintf(fmt7,"%s %s %s %s %s %s %s\n",fmt,fmt,fmt,fmt,fmt,fmt,fmt);
      if (nbody6_mode == 0)
	outstr = stropen("dat.10","w");
      else
	outstr = stropen("fort.10","w");
      for (bp=btab; bp<btab+nbody; bp++)
	fprintf(outstr,fmt7,
		Mass(bp),
		Pos(bp)[0],Pos(bp)[1],Pos(bp)[2],
		Vel(bp)[0],Vel(bp)[1],Vel(bp)[2]);
      strclose(outstr);
    }

    sprintf(runcmd,"%s < %s",exefile,parfile);
    run_sh(runcmd);

    if (nbody6_mode == 0)
      sprintf(runcmd,"cat `ls -tr conf.3_*` > OUT3; u3tos OUT3 OUT3.snap mode=6 nbody=%d ; rm OUT3",nbody);
    else
      sprintf(runcmd,"u3tos in=OUT3 out=OUT3.snap mode=1 nbody=%d",nbody);
    run_sh(runcmd);    
    
  } else {
    error("kstart=%d not supported for NBODY6(++) yet; you will need to do it manually",kstart);
    // for nbody6_mode==0
    // kstart=2:   no further lines needed
    //        3:   one line: DTADJ,DELTAT,DTADJ,TNEXT,TCRIT,QE,J,K
    //        4:   one line: ETAI,ETAR,ETAU,DTMIN,RMIN,NCRIT,NNBOPT,SMAX
    //        5:
    // for nbody6_mode==1
    // kstart=2:   standard restart
    //       >2:   read new parameters (see modify.f)
  }
}
