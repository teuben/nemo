/*
 * MKNSH96:   set up a Kuzmin disk in a special external potential 
 *
 *         See also Norman, Sellwood & Hasan (1996, ApJ 462, 114)
 *         for a description and usage of this 2D model
 *         If you follow the NSH96 paper, the NEMO potential 'nsh96'
 *         (potname=nsh96) should be used
 *          
 *     3-may-2002   initial writeup             Peter Teuben
 *    15-may-2002
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

#include <potential.h>

#include <math.h>
#include <mathfns.h>
#include <spline.h>

string defv[] = {
    "out=???\n		  output file name ",
    "nbody=1024\n	  number of particles ",
    "potname=nsh96\n      External potential (needs to be nsh96-like)",
    "potpars=\n           Potential parameters (careful: watch mkuz=,akuz=)",
    "potfile=\n           Potential data",
    "mkuz=0.75\n          Mass of Nbody (Kuzmin) disk",
    "akuz=1\n             Length scale of Nbody (Kuzmin) disk",
    "Qtoomre=1.0\n	  Toomre's Q parameter ",
    "seed=0\n      	  usual random number seed ",
    "mode=1\n             Disk mode (0,1)",
    "rmax=6\n             Cutoff radius",
    "time=0.0\n           Time of snapshot",
    "tab=f\n		  table output also? ",
    "zerocm=t\n           center the snapshot?",
    "headline=\n	  text headline for output ",
    "VERSION=0.3\n	  15-may-02 PJT",
    NULL,
};

string usage="set up a Kuzmin disk in an external potential (NSH96)";

local real r_cut, Qtoomre, z0;

local bool Qtab;
local int  ndisk;
local int  cmode;
local int  Kuzmin_only = 1;    /* sigma = M/(2.pi) (1+r^2)^{-3/2} */

local Body *disktab;
local potproc_double mypot;    /* pointer to potential calculator function */

local mkuz, akuz;


local real gdisk(real);
local void inittables(void);
local void makedisk(void);
local void centersnap(Body *btab, int nb);
local void printvec(string name, vector vec);
local void writesnap(string name, string headline);

nemo_main()
{
    int seed;

    warning("Program is under develpment - use at your own risk");

    r_cut = getdparam("rmax");
    Qtoomre = getdparam("Qtoomre");
    ndisk = getiparam("nbody");
    z0 = 0.0;   // or should we getdparam("z0"); ??
    cmode = getiparam("mode");
    seed = init_xrandom(getparam("seed"));
    Qtab = getbparam("tab");
    mypot = get_potential(getparam("potname"), 
			  getparam("potpars"), 
			  getparam("potfile"));

    if (mypot==NULL) error("Potential could not be loaded");

    inittables();
    makedisk();
    writesnap(getparam("out"), getparam("headline"));
}

#define NTAB 256

local real rcir[NTAB];
local real vcir[4*NTAB];

local real mdsk[NTAB];
local real rdsk[4*NTAB];


local void inittables()
{
    int i;

    rcir[0] = vcir[0] = 0.0;
    for (i = 1; i < NTAB; i++) {
	rcir[i] = r_cut * ((real) i) / (NTAB - 1);
	vcir[i] = sqrt(gdisk(rcir[i]) * rcir[i]);
	dprintf(2,"spline rotcur: %g %g\n",rcir[i],vcir[i]);
    }
    spline(&vcir[NTAB], &rcir[0], &vcir[0], NTAB);
}

/* gdisk = axisymmetric force, d(Phi)/d(R) */

local real gdisk(real rad)
{
  double pos[3], acc[3], pot, time;
  int maxdim = 3;
  int probe = 0;   /* along X axis */

  if (kuzmin_only)
    return mkuz*rad/qbe(sqrt(rad*rad+akuz*akuz));
  else {
    pos[0] = pos[1] = pos[2] = 0;
    pos[probe] = rad;
    (*mypot)(&maxdim,pos,acc,&pot,&time);
    return sqrt(-acc[0]*rad);
  }
}


/*
 * Kuzmin disk, a=1 M=1
 */

local void makedisk()
{
  Body *bp;
  int i, nzero=0;
  real mdsk_i, rad_i, theta_i, vcir_i, omega, Aoort, kappa;
  real mu, sig_r, sig_t, sig_z, vrad_i, veff_i, vorb_i;
  real m_cut, m, r1, r2, r3;

  /* warning: if m_cut is not close to 1.0, potential can be wrong */

  m_cut = 1.0 - 1/sqrt(1+r_cut*r_cut);
  dprintf(0,"Rcut=%g Mcut=%g\n",r_cut,m_cut);
  m = 1.0/ndisk;
  
  disktab = (Body *) allocate(ndisk * sizeof(Body));

  for (bp = disktab, i = 0; i < ndisk; bp++, i++) {
    Mass(bp) = m;
    rad_i = sqrt(1/sqr(1-xrandom(0.0,m_cut))-1);
    theta_i = xrandom(0.0, TWO_PI);
    Pos(bp)[0] = rad_i * sin(theta_i);		/* assign positions */
    Pos(bp)[1] = rad_i * cos(theta_i);
    Pos(bp)[2] = 0.0;
    // interpolate
    vcir_i = seval(rad_i, &rcir[0], &vcir[0], &vcir[NTAB], NTAB);
    omega = vcir_i / rad_i;
    Aoort = - 0.5 *
      (spldif(rad_i, &rcir[0], &vcir[0], &vcir[NTAB], NTAB) - omega);
    if (omega - Aoort < 0.0)
      printf("rad_i, omega, Aoort = %f %f %f\n", rad_i, omega, Aoort);
    kappa = 2 * sqrt(omega*omega - Aoort * omega);
    if (Kuzmin_only)
      mu = 1.0/TWO_PI/qbe(sqrt(1+rad_i*rad_i));     /* pure simple Kuzmin disk */
    else
      error("potential not implemented");
    if (cmode==1) {                 /* Straight from Josh - mkbaredisk */
      sig_r = 3.358 * Qtoomre * mu / kappa;
      sig_t = 0.5 * sig_r * kappa / omega;
      sig_z = 0.5 * sig_r;
    } else {
      sig_z = sqrt(PI * mu * z0);          /* isothermal sech sheet */
      sig_r = 2.0 * sig_z;                 /* with constant scaleheight */
      Qtoomre = sig_r * kappa / (3.358 * mu);  /* See vdKruit/Searle */
      sig_t = 0.5 * sig_r * kappa / omega;
    }

    vrad_i = grandom(0.0, sig_r);
#if 0
    /* Asymmetric drift correction */
    if (gammas > 0.0) 			/* Josh' method: averaged */
      veff_i = (vcir_i*vcir_i +
		(gammas - 3*alpha*rad_i) * sig_r*sig_r);
    else				/* a little more accurate */
      veff_i = sqr(vcir_i) - sqr(sig_r) * 
	(sqr(sig_t/sig_r) - 1.5 + 0.5*sqr(sig_z/sig_r) + 2*alpha*rad_i);
    if (veff_i < 0.0) {
      nzero++;
      veff_i = 0.0;
    } else
      veff_i = sqrt(veff_i);
#else
    veff_i = vcir_i;
#endif
    vorb_i = veff_i + grandom(0.0, sig_t);
    Vel(bp)[0] = (vrad_i * Pos(bp)[0] + vorb_i * Pos(bp)[1]) / rad_i;
    Vel(bp)[1] = (vrad_i * Pos(bp)[1] - vorb_i * Pos(bp)[0]) / rad_i;
    Vel(bp)[2] = grandom(0.0, sig_z);
    if (Qtab) {
      if (sig_r == 0.0) {
	r1 = r2 = r3 = 0.0;
      } else {
	r1 =  sig_t/sig_r;
	r2 = sig_z/sig_r;
	r3 = 1.5-(sqr(sig_t) + 0.5*sqr(sig_z))/sqr(sig_r);
      }
      printf ("%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      rad_i,vcir_i,omega,kappa,Aoort,mu,sig_r,sig_t,sig_z,veff_i,
	      Qtoomre,r1,r2,r3);
    }
  }
  if (nzero)  /* warn if something wring with asymmetric drift corrections */
    dprintf(0,"Warning: %d stars with too little orbital motion\n",nzero);
  if (getbparam("zerocm"))
    centersnap(disktab, ndisk);
}

local void centersnap(Body *btab, int nb)
{
    real mtot;
    vector cmphase[2], tmp;
    Body *bp;

    mtot = 0.0;
    CLRV(cmphase[0]);
    CLRV(cmphase[1]);
    for (bp = btab; bp < btab + nb; bp++) {
	mtot = mtot + Mass(bp);
	MULVS(tmp, Phase(bp)[0], Mass(bp));
	ADDV(cmphase[0], cmphase[0], tmp);
	MULVS(tmp, Phase(bp)[1], Mass(bp));
	ADDV(cmphase[1], cmphase[1], tmp);
    }
    MULVS(cmphase[0], cmphase[0], 1.0/mtot);
    MULVS(cmphase[1], cmphase[1], 1.0/mtot);
    for (bp = btab; bp < btab + nb; bp++) {
	SUBV(Phase(bp)[0], Phase(bp)[0], cmphase[0]);
	SUBV(Phase(bp)[1], Phase(bp)[1], cmphase[1]);
    }
    printvec("pos:",cmphase[0]);
    printvec("vel:",cmphase[1]);
}

local void printvec(string name, vector vec)
{
    dprintf(1,"%12s  %10.5f  %10.5f  %10.5f  %10.5f\n",
	   name, absv(vec), vec[0], vec[1], vec[2]);
}


local void writesnap(string name, string headline)
{
    stream outstr;
    real tzero = getdparam("time");
    int bits = MassBit | PhaseSpaceBit | TimeBit;

    if (! streq(headline, ""))
	set_headline(headline);
    outstr = stropen(name, "w");
    put_history(outstr);
    put_snap(outstr, &disktab, &ndisk, &tzero, &bits);
    strclose(outstr);
}
