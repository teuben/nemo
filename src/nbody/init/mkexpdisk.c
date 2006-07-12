/*
 * MKEXPDISK:   set up a (bare) exponential disk. (based on mkbaredisk)
 *		has table output for debugging
 *		more mode creation methods
 *	original	created		Josh Barnes
 *	15-nov-90	helpvec		PJT
 *	24-mar-94	ansi fix
 *	24-mar-97	proto fixes	pjt
 *	29-mar-97	SINGLEPREC fixed, ndisk= now nbody=	pjt
 *      29-may-01       Add time      PJT
 *       8-sep-01       gsl/xrandom
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

#include <spline.h>

string defv[] = {	/* DEFAULT INPUT PARAMETERS */
    "out=???\n		  output file name ",
    "nbody=1024\n	  number of particles ",
    "alpha=4.0\n	  inverse exponential scale length ",
    "rcut=1.25\n	  outer cutoff radius ",
    "mdisk=1.0\n	  disk mass ",
    "Qtoomre=1.0\n	  Toomre's Q parameter ",
    "gamma=1.0\n	  fudge factor for v_eff(r) if > 0 ",
    "z0=0.025\n	          vertical scaleheight (softening) ",
    "seed=12345\n	  usual random number seed ",
    "mode=1\n             creation mode: 1=josh 2=kruit/searle ",
    "time=0.0\n           tag a time, normally skipped",
    "tab=f\n		  table output also? ",
    "zerocm=t\n           center the snapshot?",
    "headline=\n	  text headline for output ",
    "VERSION=1.2d\n	  7-nov-05 PJT",
    NULL,
};

string usage="set up a (bare) exponential disk";

local real alpha, rcut, mdisk, Qtoomre, gammas, z0;

local bool Qtab;
local int  ndisk;
local int  cmode;

local Body *disktab;


local real gdisk(real);
local void inittables(void);
local void makedisk(void);
local void centersnap(Body *btab, int nb);
local void writesnap(string name, string headline);

extern double bessi0(double), bessk0(double), bessi1(double), bessk1(double);
extern double xrandom(double,double), grandom(double,double);


nemo_main()
{
    int seed;

    alpha = getdparam("alpha");
    rcut = getdparam("rcut");
    mdisk = getdparam("mdisk");
    Qtoomre = getdparam("Qtoomre");
    gammas = getdparam("gamma");
    ndisk = getiparam("nbody");
    z0 = getdparam("z0");
    cmode = getiparam("mode");
    seed = init_xrandom(getparam("seed"));
    Qtab = getbparam("tab");
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
	rcir[i] = rcut * ((real) i) / (NTAB - 1);
	vcir[i] = sqrt(- gdisk(rcir[i]) * rcir[i]);
    }
    spline(&vcir[NTAB], &rcir[0], &vcir[0], NTAB);
    for (i = 0; i < NTAB; i++) {
	rdsk[i] = rcut * pow(((real) i) / (NTAB - 1), 1.5);
	mdsk[i] = 1 - (1 + alpha * rdsk[i]) * exp(- alpha * rdsk[i]);
    }
    spline(&rdsk[NTAB], &mdsk[0], &rdsk[0], NTAB);
}

local real gdisk(real rad)
{
    real x;

    x = 0.5 * alpha * rad;
    return - alpha*alpha * mdisk * x *
	      (bessi0(x) * bessk0(x) - bessi1(x) * bessk1(x));
}

local void makedisk()
{
    Body *bp;
    int i, nzero=0;
    real mdsk_i, rad_i, theta_i, vcir_i, omega, Aoort, kappa;
    real mu, sig_r, sig_t, sig_z, vrad_i, veff_i, vorb_i;

    disktab = (Body *) allocate(ndisk * sizeof(Body));
    for (bp = disktab, i = 0; i < ndisk; bp++, i++) {
	Mass(bp) = mdsk[NTAB-1] / ndisk;
	mdsk_i = mdsk[NTAB-1] * ((real) i + 1.0) / ndisk;
	rad_i = seval(mdsk_i, &mdsk[0], &rdsk[0], &rdsk[NTAB], NTAB);
	theta_i = xrandom(0.0, TWO_PI);
	Pos(bp)[0] = rad_i * sin(theta_i);		/* assign positions */
	Pos(bp)[1] = rad_i * cos(theta_i);
	Pos(bp)[2] = grandom(0.0, 0.5 * z0);
	vcir_i = seval(rad_i, &rcir[0], &vcir[0], &vcir[NTAB], NTAB);
	omega = vcir_i / rad_i;
	Aoort = - 0.5 *
	    (spldif(rad_i, &rcir[0], &vcir[0], &vcir[NTAB], NTAB) - omega);
	if (omega - Aoort < 0.0)
	    printf("rad_i, omega, Aoort = %f %f %f\n", rad_i, omega, Aoort);
	kappa = 2 * sqrt(omega*omega - Aoort * omega);
	mu = alpha*alpha * mdisk * exp(- alpha * rad_i) / TWO_PI;
	if (cmode==1) {                 /* Straight from Josh - mkbaredisk*/
	   sig_r = 3.358 * Qtoomre * mu / kappa;
	   sig_t = 0.5 * sig_r * kappa / omega;
	   sig_z = 0.5 * sig_r;
	} else if (cmode==2) {
	   sig_z = sqrt(PI * mu * z0);          /* isothermal sech sheet */
           sig_r = 2.0 * sig_z;                 /* with constant scaleheight */
           Qtoomre = sig_r * kappa / (3.358 * mu);  /* See vdKruit/Searle */
	   sig_t = 0.5 * sig_r * kappa / omega;
        } else
	    error("illegal mode=%d",cmode);

	vrad_i = grandom(0.0, sig_r);
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
	vorb_i = veff_i + grandom(0.0, sig_t);
	Vel(bp)[0] =				/* assign velocities */
	  (vrad_i * Pos(bp)[0] + vorb_i * Pos(bp)[1]) / rad_i;
	Vel(bp)[1] =
	  (vrad_i * Pos(bp)[1] - vorb_i * Pos(bp)[0]) / rad_i;
	Vel(bp)[2] = grandom(0.0, sig_z);
	if (Qtab) {
            printf ("%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
            rad_i,mdsk_i,vcir_i,omega,kappa,Aoort,mu,sig_r,sig_t,sig_z,veff_i,
            Qtoomre,
            sig_t/sig_r,sig_z/sig_r,
            1.5-(sqr(sig_t) + 0.5*sqr(sig_z))/sqr(sig_r) );
        }
    }
    if (nzero)
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

