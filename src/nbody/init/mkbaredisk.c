/*
 * MKBAREDISK: set up a (bare) exponential disk.
 *
 *	22-jan-89	V1.1					Josh Barnes
 *	15-nov-90	V1.2 helpvec, set_xrandom, gammas	PJT
 *	24-mar-94	ansi
 *      29-mar-97       V2.0 SINGLEPREC, ndisk= now is nbody=   PJT
 *       9-sep-01       a    gsl/xrandom
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
    "gamma=1.0\n	  fudge factor for v_eff(r) ",
    "epsi=0.025\n	  softening (sets vert. scale) ",
    "seed=0\n	  	  usual random number seed ",
    "headline=\n	  text headline for output ",
    "VERSION=2.0a\n	  9-sep-01 PJT",
    NULL,
};

string usage="set up a (bare) exponential disk";

local real  alpha, rcut, mdisk, Qtoomre, gammas, epsi;
local int   ndisk;
local Body *disktab;


real gdisk(real);

extern double bessi0(double), bessk0(double), bessi1(double), bessk1(double);
extern double xrandom(double,double), grandom(double,double);


nemo_main()
{
    alpha = getdparam("alpha");
    rcut = getdparam("rcut");
    mdisk = getdparam("mdisk");
    Qtoomre = getdparam("Qtoomre");
    gammas = getdparam("gamma");
    ndisk = getiparam("nbody");
    epsi = getdparam("epsi");
    init_xrandom(getparam("seed"));
    inittables();
    makedisk();
    writesnap(getparam("out"), getparam("headline"));
}

#define NTAB 256

local real rcir[NTAB];
local real vcir[4*NTAB];

local real mdsk[NTAB];
local real rdsk[4*NTAB];

inittables()
{

    int i;

    rcir[0] = vcir[0] = 0.0;
    for (i = 1; i < NTAB; i++) {
	rcir[i] = rcut * ((real) i) / (NTAB - 1);
	vcir[i] = sqrt(- gdisk(rcir[i]) * rcir[i]);
    }
    spline(&vcir[NTAB], &rcir[0], &vcir[0], NTAB);
    for (i = 0; i < NTAB; i++) {
	rdsk[i] = rcut * pow(((double) i) / (NTAB - 1), 1.5);
	mdsk[i] = 1 - (1 + alpha * rdsk[i]) * exp(- alpha * rdsk[i]);
    }
    spline(&rdsk[NTAB], &mdsk[0], &rdsk[0], NTAB);
}

real gdisk(real rad)
{
    real x;

    x = 0.5 * alpha * rad;
    return (- alpha*alpha * mdisk * x *
	      (bessi0(x) * bessk0(x) - bessi1(x) * bessk1(x)));
}

makedisk()
{
    Body *bp;
    int i;
    real mdsk_i, rad_i, theta_i, vcir_i, omega, Aoort, kappa;
    real mu, sig_r, sig_t, sig_z, vrad_i, veff_i, vorb_i;

    disktab = (Body *) allocate(ndisk * sizeof(Body));
    for (bp = disktab, i = 0; i < ndisk; bp++, i++) {
	Mass(bp) = mdsk[NTAB-1] / ndisk;
	mdsk_i = mdsk[NTAB-1] * ((real) i + 1.0) / ndisk;
	rad_i = seval(mdsk_i, &mdsk[0], &rdsk[0], &rdsk[NTAB], NTAB);
	theta_i = xrandom(0.0, TWO_PI);
	Phase(bp)[0][0] = rad_i * sin(theta_i);		/* assign positions */
	Phase(bp)[0][1] = rad_i * cos(theta_i);
	Phase(bp)[0][2] = grandom(0.0, 0.5 * epsi);
	vcir_i = seval(rad_i, &rcir[0], &vcir[0], &vcir[NTAB], NTAB);
	omega = vcir_i / rad_i;
	Aoort = - 0.5 *
	    (spldif(rad_i, &rcir[0], &vcir[0], &vcir[NTAB], NTAB) - omega);
	if (omega - Aoort < 0.0)
	    printf("rad_i, omega, Aoort = %f %f %f\n", rad_i, omega, Aoort);
	kappa = 2 * sqrt(omega*omega - Aoort * omega);
	mu = alpha*alpha * mdisk * exp(- alpha * rad_i) / TWO_PI;
	sig_r = 3.358 * Qtoomre * mu / kappa;
	sig_t = 0.5 * sig_r * kappa / omega;
	sig_z = 0.5 * sig_r;
	vrad_i = grandom(0.0, sig_r);
	veff_i = sqrt(vcir_i*vcir_i +
			(gammas - 3*alpha*rad_i) * sig_r*sig_r);
	vorb_i = veff_i + grandom(0.0, sig_t);
	Phase(bp)[1][0] =				/* assign velocities */
	  (vrad_i * Phase(bp)[0][0] + vorb_i * Phase(bp)[0][1]) / rad_i;
	Phase(bp)[1][1] =
	  (vrad_i * Phase(bp)[0][1] - vorb_i * Phase(bp)[0][0]) / rad_i;
	Phase(bp)[1][2] = grandom(0.0, sig_z);
    }
#if 0
    centersnap(disktab, ndisk);
#else
    warning("snapshot has not been centered");
#endif
}

centersnap(btab, nb)
Body *btab;
int nb;
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

writesnap(name, headline)
string name;
string headline;
{
    stream outstr;
    real tzero = 0.0;
    int bits = MassBit | PhaseSpaceBit;

    if (! streq(headline, ""))
	set_headline(headline);
    outstr = stropen(name, "w");
    put_history(outstr);
    put_snap(outstr, &disktab, &ndisk, &tzero, &bits);
    strclose(outstr);
}
