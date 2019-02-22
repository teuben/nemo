/*
 * MKHOMSPH.C: set up a uniform sphere
 *	
 *	original version: 20-mar-89	Peter Teuben
 *		apr-89
 *		feb-90	fixed up help strings	PJT
 *		feb-92  usage, set_xrandom      PJT
 *	     23-mar-97  cleanup protos          pjt
 *           20-may-97  added vmax, power for rho=r^{-p}  PJT
 *            9-sep-01  gsl/xrandom
 *           21-mar-04  forgotten initialization  pjt
 *           25-mar-05  time was never written    pjt
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

string defv[] = {	/* DEFAULT INPUT PARAMETERS */
    "out=???\n      Output file name",
    "nbody=2048\n   Number of particles",
    "rmin=0\n       Inner cutoff radius",
    "rmax=1.2\n     Outer cutoff radius",
    "2t/w=0.0\n     Virial ratio (isotropic velocities)",
    "vmax=0\n       alternative (to virial) Vmax for isotropic velocities",
    "power=0\n      Power index for density rho= r^{-p}",
    "seed=0\n       Random number seed",
    "zerocm=t\n     Center c.o.m. ?",
    "headline=\n    Text headline for output",
    "VERSION=1.4c\n 29-aug-2018 PJT",
    NULL,
};

string usage = "create a uniform sphere of equal massive stars";

string cvsid="$Id$";


local real rmin, rmax;
local real virial, power, vmax;
local bool zerocm;

local Body *btab;
local int nbody;

extern double xrandom(double,double), grandom(double,double);

void writegalaxy(string name, string headline);
void mksphere(void);
void centersnap(body *btab, int nb);


void nemo_main(void)
{
    int seed;

    rmin = getdparam("rmin");
    rmax = getdparam("rmax");
    nbody = getiparam("nbody");
    virial = getdparam("2t/w");
    if (virial < 0.0) {
        warning("Virial=%g forced > 0",virial);
        virial = -virial;
    }
    vmax = getdparam("vmax");
    power = getdparam("power");
    if (power >= 3.0) error("Illegal power=%g, must be  < 3.0",power);
    seed = init_xrandom(getparam("seed"));
    zerocm = getbparam("zerocm");
    mksphere();
    writegalaxy(getparam("out"), getparam("headline"));
}

/*
 * WRITEGALAXY: write galaxy model to output.
 */

void writegalaxy(string name, string headline)
{
    stream outstr;
    real tsnap = 0.0;
    int bits = MassBit | PhaseSpaceBit | TimeBit;

    if (! streq(headline, ""))
	set_headline(headline);
    outstr = stropen(name, "w");
    put_history(outstr);
    put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    strclose(outstr);
}

/*
 * MKSPHERE: homogeneous sphere
 */

void mksphere(void)
{
    Body *bp;
    real rmin3, rmax3, r_i, v_i, theta_i, phi_i, mass_i, sigma = 0.0;
    real sinp, cosp, sint, cost;
    int i;

    btab = (Body *) allocate(nbody * sizeof(Body));
#ifdef OLD
    rmin3 = qbe(rmin);
    rmax3 = qbe(rmax);
#else
    rmin3 = pow(rmin,3-power);
    rmax3 = pow(rmax,3-power);
#endif
    mass_i = 1.0/nbody;
    if (virial > 0.0) {
        /* the sphere has potential energy -3GM^2/(5R)  */
        /* sigma is the 1D velocity dispersion of an isotropic one */
        sigma = sqrt(virial/(5*rmax));		/* rmin==0  and power==0 !! */
        dprintf(1,"Gaussian isotropic velocities: sigma=%g\n",sigma);
    } else if (vmax > 0.0) {
        dprintf(1,"Vmax isotropic velocities: vmax=%g\n",vmax);
    } 
    for (bp=btab, i = 0; i < nbody; bp++, i++) {
	Mass(bp) = mass_i;
#ifdef OLD
        r_i = pow ( rmin3 + i * (rmax3-rmin3)/(nbody-1),1.0/3.0);
#else
        r_i = pow ( rmin3 + xrandom(0.0,1.0)*(rmax3-rmin3), 1.0/(3-power));
#endif
        theta_i = acos(xrandom(-1.0,1.0));
        phi_i = xrandom(0.0,TWO_PI);
	Phase(bp)[0][0] = r_i * sin(theta_i) * cos(phi_i);
	Phase(bp)[0][1] = r_i * sin(theta_i) * sin(phi_i);
	Phase(bp)[0][2] = r_i * cos(theta_i);
        if (vmax > 0) {
            v_i = pow( xrandom(0.0, vmax*vmax*vmax), 1/3.0);
            theta_i = acos(xrandom(-1.0,1.0));
            phi_i = xrandom(0.0,TWO_PI);
	    Phase(bp)[1][0] = v_i * sin(theta_i) * cos(phi_i);
    	    Phase(bp)[1][1] = v_i * sin(theta_i) * sin(phi_i);
	    Phase(bp)[1][2] = v_i * cos(theta_i);
        } else {
    	    Phase(bp)[1][0] = grandom(0.0,sigma);
	    Phase(bp)[1][1] = grandom(0.0,sigma);
	    Phase(bp)[1][2] = grandom(0.0,sigma);
        }
    }
    if (zerocm)
        centersnap(btab,nbody);
}

void centersnap(Body *btab, int nb)
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

