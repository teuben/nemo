/*
 * MKVH60: set up a uniform sphere with equal mass particles
 *         according to von Hoerner's (1960) recipe - 1960ZA.....50..184V
 *	
 *         18-apr-2021  drafted                  PJT
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
    "nbody=10\n     Number of particles",
    "seed57=\n      Von Hoerners 1957 seed, must be (0.57 - 0.91)",
    "seed=0\n       NEMO's Random number seed (if used)",
    "zerocm=t\n     Center c.o.m. ?",
    "headline=\n    Text headline for output",
    "VERSION=0.1\n  18-apr-2021 PJT",
    NULL,
};

string usage = "create a uniform sphere of equal mass stars using von Hoerner's 1960 recipe";

string cvsid="$Id$";


local real virial=1, power=0;
local bool zerocm;
local bool exact = TRUE;
local bool Q57;

local Body *btab;
local int nbody;

extern double xrandom(double,double), grandom(double,double);
extern double ran_svh57(double);

void writegalaxy(string name, string headline);
void mksphere(void);
void centersnap(body *btab, int nb);
double get_vr(body *btab, int nb);
double get_ek(Body *btab, int nb);

void nemo_main(void)
{
    int seed;
    real seed57;

    nbody = getiparam("nbody");
    seed = init_xrandom(getparam("seed"));
    zerocm = getbparam("zerocm");
    Q57 = hasvalue("seed57");
    if (Q57) {
      seed57 = getrparam("seed57");
      dprintf(0,"von Hoerner seed57=%g first ran=%g\n",seed57,ran_svh57(seed57));
    }
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
    real r_i, v_i, theta_i, phi_i, mass_i, sigma = 0.0;
    real sinp, cosp, sint, cost;
    real rmax = 1.0;
    real ot = 1.0/3.0;
    int i;

    btab = (Body *) allocate(nbody * sizeof(Body));
    mass_i = 1.0/nbody;
    if (virial > 0.0) {
        /* the sphere has potential energy -3GM^2/(5R)  */
        /* sigma is the 1D velocity dispersion of an isotropic one */
        sigma = sqrt(virial/(5*rmax));
        dprintf(1,"Gaussian isotropic velocities: sigma=%g\n",sigma);
    } 
    for (bp=btab, i = 0; i < nbody; bp++, i++) {
	Mass(bp) = mass_i;

	if (Q57) {
	  error("seed57 method not implemented yet");
	} else {
	  /* positions */
	  r_i = pow(i/(nbody-1.0),ot);
	  theta_i = acos(xrandom(-1.0,1.0));
	  phi_i = xrandom(0.0,TWO_PI);
	  Phase(bp)[0][0] = r_i * sin(theta_i) * cos(phi_i);
	  Phase(bp)[0][1] = r_i * sin(theta_i) * sin(phi_i);
	  Phase(bp)[0][2] = r_i * cos(theta_i);

	  /* velocities, same procedure as positions */
	  r_i = pow(i/(nbody-1.0),ot);
	  theta_i = acos(xrandom(-1.0,1.0));
	  phi_i = xrandom(0.0,TWO_PI);
	  Phase(bp)[1][0] = r_i * sin(theta_i) * cos(phi_i);
	  Phase(bp)[1][1] = r_i * sin(theta_i) * sin(phi_i);
	  Phase(bp)[1][2] = r_i * cos(theta_i);
	}
	
    }
    double rv = get_vr(btab,nbody);
    double v2 = 1/sqrt(2*rv);
    double ek = get_ek(btab,nbody);
    dprintf(0,"virial radius: %g    v2=%g   ek=%g\n",rv,v2,ek);
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

double get_vr(Body *btab, int nb)
{
  real sum;
  Body *bp, *qp;

  sum = 0.0;
  for (bp = btab; bp < btab + nb; bp++) {
    for (qp = bp+1; qp < btab + nb; qp++) {
      sum += 1.0/distv(Phase(bp)[0], Phase(qp)[0]);
    }
  }
  sum *= (2.0) / (nb*(nb-1));
  return 1.0/sum;
}

double get_ek(Body *btab, int nb)
{
  real sum;
  Body *bp;

  sum = 0.0;
  for (bp = btab; bp < btab + nb; bp++)
    sum += dotvp(Phase(bp)[1], Phase(bp)[1]);
  return sum;
}

