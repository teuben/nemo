/*
 * MKTABDISK.C: create a disk with tabularly given density,rotation and dispersion
 *              solely for simulation observational data cubes and asymmetric drift corrections
 *	
 *       16-mar-2020   Cloned off mkdisk  - Peter Teuben
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <potential.h>
#include <spline.h>
#include <table.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

string defv[] = {
    "in=???\n           Input table (radius,density,velocity,dispersion)",
    "out=???\n		Output file name (snapshot)",
    "nbody=2048\n	Number of disk particles",
    "rmin=\n		Inner disk radius (rmin from table)",
    "rmax=\n		Outer cutoff radius (rmax from table)",
    "mass=\n		Rescale total mass to this?",
    "seed=0\n		Usual random number seed",
    "sign=1\n           Sign of Z-angular momentum vector of disk",
    "adc=t\n            Produce a table of Asymmetric Drift Corrections",
    "headline=\n	Text headline for output",
    "VERSION=0.2\n	16-may-2020 PJT",
    NULL,
};

string usage="set up a (r,d,v,s) table driven test disk";

local real rmin, rmax, mass;
local int  jz_sign;


local int ndisk;
local real frac[NDIM], vrad;
local Body *disk;


local int nrad;
local real *rad, *den, *vel, *sig, *dencoef, *velcoef, *sigcoef;


extern double xrandom(double,double), grandom(double,double);
extern int nemo_file_lines(string,int);

local void readtable(string name, bool Qadc);
local void writegalaxy(string name, string headline);
local void testdisk(real mass);
  
void nemo_main()
{
    int nfrac, seed, nz;
    real z0[2];
    bool Qadc = getbparam("adc");

    readtable(getparam("in"), Qadc);

    rmin = (hasvalue("rmin") ? getdparam("rmin") : rad[0]);
    rmax = (hasvalue("rmax") ? getdparam("rmax") : rad[nrad-1]);
    dprintf(0,"rmin/max = %g %g\n",rmin,rmax);
    
    ndisk = getiparam("nbody");

    jz_sign = getiparam("sign");
    if (ABS(jz_sign) != 1) error("%d: sign must be +1 or -1",jz_sign);

    if (hasvalue("mass"))
      mass = getdparam("mass");
    else
      mass = -1.0;

    seed = init_xrandom(getparam("seed"));
    dprintf(1,"Seed=%d\n",seed);
    testdisk(mass);
    writegalaxy(getparam("out"),getparam("headline"));
}

/* READTABLE:  read table
 *
 */

void readtable(string name, bool Qadc)
{
  int colnr[4];
  real *coldat[4];
  real adcden, adcvel, adcsig;
  real ddendr, dveldr, dsigdr;
  real sig2;
  int i, nmax;
  stream instr;
  
  nmax = nemo_file_lines(name,0);
  if (nmax<=0) error("file_lines returned %d lines in %s\n",nmax,name);
  dprintf(1,"Found nmax=%d\n",nmax);
  
  rad = (real *) allocate(nmax * sizeof(real));
  den = (real *) allocate(nmax * sizeof(real));
  vel = (real *) allocate(nmax * sizeof(real));
  sig = (real *) allocate(nmax * sizeof(real));
  dencoef = (real *) allocate(3*nmax * sizeof(real));
  velcoef = (real *) allocate(3*nmax * sizeof(real));
  sigcoef = (real *) allocate(3*nmax * sizeof(real));

  coldat[0] = rad;        colnr[0] = 1;
  coldat[1] = den;        colnr[1] = 2;
  coldat[2] = vel;        colnr[2] = 3;
  coldat[3] = sig;        colnr[3] = 4;
  instr = stropen(name,"r");
  nrad = get_atable(instr,4,colnr,coldat,nmax);
  dprintf(1,"Found nrad=%d\n",nrad);
  strclose(instr);
  
  spline(dencoef, rad, den, nrad);
  spline(velcoef, rad, vel, nrad);
  spline(sigcoef, rad, sig, nrad);

  if (Qadc) {
    for (i=1; i<nrad; i++) {    // skip first point, assuming it's (0,0)
      ddendr = spldif(rad[i], rad, den, dencoef, nrad);
      dveldr = spldif(rad[i], rad, vel, velcoef, nrad);
      dsigdr = spldif(rad[i], rad, sig, sigcoef, nrad);
      sig2 = sig[i]*sig[i];
      adcden = 1.0 * sig2 * rad[i] * ddendr / den[i];
      adcvel = 0.5 * sig2 * (1 - rad[i] * dveldr / vel[i]);
      adcsig = 2.0 * sig2 * rad[i] * sig[i] * dsigdr;
      printf("%g  %g %g %g  %g %g %g\n",rad[i],den[i],vel[i],sig[i],adcden,adcvel,adcsig);
    }
  }
}


/*
 * WRITEGALAXY: write galaxy model to output.
 */

void writegalaxy(string name, string headline)
{
    stream outstr;
    real tsnap = 0.0;
    int bits;

    if (! streq(headline, ""))
	set_headline(headline);
    outstr = stropen(name, "w");
    put_history(outstr);
    bits = MassBit | PhaseSpaceBit | TimeBit;
    put_snap(outstr, &disk, &ndisk, &tsnap, &bits);
    strclose(outstr);
}

/*
 * TESTDISK: use forces due to a potential to make a uniform
 * density test disk.  
 */

void testdisk(real totmas)
{
    Body *dp;
    real rmin2, rmax2, r_i, theta_i, vcir_i, pot_i, t;
    real  dv_r, dv_t, sint, cost, theta_0, vrandom;
    real den_i, vel_i, sig_i;
    int i, ndim=NDIM;
    double pos_d[NDIM], acc_d[NDIM], pot_d, time_d = 0.0;
    real totmas1 = 0.0;

    disk = (Body *) allocate(ndisk * sizeof(Body));
    rmin2 = rmin * rmin;
    rmax2 = rmax * rmax;
    theta_i = xrandom(0.0, TWO_PI);
    t = 0;    /* dummy time ; we do not support variable time */
    for (dp=disk, i = 0; i < ndisk; dp++, i++) {	/* loop all stars */
	if (ndisk == 1)   // ensure uniform distribution in r^2, but random in angles
	  r_i = rmin;
	else
	  r_i = sqrt(rmin2 + i * (rmax2 - rmin2) / (ndisk - 1.0));
	theta_i = xrandom(0.0, TWO_PI);
        cost = cos(theta_i);
        sint = sin(theta_i);
	
	den_i = seval(r_i, rad, den, dencoef, nrad);		      
	vel_i = seval(r_i, rad, vel, velcoef, nrad);		      
	sig_i = seval(r_i, rad, sig, sigcoef, nrad);

	Mass(dp) = den_i;
	totmas1 += Mass(dp);

	Pos(dp)[0] = r_i * cost;	          	/* set positions */
	Pos(dp)[1] = r_i * sint;
	Pos(dp)[2] = 0.0;                               /* it's a DISK ! */

	Vel(dp)[0] = -vel_i * sint * jz_sign;           // circular orbits
	Vel(dp)[1] =  vel_i * cost * jz_sign;
	Vel(dp)[0] += grandom(0.0, sig_i);              // isotropic vel dispersion
	Vel(dp)[1] += grandom(0.0, sig_i);
	Vel(dp)[2] = 0.0;
    }
    if (totmas > 0) {
      for (dp=disk, i = 0; i < ndisk; dp++, i++)
	Mass(dp) /= totmas1;
    }
}

