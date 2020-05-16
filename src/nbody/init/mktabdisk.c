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
    "VERSION=0.1\n	16-may-2020 PJT",
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

local real took(real);
local void readtable(string name, bool Qadc);
local void writegalaxy(string name, string headline);
local void testdisk(real mass);
  
local real mysech2(real z)
{
  real y = exp(z);
  return sqr(2.0/(y+1.0/y));
}

local real pick_z(real z0)
{
  real z = frandom(-6.0,6.0,mysech2) * z0;
  return z;
}

local real pick_dv(real r, real z, real z0)
{
#ifdef OLD_BURKERT  
  real dv = 1 - (1 + (z/z0) * tanh(z/z0))*(z0/r);
#else
  //real dv = tanh(z/z0);
  //dv = sqrt(1 - ABS(dv));
  real dv = ABS(z/z0);
  //dv = sqrt(exp(-dv));
  dv = exp(-dv);
#endif  
  dprintf(1,"PJT %g %g %g\n",r,z,dv);
  return dv;
}

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
    bits |= AuxBit;
    bits |= PotentialBit;
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
    real sigma_r, sigma_t, sigma_z;
    vector acc_i;
    int i, ncirc, ndim=NDIM;
    double pos_d[NDIM], acc_d[NDIM], pot_d, time_d = 0.0;
    real totmas1 = 0.0;

    disk = (Body *) allocate(ndisk * sizeof(Body));
    rmin2 = rmin * rmin;
    rmax2 = rmax * rmax;
    theta_i = xrandom(0.0, TWO_PI);
    t = 0;    /* dummy time ; we do not support variable time */
    for (dp=disk, i = 0, ncirc=0; i < ndisk; dp++, i++) {	/* loop all stars */
	Mass(dp) = mass;
	if (ndisk == 1)
	  r_i = rmin;
	else
	  r_i = sqrt(rmin2 + i * (rmax2 - rmin2) / (ndisk - 1.0));
	theta_i = xrandom(0.0, TWO_PI);
        cost = cos(theta_i);
        sint = sin(theta_i);
	Pos(dp)[0] = pos_d[0] = r_i * cost;		/* set positions */
	Pos(dp)[1] = pos_d[1] = r_i * sint;
	Pos(dp)[2] = pos_d[2] = 0.0;                    /* it's a DISK ! */
        //(*potential)(&ndim,pos_d,acc_d,&pot_d,&time_d); /* get forces    */
        SETV(acc_i,acc_d);
	vcir_i = sqrt(r_i * absv(acc_i));               /* v^2 / r = force */
#if 1
	do {                         /* iterate, if needed, to get vrandom */
	    sigma_r = grandom(0.0,frac[0]*vcir_i);
	    sigma_t = grandom(0.0,frac[1]*vcir_i);
	    sigma_z = grandom(0.0,frac[2]*vcir_i);
	    dv_t = sigma_t;
	    dv_r = sigma_r * took(r_i) ;
	    vrandom = sqrt(dv_t*dv_t + dv_r*dv_r);
	    if (vrandom > vcir_i) ncirc++;
	} while (vrandom > vcir_i);
	vcir_i = sqrt((vcir_i-vrandom)*(vcir_i+vrandom));
	dv_r += vrad;
	Vel(dp)[0] =  -vcir_i * sint * jz_sign;
	Vel(dp)[1] =   vcir_i * cost * jz_sign;
	Vel(dp)[0] += cost*dv_r - sint*dv_t;  /* add dispersions */
	Vel(dp)[1] += sint*dv_r + cost*dv_t;
#else
	sigma_r = grandom(0.0,frac[0]*vcir_i);
	sigma_t = grandom(0.0,frac[1]*vcir_i);
	sigma_z = grandom(0.0,frac[2]*vcir_i);
	dv_t = sigma_t;
	dv_r = sigma_r * took(r_i) ;

	/* Qenergy only uses radial motion: thus preserving the 
	 * guiding center for epicycles ?? (Olling 2003)
	 */

	if (Qenergy)
	  vcir_i = sqrt((vcir_i-dv_r)*(vcir_i+dv_r));
	Vel(dp)[0] =  -vcir_i * sint * jz_sign;
	Vel(dp)[1] =   vcir_i * cost * jz_sign;
	Vel(dp)[0] += cost*dv_r;
	Vel(dp)[1] += sint*dv_r;
	if (!Qenergy) {
	  Vel(dp)[0] += -sint*dv_t;
	  Vel(dp)[1] +=  cost*dv_t;
	}

#endif
	Vel(dp)[2] = sigma_z;
	/* store potential and total energy for debugging */
	Phi(dp) = pot_d;
	Aux(dp) = pot_d + 0.5*(sqr(Vel(dp)[0]) + sqr(Vel(dp)[1]) + sqr(Vel(dp)[2]));
    }
    if (ncirc) dprintf(0,"Had to respin random %d times\n",ncirc);
}


/*   
 *   2*omega/kappa; from spline interpolation.....
 */

real took(real r)
{
    return  1.0;
}
