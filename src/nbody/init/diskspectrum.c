/*
 * DISKSPECTRUM.C:      combining mkdisk + snapmass + snaprotate + snapgriof
 *                      to create a faster method for edge_gbt.sh
 *
 *          24-nov-2024   cloned off mkdisk
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <potential.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>

string defv[] = {
    "out=???\n		Output file name (spectrum)",
    "nbody=2048\n	Number of disk particles",
    
    "potname=plummer\n  Name of potential(5)",
    "potpars=\n         Parameters to potential(5); omega needed but not used",
    "potfile=\n         Optional data file with potential(5)",
    
    "rmin=0\n		Inner disk radius",
    "rmax=1\n		Outer cutoff radius",
    "model=0\n		Mass model (0=const 1=exp, 2=f... 3=dac)",
    
    "frac=0\n           Relative vel.disp w.r.t. local rotation speed",
    "seed=0\n		Usual random number seed",

    "angle=f\n          Regular angular distribution?",
    "vrad=0\n           radial velocity",
    "energy=f\n         preserve energy if random motions added?",
    "abs=f\n            Use absolute vel.disp instead of fractional?",
    "z0=0,0\n           Vertical scaleheight for density; use 2nd one for velocity dropoff if needed",
    "vloss=-1\n         Fractional loss of orbital speed at the scaleheight (<1 => Burkert)",

    "inc=30\n           Inclination to observe the disk at",

    "vbeam=5\n          FWHM of smoothing beam in velocity",
    "vrange=400\n       velocity gridding will be -vrange:vrange",
    "nvel=200\n         number of spectral pixels",
    
    "headline=\n	Text headline for output",
    "VERSION=0.1\n	27-nov-2024 PJT",
    NULL,
};

string usage="global spectrum of a rotating thin disk";

local real rmin, rmax, mass;
local bool Qangle;
local bool Qenergy;
local bool Qabs;

local int ndisk;
local real frac[NDIM], vrad;
local Body *disk;

local proc potential;

local real z0_d;    /* the old z0 */
local real z0_v;    /* dropoff in velocity */
local real vloss;
local bool Qrandom;

local real sininc = 1.0;

/* old style */
// #define OLD_BURKERT   1


extern double xrandom(double,double), grandom(double,double);

local real took(real);
local void testdisk(void);
local void spectrum(void);
  
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
    
    potential = get_potential(getparam("potname"),
                    getparam("potpars"), getparam("potfile"));
    rmin = getdparam("rmin");
    rmax = getdparam("rmax");
    vrad = getdparam("vrad");
    ndisk = getiparam("nbody");

    /* z0= is now split in a density and velocity */
    nz = nemoinpr(getparam("z0"),z0,2);
    if (nz == 0)
      z0[0] = z0[1] = 0.0;
    else if (nz == 1)
      z0[1] = 0.0;
    z0_d = z0[0];
    z0_v = z0[1];
  

    vloss = getrparam("vloss");
    if (z0_d > 0) {
      if (vloss < 0)
#ifdef OLD_BURKERT	
	dprintf(0,"Burkert et al 2010 streaming(z) model\n");
#else
        dprintf(0,"new style tanh(z) streaming model\n");
#endif
      else
	dprintf(0,"Toy streaming(z) model with vloss=%g\n",vloss);
    }
    nfrac = nemoinpr(getparam("frac"),frac,NDIM);
    switch (nfrac) {
    case 1:
      frac[1] = frac[0];
      frac[2] = 0.0;
      break;
    case 2:
      frac[2] = 0.0;
      break;
    case 3:
      break;
    default:
      error("%d: bad parsing frac=%s",nfrac,getparam("frac"));
    }
    dprintf(1,"frac: %g %g %g\n",frac[0],frac[1],frac[2]);
    Qrandom = (frac[0]>0 || frac[1]>0 || frac[2]>0);
    if (!Qrandom)
      dprintf(0,"No random motions\n");

    mass = 1.0 / ndisk;
    seed = init_xrandom(getparam("seed"));
    dprintf(1,"Seed=%d\n",seed);
    Qangle = getbparam("angle");
    Qenergy = getbparam("energy");
    Qabs = getbparam("abs");
    testdisk();
    spectrum();
}

/*
 * TESTDISK: use forces due to a potential to make a uniform
 * density test disk.  
 */

void testdisk(void)
{
    Body *dp;
    real rmin2, rmax2, r_i, theta_i, vcir_i;
    real  dv_r, dv_t, sint, cost, vrandom;
    real sigma_r, sigma_t, sigma_z;
    vector acc_i;
    int i, ncirc, ndim=NDIM;
    double pos_d[NDIM], acc_d[NDIM], pot_d, time_d = 0.0;

    disk = (Body *) allocate(ndisk * sizeof(Body));
    rmin2 = rmin * rmin;
    rmax2 = rmax * rmax;
    theta_i = xrandom(0.0, TWO_PI);
    for (dp=disk, i = 0, ncirc=0; i < ndisk; dp++, i++) {	/* loop over all stars */
	Mass(dp) = mass;
	if (ndisk == 1)
	  r_i = rmin;
	else
	  r_i = sqrt(rmin2 + i * (rmax2 - rmin2) / (ndisk - 1.0));
	if (Qangle) {
	  theta_i += TWO_PI/ndisk;
	} else {
	  theta_i = xrandom(0.0, TWO_PI);
	}
        cost = cos(theta_i);
        sint = sin(theta_i);
	Pos(dp)[0] = pos_d[0] = r_i * cost;		/* set positions */
	Pos(dp)[1] = pos_d[1] = r_i * sint;
	Pos(dp)[2] = pos_d[2] = 0.0;                    /* it's a DISK ! */
        (*potential)(&ndim,pos_d,acc_d,&pot_d,&time_d); /* get forces    */
        SETV(acc_i,acc_d);
	vcir_i = sqrt(r_i * absv(acc_i));               /* v^2 / r = force */
	/* now cheat and rotate slower away from the plane */
	if (z0_d > 0) {
	  if (vloss >= 0.0) {
	    // toy model (early sep 2017)
	    Pos(dp)[2] = grandom(0.0, 0.5 * z0_d);
	    vcir_i *= (1-vloss*ABS(Pos(dp)[2]/z0_d));
	    if (vcir_i < 0) vcir_i = 0.0; /* added dec 2017 */
	  } else {
	    // Burkert et al. 2010 model, doesn't need vloss= anymore, triggered when vloss < 0
	    Pos(dp)[2] = pick_z(z0_d);
	    vcir_i = vcir_i*vcir_i;
	    vcir_i *= pick_dv(r_i,Pos(dp)[2],z0_v);
	    if (vcir_i > 0)
	      vcir_i = sqrt(vcir_i);
	    else
	      vcir_i = 0.0;
	  }
	}

#if 1
	if (Qabs) {
	  if (Qrandom) {
	    sigma_r = grandom(0.0,frac[0]);
	    sigma_t = grandom(0.0,frac[1]);
	    sigma_z = grandom(0.0,frac[2]);
	  } else
	    sigma_r = sigma_t = sigma_z = 0.0;
	  Vel(dp)[0] =  -vcir_i * sint ;
	  Vel(dp)[1] =   vcir_i * cost ;
	  Vel(dp)[0] += cost*sigma_r - sint*sigma_t;  /* add dispersions */
	  Vel(dp)[1] += sint*sigma_r + cost*sigma_t;
	} else {
	  do {                         /* iterate, if needed, to get vrandom */
	    if (Qrandom) {
	      sigma_r = grandom(0.0,frac[0]*vcir_i);
	      sigma_t = grandom(0.0,frac[1]*vcir_i);
	      sigma_z = grandom(0.0,frac[2]*vcir_i);
	    } else
	      sigma_r = sigma_t = sigma_z = 0.0;
	    dv_t = sigma_t;
	    dv_r = sigma_r * took(r_i) ;
	    vrandom = sqrt(dv_t*dv_t + dv_r*dv_r);
	    if (vrandom > vcir_i) ncirc++;
	  } while (Qenergy &&  vrandom > vcir_i);
	  vcir_i = sqrt((vcir_i-vrandom)*(vcir_i+vrandom));
	  dv_r += vrad;
	  Vel(dp)[0] =  -vcir_i * sint ;
	  Vel(dp)[1] =   vcir_i * cost ;
	  Vel(dp)[0] += cost*dv_r - sint*dv_t;  /* add dispersions */
	  Vel(dp)[1] += sint*dv_r + cost*dv_t;
	}
#else
	if (Qrandom) {
	  sigma_r = grandom(0.0,frac[0]*vcir_i);
	  sigma_t = grandom(0.0,frac[1]*vcir_i);
	  sigma_z = grandom(0.0,frac[2]*vcir_i);
	} else
	  sigma_r = sigma_t = sigma_z = 0.0;
	dv_t = sigma_t;
	dv_r = sigma_r * took(r_i) ;

	/* Qenergy only uses radial motion: thus preserving the 
	 * guiding center for epicycles ?? (Olling 2003)
	 */

	if (Qenergy)
	  vcir_i = sqrt((vcir_i-dv_r)*(vcir_i+dv_r));
	Vel(dp)[0] =  -vcir_i * sint ;
	Vel(dp)[1] =   vcir_i * cost ;
	Vel(dp)[0] += cost*dv_r;
	Vel(dp)[1] += sint*dv_r;
	if (!Qenergy) {
	  Vel(dp)[0] += -sint*dv_t;
	  Vel(dp)[1] +=  cost*dv_t;
	}

#endif
	Vel(dp)[2] = sigma_z;
	Vel(dp)[1] *=  sininc;
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


void spectrum(void)
{
  double mtot = 0.0;
  Body *dp;
  int i;
  
  for (dp=disk, i = 0; i < ndisk; dp++, i++) {	/* loop over all stars */
    mtot += Mass(dp);
  }
  printf("Total mass: %g   ndisk=%d\n", mtot, ndisk);
}
