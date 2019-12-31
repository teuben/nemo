/*
 * MKDISKFEL.C: set up a test disk in a potential filling f(E,Lz)
 *              try to fill f(E,L) to some degree
 *	
 *       V0.1   - 31-dec-2019   drafted from mkdisk        PJT
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
#include <snapshot/put_snap.c>

string defv[] = {
    "out=???\n		Output file name (snapshot)",
    "nbody=2048\n	Number of disk particles ",
    "potname=plummer\n  Name of potential",
    "potpars=\n         Parameters to potential",
    "potfile=\n         Optional data file with potential",
    "rmin=0\n		Inner disk radius",
    "rmax=1\n		Outer cutoff radius",
    "emin=\n            Emin, if Used",
    "emax=\n            Emax, if Used",

    "mass=0\n		Total mass of disk (0 means no masses supplied)",
    "seed=0\n		Usual random number seed",
    "sign=1\n           Sign of Z-angular momentum vector of disk",
    "launch=y\n         Launch from Y or X axis",
    "headline=\n	Text headline for output",
    "VERSION=0.1\n	31-dec-2019 PJT",
    NULL,
};

string usage="set up a test disk in a spherical potential filling f(E,Lz)";

local real rmin, rmax, mass;
local real emin, emax;
local int  jz_sign;
local bool Qangle;
local bool Qenergy;
local bool Qabs;

local int ndisk;
local real frac[NDIM], vrad;
local Body *disk;

local proc potential;


extern double xrandom(double,double), grandom(double,double);

local real took(real);
local void writegalaxy(string name, string headline, bool Qmass);
local void testdisk(void);
  
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
    bool Qmass;
    int nfrac, seed, nz, ndim=3;
    real z0[2];
    real pos[3], acc[3], pot;
    real time = 0.0;
    int ilaunch, jlaunch;
    string launch = getparam("launch");

    if (streq(launch,"x"))
      ilaunch = 0;
    else if (streq(launch,"y"))
      ilaunch = 1;
    else
      error("Launch must be x or y");
    jlaunch = 1 - ilaunch;    // swap X and Y
    
    rmin = getdparam("rmin");
    rmax = getdparam("rmax");
    ndisk = getiparam("nbody");

    potential = get_potential(getparam("potname"),
                    getparam("potpars"), getparam("potfile"));
    pos[0] = pos[1] = pos[2] = 0.0;
    pos[ilaunch] = rmin;
    (*potential)(&ndim,pos,acc,&pot,&time);
    dprintf(0,"Potential at rmin=%g: %g\n",rmin,pot);    
    pos[ilaunch] = rmax;
    (*potential)(&ndim,pos,acc,&pot,&time);
    dprintf(0,"Potential at rmax=%g: %g\n",rmax,pot);    


    jz_sign = getiparam("sign");
    if (ABS(jz_sign) != 1) error("%d: sign must be +1 or -1",jz_sign);

    mass = getdparam("mass") / ndisk;
    if (mass==0.0)  {
	Qmass=FALSE;
	warning("mass=0 -- No masses created in output snapshot");
    } else
        Qmass=TRUE;
    seed = init_xrandom(getparam("seed"));
    dprintf(1,"Seed=%d\n",seed);
    Qangle = FALSE;
    Qenergy = FALSE;
    Qabs = FALSE;
    testdisk();
    writegalaxy(getparam("out"), getparam("headline"), Qmass);
}

/*
 * WRITEGALAXY: write galaxy model to output.
 */

void writegalaxy(string name, string headline, bool Qmass)
{
    stream outstr;
    real tsnap = 0.0;
    int bits;

    if (! streq(headline, ""))
	set_headline(headline);
    outstr = stropen(name, "w");
    put_history(outstr);
    if (Qmass)
        bits = MassBit | PhaseSpaceBit | TimeBit;
    else
        bits = PhaseSpaceBit | TimeBit;
    bits |= AuxBit;
    bits |= PotentialBit;
    put_snap(outstr, &disk, &ndisk, &tsnap, &bits);
    strclose(outstr);
}

/*
 * TESTDISK: use forces due to a potential to make a uniform
 * density test disk.  
 */

void testdisk(void)
{
    Body *dp;
    real rmin2, rmax2, r_i, theta_i, vcir_i, pot_i, t;
    real  dv_r, dv_t, sint, cost, theta_0, vrandom;
    real sigma_r, sigma_t, sigma_z;
    vector acc_i;
    int i, ncirc, ndim=NDIM;
    double pos_d[NDIM], acc_d[NDIM], pot_d, time_d = 0.0;


    

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

#if 1
	if (Qabs) {
	  sigma_r = grandom(0.0,frac[0]);
	  sigma_t = grandom(0.0,frac[1]);
	  sigma_z = grandom(0.0,frac[2]);
	  Vel(dp)[0] =  -vcir_i * sint * jz_sign;
	  Vel(dp)[1] =   vcir_i * cost * jz_sign;
	  Vel(dp)[0] += cost*sigma_r - sint*sigma_t;  /* add dispersions */
	  Vel(dp)[1] += sint*sigma_r + cost*sigma_t;
	} else {
	  do {                         /* iterate, if needed, to get vrandom */
	    sigma_r = grandom(0.0,frac[0]*vcir_i);
	    sigma_t = grandom(0.0,frac[1]*vcir_i);
	    sigma_z = grandom(0.0,frac[2]*vcir_i);
	    dv_t = sigma_t;
	    dv_r = sigma_r * took(r_i) ;
	    vrandom = sqrt(dv_t*dv_t + dv_r*dv_r);
	    if (vrandom > vcir_i) ncirc++;
	  } while (Qenergy &&  vrandom > vcir_i);
	  vcir_i = sqrt((vcir_i-vrandom)*(vcir_i+vrandom));
	  dv_r += vrad;
	  Vel(dp)[0] =  -vcir_i * sint * jz_sign;
	  Vel(dp)[1] =   vcir_i * cost * jz_sign;
	  Vel(dp)[0] += cost*dv_r - sint*dv_t;  /* add dispersions */
	  Vel(dp)[1] += sint*dv_r + cost*dv_t;
	}
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
