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
    "rmax=2\n		Outer cutoff radius",
    "emin=\n            Emin, if Used",
    "emax=\n            Emax, if Used",

    "mass=0\n		Total mass of disk (0 means no masses supplied)",
    "seed=0\n		Usual random number seed",
    "sign=1\n           Sign of Z-angular momentum vector of disk (1=counterclock)",
    "launch=y\n         Launch from Y or X axis",
    "maxlz=f\n          Try and find only the maxlz orbits per energy/radius",
    "headline=\n	Text headline for output",
    "VERSION=0.2\n	31-dec-2019 PJT",
    NULL,
};

string usage="set up a test disk in a spherical potential filling f(E,Lz)";

local real rmin, rmax, mass;
local real pmin, pmax;
local real emin, emax;
local int  jz_sign;
local bool Qangle;
local bool Qenergy;
local bool Qabs;

local int ndisk;
local Body *disk;

local proc potential;


extern double xrandom(double,double), grandom(double,double);

local void writegalaxy(string name, string headline, bool Qmass);
local void testdisk0(int ilaunch);
local void testdisk1(int ilaunch);

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
    bool Qmaxlz = getbparam("maxlz");

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
			      getparam("potpars"),
			      getparam("potfile"));
    pos[0] = pos[1] = pos[2] = 0.0;
    pos[ilaunch] = rmin;
    (*potential)(&ndim,pos,acc,&pmin,&time);
    dprintf(0,"Potential at rmin=%g: %g\n",rmin,pmin);
    pos[ilaunch] = rmax;
    (*potential)(&ndim,pos,acc,&pmax,&time);
    dprintf(0,"Potential at rmax=%g: %g\n",rmax,pmax);


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
    if (Qmaxlz) {
      warning("Special test: sampling only maxlz orbits");
      testdisk0(ilaunch);
    } else
      testdisk1(ilaunch);      
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

void testdisk0(int ilaunch)
{
    Body *dp;
    real rmin2, rmax2, r_i, vcir_i, pot_i, t;
    vector acc_i;
    int i, ndim=NDIM;
    double pos_d[NDIM], acc_d[NDIM], pot_d, time_d = 0.0;



    disk = (Body *) allocate(ndisk * sizeof(Body));
    rmin2 = rmin * rmin;
    rmax2 = rmax * rmax;
    t = 0;    /* dummy time ; we do not support variable time */
    pos_d[0] =  pos_d[1] =  pos_d[2] = 0.0;
    for (dp=disk, i = 0; i < ndisk; dp++, i++) {	/* loop all stars */
	Mass(dp) = mass;
	if (ndisk == 1)
	  r_i = rmin;
	else
	  r_i = sqrt(rmin2 + i * (rmax2 - rmin2) / (ndisk - 1.0));
	pos_d[ilaunch] = r_i;
        (*potential)(&ndim,pos_d,acc_d,&pot_d,&time_d);
        SETV(acc_i,acc_d);
	vcir_i = sqrt(r_i * absv(acc_i));               /* v^2 / r = force */

	Pos(dp)[0]         = pos_d[0];
	Pos(dp)[1]         = pos_d[1];
	Pos(dp)[2]         = pos_d[2];
	Vel(dp)[ilaunch]   = 0.0;
	Vel(dp)[1-ilaunch] = vcir_i * jz_sign * (1 - 2*ilaunch);   // 
	Vel(dp)[2]         = 0.0;
	/* store potential and total energy for debugging */
	Phi(dp) = pot_d;
	Aux(dp) = pot_d + 0.5*(sqr(Vel(dp)[0]) + sqr(Vel(dp)[1]) + sqr(Vel(dp)[2]));
    }
}


void testdisk1(int ilaunch)
{
    Body *dp;
    real r_i, p_j, vcir_i, pot_i, t, v2;
    int  ndisk2, nout;
    vector acc_i;
    
    int i, j, ndim=NDIM;
    double pos_d[NDIM], acc_d[NDIM], pot_d, time_d = 0.0;

    dprintf(0,"Pmin/max: %g %g\n",pmin,pmax);
    ndisk2 = (int) sqrt((double)ndisk);
    dprintf(0,"ndisk by side: %d\n",ndisk2);

    disk = (Body *) allocate(ndisk * sizeof(Body));
    t = 0;   
    pos_d[0] =  pos_d[1] =  pos_d[2] = 0.0;

    nout = 0;
    dp = disk;
    for (i = 0; i < ndisk2; ++i) {	/* loop over all radii */
      dprintf(0,"p_j=%g r_i=%g\n",p_j,r_i);
      r_i = rmin + i * (rmax - rmin) / (ndisk2 - 1.0);
      for (j = 0; j < ndisk2; ++j) {	/* loop over all potentials */
	p_j = pmin + j * (pmax - pmin) / (ndisk2 - 1.0);
	
	pos_d[ilaunch] = r_i;
        (*potential)(&ndim,pos_d,acc_d,&pot_d,&time_d);
	v2 = 2*(p_j - pot_d);
	if (v2 < 0) continue;
	nout++;

	Pos(dp)[0]         = pos_d[0];
	Pos(dp)[1]         = pos_d[1];
	Pos(dp)[2]         = pos_d[2];
	Vel(dp)[ilaunch]   = 0.0;
	Vel(dp)[1-ilaunch] = sqrt(v2) * jz_sign * (1 - 2*ilaunch);   // barbatruuk
	Vel(dp)[2]         = 0.0;
	/* store potential and total energy for debugging */
	Phi(dp) = pot_d;
	Aux(dp) = pot_d + 0.5*(sqr(Vel(dp)[0]) + sqr(Vel(dp)[1]) + sqr(Vel(dp)[2]));
	dp++;
      }
    }
    dprintf(0,"Found %d/%d within CZV\n",nout,ndisk);
    ndisk = nout;
}
