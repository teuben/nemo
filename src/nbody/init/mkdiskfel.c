/*
 * MKDISKFEL.C: set up a test disk in a potential filling f(E,Lz)
 *              try to fill f(E,L) to some degree
 *	
 *       V0.1   - 31-dec-2019   drafted from mkdisk        PJT
 *        0.4   -  1-jan-2020   sampling in E-r space      PJT
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
    "VERSION=0.5\n	2-jan-2020 PJT",
    NULL,
};

string usage="set up a test disk in a spherical potential filling f(E,Lz)";

local real rmin, rmax, mass;
local real pmin, pmax;
local real emin, emax;
local int  jz_sign;

local int ndisk;
local Body *disk;

local proc potential;


extern double xrandom(double,double), grandom(double,double);

local void writegalaxy(string name, string headline, bool Qmass);
local void testdisk0(int ilaunch);
local void testdisk1(int ilaunch);
local void testdisk2(int ilaunch);

void nemo_main()
{
    bool Qmass, Qfixe;
    int seed, ndim=3;
    real pos[3], acc[3], pot;
    real time = 0.0;
    int ilaunch;
    string launch = getparam("launch");
    bool Qmaxlz = getbparam("maxlz");

    if (streq(launch,"x"))
      ilaunch = 0;
    else if (streq(launch,"y"))
      ilaunch = 1;
    else
      error("Launch must be x or y");
    
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

    if (hasvalue("emin")) {
      emin = getrparam("emin");
      emax = getrparam("emax");
      if (emin < emax)
	warning("Different emin/emax not yet supported, assuming emin=%g",emin);
      if (emin < pmin) error("emin=%g  <  pmin=%g",emin,pmin);
      if (emin > pmax) error("emin=%g  >  pmax=%g",emin,pmax);
      Qfixe = TRUE;
    } else
      Qfixe = FALSE;


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
    if (Qmaxlz)
      testdisk0(ilaunch);
    else if (Qfixe)
      testdisk2(ilaunch);
    else
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
    bits |= AccelerationBit;
    bits |= PotentialBit;
    put_snap(outstr, &disk, &ndisk, &tsnap, &bits);
    strclose(outstr);
}

//
//  only 'maxlz' orbits
//

void testdisk0(int ilaunch)
{
    Body *dp;
    real r_i, vcir_i, pot_i, t= 0.0;
    vector acc_i;
    int i, ndim=NDIM;
    double pos_d[NDIM], acc_d[NDIM], pot_d, time_d = 0.0;

    disk = (Body *) allocate(ndisk * sizeof(Body));
    pos_d[0] =  pos_d[1] =  pos_d[2] = 0.0;
    
    for (dp=disk, i = 0; i < ndisk; dp++, i++) {	/* loop all stars */
	Mass(dp) = mass;
	r_i = rmin + i * (rmax - rmin) / (ndisk - 1.0);
	pos_d[ilaunch] = r_i;
        (*potential)(&ndim,pos_d,acc_d,&pot_d,&time_d);
        SETV(acc_i,acc_d);
	vcir_i = sqrt(r_i * ABS(acc_d[ilaunch]));

	Pos(dp)[0]         = pos_d[0];
	Pos(dp)[1]         = pos_d[1];
	Pos(dp)[2]         = pos_d[2];
	Vel(dp)[ilaunch]   = 0.0;
	Vel(dp)[1-ilaunch] = vcir_i * jz_sign * (1 - 2*ilaunch);   // 
	Vel(dp)[2]         = 0.0;
	Phi(dp)            = pot_d;
	Acc(dp)[0]         = acc_d[0];
	Acc(dp)[1]         = acc_d[1];
	Acc(dp)[2]         = acc_d[2];
    }
}

//
//  full sampling of f(E,Lz) but via E,
//

void testdisk1(int ilaunch)
{
    Body *dp;
    real r_i, p_j, v2, t = 0.0;
    int  ndisk2, nout;
    vector acc_i;
    
    int i, j, ndim=NDIM;
    double pos_d[NDIM], acc_d[NDIM], pot_d, time_d = 0.0;

    dprintf(1,"Pmin/max: %g %g\n",pmin,pmax);
    ndisk2 = (int) sqrt((double)ndisk);
    dprintf(0,"Sampling uniformly in E-r 2D ndisk2=%d\n",ndisk2);

    disk = (Body *) allocate(ndisk * sizeof(Body));
    pos_d[0] =  pos_d[1] =  pos_d[2] = 0.0;

    nout = 0;
    dp = disk;
    for (i = 0; i < ndisk2; ++i) {	/* loop over all radii */
      dprintf(1,"p_j=%g r_i=%g\n",p_j,r_i);
      r_i = rmin + i * (rmax - rmin) / (ndisk2 - 1.0);
      for (j = 0; j < ndisk2; ++j) {	/* loop over all potentials */
	p_j = pmin + j * (pmax - pmin) / (ndisk2 - 1.0);
	
	pos_d[ilaunch] = r_i;
        (*potential)(&ndim,pos_d,acc_d,&pot_d,&time_d);
	v2 = 2*(p_j - pot_d);
	if (v2 < 0) continue;     // outside CZV
	nout++;

	Pos(dp)[0]         = pos_d[0];
	Pos(dp)[1]         = pos_d[1];
	Pos(dp)[2]         = pos_d[2];
	Vel(dp)[ilaunch]   = 0.0;
	Vel(dp)[1-ilaunch] = sqrt(v2) * jz_sign * (1 - 2*ilaunch);   // barbatruuk
	Vel(dp)[2]         = 0.0;
	Phi(dp)            = pot_d;
	Acc(dp)[0]         = acc_d[0];
	Acc(dp)[1]         = acc_d[1];
	Acc(dp)[2]         = acc_d[2];
	dp++;
      }
    }
    dprintf(0,"Found %d/%d within CZV\n",nout,ndisk);
    ndisk = nout;
}

//
//  All orbits with given energy (emin)
//

void testdisk2(int ilaunch)
{
    Body *dp;
    real r_i, vcir_i, pot_i, v2, t= 0.0;
    vector acc_i;
    int i, nout, ndim=NDIM;
    double pos_d[NDIM], acc_d[NDIM], pot_d, time_d = 0.0;

    warning("new option");

    disk = (Body *) allocate(ndisk * sizeof(Body));
    pos_d[0] =  pos_d[1] =  pos_d[2] = 0.0;

    nout = 0;
    for (dp=disk, i = 0; i < ndisk; dp++, i++) {	/* loop all stars */
	Mass(dp) = mass;
	r_i = rmin + i * (rmax - rmin) / (ndisk - 1.0);
	pos_d[ilaunch] = r_i;
        (*potential)(&ndim,pos_d,acc_d,&pot_d,&time_d);

	v2 = 2*(emin- pot_d);
	if (v2 < 0) break;     // outside CZV
	nout++;

	Pos(dp)[0]         = pos_d[0];
	Pos(dp)[1]         = pos_d[1];
	Pos(dp)[2]         = pos_d[2];
	Vel(dp)[ilaunch]   = 0.0;
	Vel(dp)[1-ilaunch] = sqrt(v2) * jz_sign * (1 - 2*ilaunch);   // 
	Vel(dp)[2]         = 0.0;
	Phi(dp)            = pot_d;
	Acc(dp)[0]         = acc_d[0];
	Acc(dp)[1]         = acc_d[1];
	Acc(dp)[2]         = acc_d[2];
    }
    dprintf(0,"Found %d/%d within CZV\n",nout,ndisk);
    ndisk = nout;
}
