/*
 * MKFLOWDISK.C: set up a disk with initial conditions taken from 
 *               a potential flow (see flowcode)
 *
 *	
 *	original version: 18-nov-03	Peter Teuben - Maryland
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <potential.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

string defv[] = {
    "out=???\n		  output file name",
    "nbody=2048\n	  number of particles",
    "potname=vrt\n        name of flow potential",
    "potpars=\n           optional parameters to potential",
    "potfile=\n           optional data file with potential",
    "rmin=0.0\n	          inner disk radius",
    "rmax=1.0\n		  outer cutoff radius",
    "mass=0\n             total mass of disk",
    "uniform=t\n          uniform surface density, or use density from potfile?",
    "k=1\n                spiral wavenumber  (> 0 trailing spirals)",
    "phase=0\n            phase offset of spiral at rmax",
    "seed=0\n		  random number seed",
    "nmodel=1\n           number of models",
    "headline=\n	  text headline for output ",
    "VERSION=1.0\n	  18-nov-03 PJT",
    NULL,
};

string usage = "toy spiral density perturbation in a uniform disk";


local real rmin, rmax, sigma[3];
local int  ndisk, nmodel;
local real SPa, SPk, SPw;     /* spiral parameters */
local real width;
local real totmass;

local Body *disk = NULL;

local proc potential;

extern double xrandom(double, double);
extern double grandom(double, double);

nemo_main()
{
    stream outstr;
    
    potential = get_potential(getparam("potname"),
                    getparam("potpars"), getparam("potfile"));
    rmin = getdparam("rmin");
    rmax = getdparam("rmax");
    ndisk = getiparam("nbody");
    nmodel = getiparam("nmodel");
    totmass = getdparam("mass");
    
    SPa = 1.0;   /* getdparam("a") */
    SPk = - getdparam("k");	/* corrected for rot counter clock wise */
    
    SPw = 1.0;  /* SPw = getdparam("w"); */

    init_xrandom(getparam("seed"));

    outstr = stropen(getparam("out"), "w");
    put_history(outstr);
    if (hasvalue("headline"))
	set_headline(getparam("headline"));
    
    while (nmodel--) {
        testdisk(nmodel);
        writegalaxy(outstr);
    }
    strclose(outstr);
    free(disk);
}

/*
 * WRITEGALAXY: write galaxy model to output.
 */

writegalaxy(stream outstr)
{
    real tsnap = 0.0;
    int bits = MassBit | PhaseSpaceBit | TimeBit;

    put_snap(outstr, &disk, &ndisk, &tsnap, &bits);
}

/*
 * TESTDISK: use forces due to a potential to make a uniform
 * density test disk.  
 */

testdisk(int n)
{
    Body *dp;
    real rmin2, rmax2, r_i, theta_i, msph_i, vcir_i, pot_i, mass_i, sigma_i;
    real uni, gau, unifrac, f;
    real cost, sint;
    bool Qsigma = (sigma[0]>0.0 || sigma[1]>0.0);
    int i, ndim=NDIM;
    double pos_d[NDIM], acc_d[NDIM], vel_d[NDIM], pot_d, time_d = 0.0;

    if (disk == NULL) disk = (Body *) allocate(ndisk * sizeof(Body));
    rmin2 = rmin * rmin;
    rmax2 = rmax * rmax;
    mass_i = 1.0 / (real)ndisk;

    for (dp=disk, i = 0; i < ndisk; dp++, i++) {
	Mass(dp) = mass_i;
        r_i = sqrt(rmin2 + xrandom(0.0,1.0) * (rmax2 - rmin2));
        theta_i = grandom(0.0,TWO_PI);
        theta_i -= SPk * r_i * TWO_PI;    /* positive SPk is trailing SP  */
        cost = cos(theta_i);
        sint = sin(theta_i);
	Phase(dp)[0][0] = pos_d[0] = r_i * sint;        /* set positions      */
	Phase(dp)[0][1] = pos_d[1] = r_i * cost;
	Phase(dp)[0][2] = pos_d[2] = 0.0;
        (*potential)(&ndim,pos_d,vel_d,&pot_d,&time_d);      /* get flow    */
	Phase(dp)[1][0] = vel_d[0];
	Phase(dp)[1][1] = vel_d[1];
	Phase(dp)[1][2] = 0.0;
    }
}
