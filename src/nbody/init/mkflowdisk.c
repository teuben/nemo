/*
 * MKFLOWDISK.C: set up a disk with initial conditions taken from 
 *               a potential flow (see also flowcode)
 *
 *	original version: 18-nov-03	Peter Teuben - Maryland
 *                        20-nov-03     changed meaning of angles, and thus signs of k/pitch
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <potential.h>
#include <mathfns.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

string defv[] = {
    "out=???\n		  output file name",
    "nbody=2048\n	  number of particles",
    "potname=vrtm51\n     name of flow potential",
    "potpars=\n           optional parameters to potential",
    "potfile=\n           optional data file with potential",
    "rmin=0.0\n	          inner disk radius",
    "rmax=1.0\n		  outer cutoff radius",
    "rref=\n              reference radius for spiral phase (defaults to rmax)",
    "mass=0\n             total mass of disk",
    "uniform=t\n          uniform surface density, or use density from potfile?",
    "k=\n                 spiral wavenumber for linear spirals (> 0 trailing spirals)",
    "pitch=\n             pitch angle for log spirals",
    "phase=0\n            phase offset of spiral at rmax (in degrees)",
    "seed=0\n		  random number seed",
    "nmodel=1\n           number of models",
    "test=f\n             test shape of spiral",
    "headline=\n	  text headline for output ",
    "VERSION=1.1\n	  20-nov-03 PJT",
    NULL,
};

string usage = "toy spiral density perturbation in a uniform disk";


local real rmin, rmax, rref;
local int  ndisk, nmodel;
local real SPk;     /* spiral parameters */
local real pitch;
local real totmass;
local real offset;
local bool Qtest;
local bool Qlinear;

local Body *disk = NULL;
local real theta[361], dens[361];

local proc potential;

extern double xrandom(double, double);
extern double grandom(double, double);

void nemo_main()
{
    stream outstr;
    
    potential = get_potential(getparam("potname"),
                    getparam("potpars"), getparam("potfile"));
    rmin = getdparam("rmin");
    rmax = getdparam("rmax");
    rref = rmax;
    if (hasvalue("rref")) rref = getdparam("rref");

    ndisk = getiparam("nbody");
    nmodel = getiparam("nmodel");
    totmass = getdparam("mass");
    offset = getdparam("phase") * PI / 180.0;    
    Qtest = getbparam("test");


    Qlinear = hasvalue("k");
    if (Qlinear)
      SPk = -getdparam("k");	/* corrected for rot counter clock wise */
    else if (hasvalue("pitch"))
      pitch = -getdparam("pitch"); /* corrected for rot counter clock wise */
    else
      error("Either k= (linear) or pitch= (logarithmic) spiral indicator needed");
    
    init_xrandom(getparam("seed"));

    outstr = stropen(getparam("out"), "w");
    put_history(outstr);
    if (hasvalue("headline"))
	set_headline(getparam("headline"));
    setdensity();
    
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

setdensity(void) 
{
  int i, ndim=NDIM;
  double pos_d[NDIM], vel_d[NDIM], den_d, time_d = 0.0;
  double t;


  dprintf(1,"setdensity - 0:360:1 steps at rmax=%g\n",rmax);
  for (i=0; i<=360; i++) {
    theta[i] = i;
    t = i * PI/180.0;
    pos_d[0] = rmax * cos(t);
    pos_d[1] = rmax * sin(t);
    pos_d[2] = 0.0;
    (*potential)(&ndim,pos_d,vel_d,&den_d,&time_d);
    dens[i] = den_d;
    dprintf(1,"DEN: %g   %g %g  %g %g   %g\n",
	    theta[i],pos_d[0],pos_d[1],vel_d[0],vel_d[1],den_d);
  }
}

/* see also: spline.c(interval) */
 
local int binsearch(real u, real *x, int n)
{
    int i, j, k;
 
    if (u < x[0])                   /* check below left edge */
        return 0;
    else if (u >= x[n-1])           /* and above right edge */
        return n;
    else {
        i = 0;
        k = n;
        while (i+1 < k) {
            j = (i+k)/2;
            if (x[j] <= u)
                i = j;
            else
                k = j;
        }
        return i;
    }
}

double density(double t)
{
  int i;

  i = binsearch(t, theta, 361);
  return dens[i];
}



/*
 * TESTDISK: use forces due to a potential to make a uniform
 * density test disk.  
 */

testdisk(int n)
{
    Body *dp;
    real rmin2, rmax2, r_i, theta_i, mass_i;
    real cost, sint;
    int i, ndim=NDIM;
    double pos_d[NDIM], vel_d[NDIM], pot_d, time_d = 0.0;
    double tani = tan(pitch*PI/180.0);

    if (disk == NULL) disk = (Body *) allocate(ndisk * sizeof(Body));
    rmin2 = rmin * rmin;
    rmax2 = rmax * rmax;
    mass_i = 1.0 / (real)ndisk;

    for (dp=disk, i = 0; i < ndisk; dp++, i++) {
	Mass(dp) = mass_i;
        r_i = sqrt(rmin2 + xrandom(0.0,1.0) * (rmax2 - rmin2));
	if (Qtest)
	  theta_i = 0.0;
	else 
	  theta_i = frandom(0.0,360.0,density) * PI / 180.0;
	theta_i += offset;
	if (Qlinear) {
	  theta_i += SPk * (r_i-rref) * TWO_PI;    /* positive SPk is trailing SP  */
	} else {
	  theta_i += log(r_i/rref)/tani;
	}
        cost = cos(theta_i);
        sint = sin(theta_i);
	Phase(dp)[0][0] = pos_d[0] = r_i * cost;        /* set positions      */
	Phase(dp)[0][1] = pos_d[1] = r_i * sint;
	Phase(dp)[0][2] = pos_d[2] = 0.0;
        (*potential)(&ndim,pos_d,vel_d,&pot_d,&time_d);      /* get flow    */
	Phase(dp)[1][0] = vel_d[0];
	Phase(dp)[1][1] = vel_d[1];
	Phase(dp)[1][2] = 0.0;
    }
}
