/*
 * MKFLOWDISK.C: set up a disk with initial conditions taken from 
 *               a potential flow (see also flowcode)
 *
 *   18-nov-03   1.0 cloned off mkspiral - 	Peter Teuben - Maryland
 *   20-nov-03   1.1 changed meaning of angles, and thus signs of k/pitch
 *   21-nov-03   1.2 changed signs once again, make is consistent w/ vrtm51.c
 *   23-nov-03   1.3 add key=
 *   26-nov-03   1.3b   fixed final bugs in sign errors and indexing in binsearch()
 *    3-nov-03   1.3c   implemented uniform=
 *   13-dec-03   1.4 
 *    1-jan-04   1.5 constant=
 *
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
    "uniform=f\n          uniform surface density, or use density from potfile?",
    "k=\n                 spiral wavenumber for linear spirals (> 0 trailing arms)",
    "pitch=\n             pitch angle for log spirals (>0 trailing arms)",
    "phase=0\n            phase offset of spiral at rmax (in degrees)",
    "key=\n               Add a key, if present",
    "seed=0\n		  random number seed",
    "nmodel=1\n           number of models",
    "sign=1\n             Change sign of Z-angular momentum of the disk",
    "test=f\n             test shape of spiral (all particles at 0 phase offset)",
    "constant=f\n         force vt constant in rotating frame ?",
    "headline=\n	  text headline for output ",
    "VERSION=1.5\n	  1-jan-04 PJT",
    NULL,
};

string usage = "toy spiral density perturbation in a uniform disk";


local real rmin, rmax, rref;
local int  ndisk, nmodel, key;
local real SPk;     /* linear spiral wavenumber */
local real pitch;   /* log spiral pitch */
local real totmass;
local real offset;  /* phase offset at rref */
local bool Qtest;
local bool Qlinear;
local bool Quniform;
local bool Qkey;
local bool Qconstant;
local int  jz_sign;

local Body *disk = NULL;
local real theta[361], dens[361];

local proc potential;


void nemo_main()
{
    stream outstr;
    
    potential = get_potential(getparam("potname"),
                    getparam("potpars"), getparam("potfile"));
    rmin = getdparam("rmin");
    rmax = getdparam("rmax");
    rref = rmax;
    if (hasvalue("rref")) rref = getdparam("rref");
    Qkey = hasvalue("key");
    if (Qkey) key = getiparam("key");
    Qconstant = getbparam("constant");

    ndisk = getiparam("nbody");
    jz_sign = getiparam("sign");
    nmodel = getiparam("nmodel");
    totmass = getdparam("mass");
    offset = getdparam("phase") * PI / 180.0;    
    Qtest = getbparam("test");
    Quniform = getbparam("uniform");
    if (ABS(jz_sign) != 1) error("sign must be +1 or -1");

    Qlinear = hasvalue("k");
    if (Qlinear)
      SPk = getdparam("k");	/* corrected for rot counter clock wise */
    else if (hasvalue("pitch"))
      pitch = getdparam("pitch"); /* corrected for rot counter clock wise */
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
    if (Qkey) bits |= KeyBit;

    put_snap(outstr, &disk, &ndisk, &tsnap, &bits);
}

setdensity(void) 
{
  int i, ndim=NDIM;
  double pos_d[NDIM], vel_d[NDIM], den_d, time_d = 0.0;
  double t, tref;
  double rmean = 0.5*(rmin+rmax);
  double tanp = tan(pitch*PI/180.0);

  if (Qlinear)
    tref = 0;
  else
    tref = log(rmean/rref)/tanp - offset;
  dprintf(1,"setdensity - 0:360:1 steps at rref=%g tref=%g\n",rref,tref);
  for (i=0; i<=360; i++) {
    theta[i] = i;
    t = i * PI/180.0;
    pos_d[0] = rmean * cos(t-tref);
    pos_d[1] = rmean * sin(t-tref);
    pos_d[2] = 0.0;
    (*potential)(&ndim,pos_d,vel_d,&den_d,&time_d);
    dens[i] = den_d;
    dprintf(2,"DEN: %g   %g %g  %g %g   %g\n",
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
        return i+1;
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
    real rmin2, rmax2, r_i, theta_i, mass_i, rring;
    real cost, sint;
    int i, ndim=NDIM;
    double pos_d[NDIM], vel_d[NDIM], pot_d, time_d = 0.0;
    double tani = tan(pitch*PI/180.0);

    if (disk == NULL) disk = (Body *) allocate(ndisk * sizeof(Body));
    rmin2 = rmin * rmin;
    rmax2 = rmax * rmax;
    rring = 0.5*(rmin+rmax);
    mass_i = 1.0 / (real)ndisk;

    for (dp=disk, i = 0; i < ndisk; dp++, i++) {
	Mass(dp) = mass_i;
	if (Qkey) Key(dp) = key;
        r_i = sqrt(rmin2 + xrandom(0.0,1.0) * (rmax2 - rmin2));
	if (Qtest)
	  theta_i = 0.0;
	else if (Quniform)
	  theta_i = xrandom(0.0,TWO_PI);
	else 
	  theta_i = frandom(0.0,360.0,density) * PI / 180.0;
	theta_i += offset;
	if (Qlinear) {
	  theta_i -= SPk * (r_i-rref) * TWO_PI;    /* linear spiral arm */
	} else {
	  theta_i -= log(r_i/rref)/tani;           /* logaarithmic spiral arm */
	}
        cost = cos(theta_i);
        sint = sin(theta_i);
	Phase(dp)[0][0] = pos_d[0] = r_i * cost;        /* set positions      */
	Phase(dp)[0][1] = pos_d[1] = r_i * sint;
	Phase(dp)[0][2] = pos_d[2] = 0.0;
        (*potential)(&ndim,pos_d,vel_d,&pot_d,&time_d);      /* get flow    */
	Phase(dp)[1][0] = vel_d[0] * jz_sign;
	Phase(dp)[1][1] = vel_d[1] * jz_sign;
	Phase(dp)[1][2] = 0.0;
    }
}
