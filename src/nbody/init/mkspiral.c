/*
 * MKSPIRAL.C: set up a uniform-density test disk  with spiral perturbation
 *	in a spherical potential(5) - based on MKDISK
 *	
 *	original version: 13-Mar-89	Peter Teuben - Maryland
 *	 7-apr-90 V1.2 ???					PJT
 *      25-nov-90 V1.3 helpvec etc. 				PJT
 *	20-feb-92 V1.3a usage, nemo_main			PJT
 *	23-feb-93 V1.3b fixed potential() arglist bug		PJT
 *       6            c now writing history			PJT
 *      26-feb-93 V1.4  potname,potpars,potfile official potential keywords
 *	15-sep-95 V1.5  allow many models nmodel=, add sign=	PJT
 *                      radii now fully random fromrmin to rmax
 *	23-mar-97 V1.5b fixed for SINGLEPREC			pjt
 *      27-mar-97 V1.6  add mass= for optional total mass       pjt
 *                      but actually was still broken for SINGLEPREC
 *      10-mar-01 V1.7  allow sigma= to have 3 parameters
 *       8-sep-01       a   init_xrandom
 *       8-apr-03       b   timebit
 *      22-apr-09 V1.8  linear/logarithm spiral option
 *      11-aug-09 V1.9  add sigmator= to control anisotropy     PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <potential.h>
#include <spline.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

string defv[] = {
    "out=???\n		  output file name",
    "nbody=2048\n	  number of particles",
    "potname=plummer\n    name of background potential",
    "potpars=\n           optional parameters to potential",
    "potfile=\n           optional data file with potential",
    "rmin=0.0\n	          inner disk radius",
    "rmax=1.0\n		  outer cutoff radius",
    "mass=0\n             total mass of disk",
    "a=1\n                relative spiral perturbation w.r.t. disk",
    "k=1\n                spiral wavenumber  (> 0 trailing spirals)",
    "w=15\n               gaussian width of spiral - in degrees",
    "angular=t\n	  constant angular width",
    "seed=0\n		  random number seed",
    "sign=1\n             sign of angular momentum of disk",
    "sigma=0\n            velocity dispersion, plus optional central offset and exp scalelength",
    "sigmator=1\n         Ratio sigma_t/sigma_r:  use < 0 to compute it from kappa/2.omega",
    "nmodel=1\n           number of models",
    "headline=\n	  text headline for output ",
    "linear=t\n           Linear or Logarithmic spiral?",
    "VERSION=1.9a\n	  13-jul-2020 PJT",
    NULL,
};

string usage = "toy spiral density perturbation in a uniform disk";

string cvsid="$Id$";


local real rmin, rmax, sigma[3], sigmator;
local int  ndisk, nmodel, nsigma, sign;
local real SPa, SPk, SPw;     /* spiral parameters */
local real width;
local real totmass;
local bool Qlinear;

local Body *disk = NULL;

local proc potential;

extern double xrandom(double, double);
extern double grandom(double, double);

local void inittables(void);
local void writegalaxy(stream outstr);
local void testdisk(int n);

void nemo_main(void)
{
    stream outstr;
    
    potential = get_potential(getparam("potname"),
                    getparam("potpars"), getparam("potfile"));
    rmin = getdparam("rmin");
    rmax = getdparam("rmax");
    ndisk = getiparam("nbody");
    nmodel = getiparam("nmodel");
    sign = getiparam("sign");
    nsigma = nemoinpr(getparam("sigma"),sigma,3);
    sigmator = getrparam("sigmator");
    totmass = getdparam("mass");
    Qlinear = getbparam("linear");
    
    SPa = getdparam("a");
    SPk = - getdparam("k");	/* corrected for rot counter clock wise */
    if (sign < 0) SPk *= -1.0;  /* flip K if spinning the other way around */
    SPw = getdparam("w");
    if (getbparam("angular")) {
	width = -1.0;
	SPw /= 180 * PI;		/* convert degrees to radians */
    } else
        width = SPw;
    init_xrandom(getparam("seed"));
    if (sigmator<1) inittables();

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

#define NTAB 256

local real rcir[NTAB];
local real vcir[4*NTAB];


local void inittables()
{

  int i, ndim=3;
  double pos_d[3], acc_d[3], pot_d, time_d = 0.0;

  pos_d[0] = pos_d[1] = pos_d[2] = 0.0;

  rcir[0] = vcir[0] = 0.0;
  for (i = 1; i < NTAB; i++) {
    pos_d[0] = rcir[i] = rmax * ((real) i) / (NTAB - 1);
    (*potential)(&ndim,pos_d,acc_d,&pot_d,&time_d);      /* get forces    */	
    vcir[i] = -acc_d[0] * rcir[i];
    if (vcir[i] < 0) error(" (inittables) force not inward bound at r=%g",rcir[i]);
    vcir[i] = sqrt(vcir[i]);
  }
  spline(&vcir[NTAB], &rcir[0], &vcir[0], NTAB);
}


/*
 * WRITEGALAXY: write galaxy model to output.
 */

void writegalaxy(stream outstr)
{
    real tsnap = 0.0;
    int bits = MassBit | PhaseSpaceBit | TimeBit;

    put_snap(outstr, &disk, &ndisk, &tsnap, &bits);
}

/*
 * TESTDISK: use forces due to a potential to make a uniform
 * density test disk.  
 */

void testdisk(int n)
{
    Body *dp;
    real rmin2, rmax2, r_i, theta_i, msph_i, vcir_i, pot_i, mass_i, sigma_i;
    real uni, gau, unifrac, f;
    real cost, sint, Aoort, omega, kappa, sigma_t;
    vector acc_i;
    bool Qsigma = (sigma[0]>0.0 || sigma[1]>0.0);
    int i, ndim=NDIM;
    double pos_d[NDIM], acc_d[NDIM], pot_d, time_d = 0.0;
    

    if (disk == NULL) disk = (Body *) allocate(ndisk * sizeof(Body));
    rmin2 = rmin * rmin;
    rmax2 = rmax * rmax;
    mass_i = 1.0 / (real)ndisk;
    uni = PI;                            /* area of uniform */
    if (width<0)                    /* constant angle spiral */
        gau = SPa * SPw * sqrt(TWO_PI);    /* area of gauss - erf??? */
    else                            /* fxied width spiral */
        gau = SPa * atan(2*SPw/(rmin+rmax)) * sqrt(TWO_PI);
    unifrac = uni / (uni+gau);
    dprintf(1,"fraction of particles in disk = %f (%d)\n",unifrac,n);
    if (Qsigma > 0)  
      warning("hot disk, sigma=%g,%g,%g",sigma[0],sigma[1],sigma[2]);

    for (dp=disk, i = 0; i < ndisk; dp++, i++) {
	Mass(dp) = mass_i;
        r_i = sqrt(rmin2 + xrandom(0.0,1.0) * (rmax2 - rmin2));
        f = xrandom(0.0,1.0);			   /* pick random disk/spiral */
	if (f<unifrac)                                     /* disk particle   */
	    theta_i = xrandom(0.0, TWO_PI);	/* pick random angle          */
	else {                                           /* spiral particle   */
            if (width>0)                         /* variable angular width    */
                SPw = atan(width/r_i);              
            theta_i = grandom(0.0,SPw);		/* get gaussian angle         */
            f = xrandom(0.0,1.0);		/* pick which spiral arm      */
            if (f<0.5)                     
                theta_i += PI;
	    if (Qlinear)
	      theta_i -= SPk * r_i * TWO_PI;    /* positive SPk is trailing SP  */
	    else
	      theta_i -= SPk * log(r_i) * TWO_PI;    /* positive SPk is trailing SP  */
        }
        cost = cos(theta_i);
        sint = sin(theta_i);
	Phase(dp)[0][0] = pos_d[0] = r_i * cost;        /* set positions      */
	Phase(dp)[0][1] = pos_d[1] = r_i * sint;
	Phase(dp)[0][2] = pos_d[2] = 0.0;
        (*potential)(&ndim,pos_d,acc_d,&pot_d,&time_d);      /* get forces    */
        SETV(acc_i,acc_d);
	vcir_i = sqrt(r_i * absv(acc_i));
        dprintf(1,"r=%g vcir=%g acc=%g %g\n",r_i,vcir_i,acc_i[0],acc_i[1]);
        if (Qsigma) {
	    sigma_i = sigma[0] + sigma[1]*exp(-r_i/sigma[2]);
	    if (sigmator < 0) {
	      omega = vcir_i / r_i;
	      Aoort = -0.5 * (spldif(r_i, &rcir[0], &vcir[0], &vcir[NTAB], NTAB) - omega);
	      if (omega - Aoort < 0.0)
		dprintf(0,"R, omega, Aoort = %f %f %f\n", r_i, omega, Aoort);
	      kappa = 2 * sqrt(omega*omega - Aoort * omega);
	      sigma_t = 0.5 * sigma_i * kappa / omega;
	    } else
	      sigma_t = sigmator * sigma_i;
	    sigma_i = grandom(0.0,sigma_i);
	    sigma_t = grandom(0.0,sigma_t);
        } else
	  sigma_i = sigma_t = 0.0;
	/* sign > 0 makes for counter clock wise */
	Phase(dp)[1][0] = - sign * (vcir_i+sigma_t) * sint + sigma_i * cost;
	Phase(dp)[1][1] =   sign * (vcir_i+sigma_t) * cost + sigma_i * sint;
	Phase(dp)[1][2] = 0.0;
    }
}
