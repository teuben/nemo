/*
 * FLOWCODE.C: main routines for 'flow' orbit integrator code
 *	
 *	10-apr-96 V0.1 cloned off potcode
 *      20-jan-04  0.3 diffusion/sigma
 *       3-feb-04  0.4 major CVS version skew fix
 *
 * To improve:  use allocate() for number of particles; not static
 *
 * See also:
 *    http://www.amara.com/ftpstuff/streamlines1.txt
 *    http://www.amara.com/ftpstuff/streamlines2.txt
 */

#include "defs.h"

string defv[] = {	
    "in=???\n		  input file name (snapshot)",
    "out=\n		  output file name (snapshot)",
    "potname=???\n        name of flow_potential(5)",
    "potpars=\n           parameters to flow_potential",
    "potfile=\n           optional filename to flow_potential",
    "save=\n		  state file name",
    "freq=64.0\n	  fundamental frequency (inv delta-t)",
    "mode=0\n		  integrator: 0=> Euler 1 => RK, 2 => PC, 3 => PC1 4=>RK4 5=> Leapfrog",
    "tstop=2.0\n	  time to stop integration",
    "freqout=4.0\n	  major data-output frequency",
    "minor_freqout=32.0\n minor data-output frequency",
    "options=\n	          misc options {mass,phi,acc,reset_time,kappa,angle}",
    "eta=0\n		  fractional dissipation per timestep",
    "cell=0.1\n		  cell size",
    "rmax=-1\n            Maximum allowed gridsize -rmax:rmax in all directions",
    "fheat=0\n            diffusion/dissipation ratio",
    "sigma=0\n            diffusion angle (degrees) per timestep",
    "seed=0\n		  random seed",
    "headline=\n          random verbiage",
    "VERSION=0.4\n	  3-feb-04 PJT",
    NULL,
};

string usage = "Analytical flow_potential orbit integrator";

local proc pot;

void nemo_main()
{
    void force(), setparams();

    setparams();
    inputdata();
    initoutput();
    initstep(bodytab, nbody, &tnow, force);
    output();
    while (tnow + 0.1/freq < tstop) {
	orbstep(bodytab, nbody, &tnow, force, 1.0/freq, mode);
#if 0
	/* totally useless here */
	dissipate(bodytab, nbody, NDIM, dr, eta, rmax, fheat);
	diffuse(bodytab, nbody, NDIM, sigma);
#endif
	output();
    }
    stopoutput();
}

void setparams()
{
    infile = getparam("in");
    outfile = getparam("out");
    savefile = getparam("save");

    pot = get_potential (getparam("potname"),
       			 getparam("potpars"), 
			 getparam("potfile"));
    ome = get_pattern();
    ome2 = ome*ome;
    half_ome2 = 0.5 * ome2;
    two_ome = 2.0 * ome;
    freq = getdparam("freq");
    mode = getiparam("mode");
    tstop = getdparam("tstop");
    freqout = getdparam("freqout");
    minor_freqout = getdparam("minor_freqout");
    options = getparam("options");
    eta = getdparam("eta");
    sigma = getdparam("sigma") * PI/180.0;	
    dr[0] = getdparam("cell");
    dr[1] = dr[2] = dr[0];		/* square cells in dissipate */
    rmax = getdparam("rmax");
    fheat = getdparam("fheat");
    headline = getparam("headline");
    set_xrandom(getiparam("seed"));
}

/*
 * FORCE: 'force' calculation routine.
 *		Note: we simply copy the computed flow at the new
 *		position into the force array for further processing later
 */

void force(btab, nb, time)
bodyptr btab;			/* array of bodies */
int nb;				/* number of bodies */
real time;			/* current time */
{
    bodyptr p;
    vector lacc,lpos;
    real   lphi;
    int    ndim=NDIM;


    diffuse(bodytab, nb, NDIM, sigma);

    for (p = btab; p < btab+nb; p++) {		/* loop over bodies */
        SETV(lpos,Pos(p));
        (*pot)(&ndim,lpos,lacc,&lphi,&time);
        Phi(p) = lphi;
        SETV(Acc(p),lacc);
    }
    if (sigma > 0) 
      diffuse(btab, nb, NDIM, sigma);


}
