/*
 * POTCODE.C: main routines for fixed potential orbit integrator code
 *
 *	 feb-89	V1.0 Original code				PJT
 *	 apr-90	V2.0 new potential(5), renamed some keywords	PJT
 *    21-nov-90 V2.1 void/int repair				PJT
 *	1-nov-91 V2.2 allow absent masses			PJT
 *	6-jun-92 V3.0 added rotating potentials                 PJT
 *     17-jun-92 V3.1 allow energy conserved in 'dissipate()'   PJT
 *     19-jun-92 V3.2 added sigma= diffusion term               PJT
 *                    renamed internal variable eps to eta 
 *     21-jun-92 V3.2a  fixed the integrator error when faking diss/diff   PJT
 *     29-sep-92 V4.0 allow max grid size                       PJT
 *		      options= fix in code_io.c
 *     16-feb-03 V4.0a use get_pattern()			pjt
 *
 * To improve:  use allocate() for number of particles; not static
 */

#include "defs.h"

string defv[] = {	
    "in=???\n		  input file name (snapshot)",
    "out=\n		  output file name (snapshot)",
    "potname=???\n        name of potential(5)",
    "potpars=\n           parameters to potential",
    "potfile=\n           optional filename to potential",
    "save=\n		  state file name",
    "freq=64.0\n	  fundamental frequency (inv delta-t)",
    "mode=4\n		  integrator: 0=> Euler 1 => RK, 2 => PC, 3 => PC1 4=>RK4 5=> Leapfrog",
    "tstop=2.0\n	  time to stop integration",
    "freqout=4.0\n	  major data-output frequency",
    "minor_freqout=32.0\n minor data-output frequency",
    "options=\n		  misc options {mass,phi,acc,reset_time}",
    "eta=0\n		  fractional dissipation per timestep",
    "cell=0.1\n		  cell size",
    "rmax=-1\n            Maximum allowed gridsize -rmax:rmax in all directions",
    "sigma=0\n            diffusion angle (degrees) per timestep",
    "seed=0\n		  random seed",
    "headline=PotCode\n   random mumble for humans",
    "VERSION=4.0a\n	  15-feb-03 PJT",
    NULL,
};

string usage = "Analytical potential orbit integrator";

local proc  pot;

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
	dissipate(bodytab, nbody, NDIM, dr, eta, rmax);
	diffuse(bodytab, nbody, NDIM, sigma);
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
    ome = get_pattern();     /* pattern speed first par of potential */
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
    headline = getparam("headline");
    set_xrandom(getiparam("seed"));
}

/*
 * FORCE: force calculation routine.
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

    for (p = btab; p < btab+nb; p++) {		/* loop over bodies */
        SETV(lpos,Pos(p));
        (*pot)(&ndim,lpos,lacc,&lphi,&time);

	if (ome!=0.0) {
	   lphi -= half_ome2*(sqr(lpos[0])+sqr(lpos[1]));
           lacc[0] += ome2*lpos[0] + two_ome*Vel(p)[1];
           lacc[1] += ome2*lpos[1] - two_ome*Vel(p)[0];
        }

        Phi(p) = lphi;
        SETV(Acc(p),lacc);
    }
}
