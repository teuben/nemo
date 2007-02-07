/*
 * FLOWCODE.C: main routines for 'flow' orbit integrator code
 *	
 *	10-apr-96 V0.1 cloned off potcode
 *      20-jan-04  0.3 diffusion/sigma
 *       3-feb-04  0.4 major CVS version skew fix
 *      24-dec-04  0.6b global fix for MacOSX 
 *       7-feb-07  0.6c prototype fix for gcc4
 *
 * To improve:  use allocate() for number of particles; not static
 *
 * See also:
 *    http://www.amara.com/ftpstuff/streamlines1.txt
 *    http://www.amara.com/ftpstuff/streamlines2.txt
 */

#define global
#include "defs.h"

string defv[] = {	
    "in=???\n		  input file name (snapshot)",
    "out=\n		  output file name (snapshot)",
    "potname=???\n        name of flow_potential(5)",
    "potpars=\n           parameters to flow_potential",
    "potfile=\n           optional filename to flow_potential",
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
    "freqdiff=\n          frequency of diffusion [freq]",
    "seed=0\n		  random seed",
    "headline=\n          random verbiage",
    "VERSION=0.6c\n	  7-feb-07 PJT",
    NULL,
};

string usage = "Analytical flow_potential orbit integrator";

local proc pot;

void nemo_main()
{
    setparams();
    inputdata();
    initoutput();
    initstep(bodytab, nbody, &tnow, force);
    output();
    while (tnow + 0.1/freq < tstop) {
	orbstep(bodytab, nbody, &tnow, force, 1.0/freq, mode);
	output();
    }
    stopoutput();
}

void setparams(void)
{
    infile = getparam("in");
    outfile = getparam("out");

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
    if (hasvalue("freqdiff"))
      freqdiff = getdparam("freqdiff");
    else
      freqdiff = freq;
    if (freqdiff > freq) warning("Cannot use freqdiff > freq");
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
 *              (see orbstep)
 */

void force(bodyptr btab, int nb, real time, bool Qnew)
{
    bodyptr p;
    vector lacc,lpos;
    real   lphi;
    int    ndim=NDIM;
    static real time_next_diff = 0.0;
    static int count=0;

    count += nb;
    dprintf(1,"force %d %d \n",nb,count);

    for (p = btab; p < btab+nb; p++) {		/* loop over bodies */
        SETV(lpos,Pos(p));
        (*pot)(&ndim,lpos,lacc,&lphi,&time);
        Phi(p) = lphi;
        SETV(Acc(p),lacc);
    }
    if (sigma > 0.0) rotate_aux(btab,nb);
    if (sigma == 0.0) return;

    if (time >= time_next_diff && Qnew) {
      diffuse(btab, nb, NDIM, sigma, FALSE);       /* get new diffusion angles */  
      time_next_diff += 1.0/freqdiff;
    }
}
