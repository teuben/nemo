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
 *     24-feb-03      initial  framework for epicycle mode      pjt
 *      5-mar-03 V5.0 analytical epicyclic orbits               pjt
 *      6-may-03     a  missing factor 2 for A ???              pjt
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
    "VERSION=5.0a\n       6-may-03 PJT",
    NULL,
};

string usage = "Analytical potential orbit integrator";

local proc  pot;

void nemo_main()
{
    void force(), force1(), setparams();

    setparams();
    inputdata();
    initoutput();
    if (mode < 0) 
      force1(bodytab, nbody, &tnow);             /* epicycle "integration" constants */
    else
      initstep(bodytab, nbody, &tnow, force);    /* forces from potential */
    output();
    while (tnow + 0.1/freq < tstop) {            /* integration loop */
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


/*
 * FORCE: 'force' calculation routine for epicyclic orbits
 *        where we assume that the particles are:
 *        - planar 2D orbits
 *        - (x0,y0) is the guiding center, (u,v) can be the deviant into the epicycle 
 *          and a (u0,v0) is the deviation from circular motion
 *    *** for now the Z oscillations are not solved for ***
 *
 *   Recall:
 *      A = -1/2 R d(Omega)/dR
 *      Kappa^2 = R d(Omega^2)/dR + 4 Omega^2
 */


void force1(btab, nb, time)
bodyptr btab;			/* array of bodies */
int nb;				/* number of bodies */
real time;			/* current time */
{
  bodyptr p;
  vector lpos,a1,a2;
  real   lphi,r2,r,ome1,ome2, eps, f1,f2,kap2,kap1, A, B, oldkap, tol, tolmin;
  real   vr, vt, vz, lz, eps2, tol2, nu, oldnu;
  int    iter1, iter2, ndim=NDIM;
  
  if (ome!=0.0) error("Force-1 does not work in rotating potentials");

  printf("#   r    omega    kappa      A         B   omega-kappa/2 nu    it it   tol     tol    vr        vt-ome/r   vz      lz\n");

  for (p = btab; p < btab+nb; p++) {		/* loop over bodies */
    r2 = sqr(Pos(p)[0]) + sqr(Pos(p)[1]);
    r = sqrt(r2);
    SETV(lpos,Pos(p));
    (*pot)(&ndim,lpos,a1,&lphi,&time);
    ome2 = sqrt(sqr(a1[0]) + sqr(a1[1]))/r;
    ome1 = sqrt(ome2);
    vr = (Pos(p)[0]*Vel(p)[0] + Pos(p)[1]*Vel(p)[1])/r;
    vt = sqrt( (Vel(p)[0]*Vel(p)[0] + Vel(p)[1]*Vel(p)[1])-vr*vr);
    lz = Pos(p)[0]*Vel(p)[1] - Pos(p)[1]*Vel(p)[0];
    vt -= ome1*r;
    vz = Vel(p)[2];

    /* should do iteration until converging, for now slightly asymmetric results */
    tolmin = 1e-7;   /*  tolerance to achieve */
    eps = 0.01;       /*  first small fractional step in radius */
    for (iter1=0; iter1<15; iter1++) {                  /* iterate to find kappa */
      MULVS(lpos,Pos(p), 1+eps);
      (*pot)(&ndim,lpos,a1,&lphi,&time);
      MULVS(lpos,Pos(p), 1-eps);
      (*pot)(&ndim,lpos,a2,&lphi,&time);
      
      f1 = sqrt(sqr(a1[0]) + sqr(a1[1]))/(r*(1+eps));
      f2 = sqrt(sqr(a2[0]) + sqr(a2[1]))/(r*(1-eps));
      kap2 = (f1-f2)/(2*eps) + 4*ome2;     /* Kappa^2 = R d(Omega^2)/dR + 4 Omega^2 */
      kap1 = sqrt(kap2);
      
      if (iter1==0) {
	oldkap = kap1;
	eps /= 2.0;
	continue;
      } else {
	tol = (kap1-oldkap)/kap1;
	tol = ABS(tol);
	dprintf(1,"%g : iter1 %d %g %g %g\n",r,iter1,eps,kap1,tol);
	if (iter1>1 && tol < tolmin) break;
	eps /= 2.0;
      }
    }
    eps2 = 0.01;
    SETV(lpos,Pos(p));
    for (iter2=0; iter2<10; iter2++) {                  /* iterate to find nu */
      lpos[2] = eps2*r;
      (*pot)(&ndim,lpos,a1,&lphi,&time);
      nu = a1[2]/lpos[2];
      nu = sqrt(ABS(nu));
      if (iter2==0) {
	oldnu = nu;
	eps2 /= 2.0;
	continue;
      } else {
	tol2 = (nu-oldnu)/nu;
	tol2 = ABS(tol2);
	dprintf(1,"%g : iter2 %d %g %g %g\n",r,iter2,eps2,nu,tol2);
	if (iter2>1 && tol2 < tolmin) break;
	eps2 /= 2.0;
      }
    }
    A = -0.25*(f1-f2)/(2*eps)/ome1;     /*  A = -1/2 R d(Omega)/dR   */
#if 1
    A /= 2;  /* ??? */
#endif
    B = A - ome1;
    printf("%g %g %g %g %g %g %g %d %d %g %g  %g %g %g %g\n",
	   r,ome1,kap1,A,B,ome1-kap1/2.0,nu,iter1,iter2,tol,tol2,vr,vt,vz,lz);
    p->A = A;
    p->B = B;
    p->kappa = kap1;
    p->nu = nu;
    p->xiv0 = -vr;
    p->etav0 =  vt;             /* check sign  */
    p->zetav0 =  vz;
    Acc(p)[0] = Pos(p)[0];      /* Guiding Center is stored in the acc's */
    Acc(p)[1] = Pos(p)[1];
    Acc(p)[2] = 0.0;            /* force guiding center in the galactic plane */
    Phi(p) = -1.0;
    Key(p) = p-btab+1;          /* label particles 1....nb */
  }
  printf("### This section of the code is under development, don't believe anything it does\n");
}
