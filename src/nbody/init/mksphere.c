/*
 *  mksphere.c: general spherical distribution with user (table) supplied
 *		density profile
 *
 *	 6-sep-95       V1.0 toy, created for NHK video		pjt
 *       9-sep-01       gsl/xrandom
 */


#include  <stdinc.h>
#include  <getparam.h>
#include  <vectmath.h>
#include  <filestruct.h>

#include <snapshot/snapshot.h>  
#include <snapshot/body.h>
#include <snapshot/put_snap.c>


string  defv[] = {
    "out=???\n		      Output file name",
    "nbody=???\n	      Number of particles",
    "radii=\n                 Radii in sphere",
    "density=\n               Associated densities at given radii",
    "seed=0\n                 Seed for the random number generator",
    "zerocm=t\n               Centrate snapshot (t/f)?",
    "rmin=\n                  minimum radius, if to override from radii",
    "rmax=\n                  maximum radius, if to override from radii",
    "headline=\n	      Verbiage for output",
    "VERSION=1.0b\n           9-sep-01 PJT",
    NULL,
};

string usage="construct an arbitrary spherical mass distribution";

#ifndef MAXTAB
# define MAXTAB  1024
#endif

real rad[MAXTAB], den[MAXTAB], mass[MAXTAB], rmin, rmax;
int ntab;

local double mr(double);
local Body *mkplummer(int , int , bool);

extern double frandom(double, double, rproc);
extern double xrandom(double, double);

nemo_main()
{
    bool    zerocm;
    int     nbody, seed, bits, i, n;
    real    snap_time = 0.0;
    Body    *btab;
    stream  outstr;
    string  headline;
    char    hisline[128];

    nbody = getiparam("nbody");
    seed = init_xrandom(getparam("seed"));
    zerocm = getbparam("zerocm");

    ntab = nemoinpr(getparam("radii"), rad, MAXTAB);
    if (ntab<0) error("Parsing error radii");
    if (ntab<2) error("Need at least two radii");
    rmin = hasvalue("rmin") ? getdparam("rmin") : rad[0];
    rmax = hasvalue("rmax") ? getdparam("rmax") : rad[ntab-1];
    n = nemoinpr(getparam("density"), den, MAXTAB);
    if (n<1) error("Need at least 1 value for density=");
    if (n>ntab) warning("Last %d values of density= ignored",n-ntab);
    if (n<ntab)
        for (i=n; i<ntab; i++)  den[i] = den[i-1];
    make_mass();

    dprintf(0,"seed=%d\n",seed);

    outstr = stropen(getparam("out"), "w");

    btab = mkplummer(nbody, seed, zerocm);
    bits = (MassBit | PhaseSpaceBit | TimeBit);
    sprintf(hisline,"init_xrandom: seed used %d",seed);
    app_history(hisline);
    if (hasvalue("headline")) 
        set_headline(getparam("headline"));
    put_history (outstr);           /* update history */
    put_snap (outstr, &btab, &nbody, &snap_time, &bits);
    strclose(outstr);
}


local Body *mkplummer(int nbody, int seed, bool zerocm)
{
    int  i;
    real  mtot;
    real  radius;		/* absolute value of position vector      */
    real  velocity;		/* absolute value of velocity vector      */
    real  theta, phi;		/* direction angles of above vectors      */
    real  x, y;		        /* for use in rejection technique         */
    real  scalefactor;          /* for converting between different units */
    real  inv_scalefactor;      /* inverse scale factor                   */
    real  sqrt_scalefactor;     /* sqare root of scale factor             */
    real  mrfrac;               /* m( rfrac )                             */
    real  m_min, m_max;         /* mass shell limits for quiet=1          */
    real   m_med;		/* mass shell value for quiet=2           */
    real   w_sum;               /* temporary storage for c.o.m. calc      */
    vector w_pos, w_vel;        /* temporary storage for c.o.m. calc      */
    Body  *btab;                /* pointer to the snapshot                */
    Body  *bp;                  /* pointer to one particle                */

    btab = (Body *) allocate (nbody * sizeof(Body));

    mtot = 0.0;
    for (i = 0, bp=btab; i < nbody; i++, bp++) {
        radius = frandom(rmin, rmax, mr);
	theta = xrandom(-1.0, 1.0);
        theta = acos(theta);
	phi = xrandom(0.0, TWO_PI);
	Pos(bp)[0] = radius * sin( theta ) * cos( phi );
	Pos(bp)[1] = radius * sin( theta ) * sin( phi );
        Pos(bp)[2] = radius * cos( theta );

        Mass(bp) = 1.0/ (real) nbody;

	Vel(bp)[0] = 0.0;
	Vel(bp)[1] = 0.0;
	Vel(bp)[2] = 0.0;
    }

    if (zerocm) {       /* False for Masspectrum */
        w_sum = 0.0;
        CLRV(w_pos);
        CLRV(w_vel);
        for (i = 0, bp = btab; i < nbody; i++, bp++) {
            w_sum = w_sum + 1.0;        /* all bodies same mass */
            ADDV(w_pos, w_pos, Pos(bp));
            ADDV(w_vel, w_vel, Vel(bp));
        }
        DIVVS(w_pos, w_pos, w_sum);
        DIVVS(w_vel, w_vel, w_sum);
        for (i = 0, bp = btab; i < nbody; i++, bp++) {
            SUBV(Pos(bp), Pos(bp), w_pos);
            SUBV(Vel(bp), Vel(bp), w_vel);
        }
    }

    return btab; 
}




make_mass()
{
    int i;

    for (i=0; i<ntab; i++) {
        mass[i] = rad[i] * rad[i] * den[i];
        dprintf(0,"%d: %g %g %g\n",i+1,rad[i],den[i],mass[i]);
        if (i > 0) 
            if (rad[i] < rad[i-1]) error("Radii not in ascending order");
    }        
}

local double mr(double r)
{
    int i, j, k;

    if (r < rad[0] || r > rad[ntab-1]) return 0.0;

    i=0;
    j=ntab-1;
    while (j-i > 1) {
        k = (i+j)/2;
        if (r < rad[k]) 
            j = k;
        else
            i = k;
    }
    if (i==j) {
        warning("mr: Fixing a Programming error");
        j = i+1;
    }
    return mass[i] + (r - rad[i])/(rad[j] - rad[i]) * (mass[j] - mass[i]);
}
