/*
 *  mkhomsph.c: creates an equal-mass homogeneous spherical N-body system
 *.............................................................................
 *    version 1:  May 1989   Piet Hut               email: piet@iassns.bitnet
 *                           Institute for Advanced Study, Princeton, NJ, USA
 *    V1.1   sep-91          Adapted for NEMO2             Peter Teuben
 *    V1.2   jul-91          set_xrandom(seed)             PJT
 *    V1.2a  feb-92          usage,                        PJT
 *    V1.2b  23-mar-97       cleanup protos		   pjt
 *           8-sep-01	     init_xrandom
 *.............................................................................
 *     A homogeneous spherical system cannot be both isotropic and in dynamical
 *  equilibrium.  Therefore several options are presented for providing an
 *  appropriate velocity distribution
 *.............................................................................
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <math.h>

#include <filestruct.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

extern double xrandom(double,double);

/*-----------------------------------------------------------------------------
 *  main  --
 *
 *    options:
 *            The following options are allowed; one for chosing the shape
 *            of the velocity distribution, one for chosing the extent of 
 *            the velocity distribution, and one for the seed of the random
 *            number generator:
 *      shape:
 *            i    isotropic distribution,
 *            t    tangential velocity distribution, in dynamic equilibrium
 *            r    radial velocity distribution.
 *     extent:
 *            b    bound system, but with dynamically unlimited extent
 *            l    bound system, limited dynamically to the original sphere
 *       seed:
 *            #    where # stands for an integer either in the range
 *                  [1,2147483647], or # = 0 which causes a seed to be
 *                  chosen as the value of the UNIX clock in seconds; this
 *                  guarantees that no two calls will give the same value for
 *                  the seed if they are more than 2.0 seconds apart.
 * 
 *               Apart from specifying the seed for the random number generator
 *            at most two of these options can be given, one each from
 *            {i,r,t} and {b,l}.  Every combination is legal, but  -tb  has
 *            the same effect as  -tl , since a tangential distribution is
 *            automatically limited to its original region. 
 *               The default values are: -i -l -s0
 *-----------------------------------------------------------------------------
 */


string defv[] = {
    "out=???\n          Output file (snapshot)",
    "nbody=256\n        Number of bodies",
    "seed=0\n           Random seed",
    "shape=isotropic\n  Velocity shape {isotropic,tangential,radial}",
    "extent=l\n         Boundness",
    "headline=\n        Random extra text",
    "VERSION=1.2b\n     23-mar-97 PJT",
    NULL,
};

string usage = "create an equal-mass homogeneous spherical N-body system";

nemo_main()
{
    int  n, seed, bits;
    real time;
    bool  shape_flag = FALSE;
    bool  extent_flag = FALSE;
    bool  c_flag = FALSE;
    bool  n_flag = FALSE;
    bool  s_flag = FALSE;
    bool  i_flag = FALSE;
    bool  t_flag = FALSE;
    bool  r_flag = FALSE;
    bool  b_flag = FALSE;
    bool  l_flag = FALSE;
    char  *shape, *extent;
    Body *btab;
    stream outstr;

    outstr = stropen(getparam("out"),"w");
    n = getiparam("nbody");
    seed = init_xrandom(getparam("seed"));    
    shape = getparam("shape");  /* default: i */
    switch (*shape) {
        case 'i':   i_flag = TRUE; break;
        case 't':   t_flag = TRUE; break;
        case 'r':   r_flag = TRUE; break;
        default :   error("Illegal shape. Choose one of {iso,rad,tan}");
    }
    extent = getparam("extent");    /* default: l */
    switch (*extent) {
        case 'b':   b_flag = TRUE; break;
        case 'l':   l_flag = TRUE; break;
        default :   error("Illegal extent. Choose one of {b,l}");
    }

    btab = (body *) allocate(n * sizeof(Body));
    mkhomsph(n, btab, i_flag, r_flag, t_flag, b_flag, l_flag);
    time = 0.0;
    bits = TimeBit | MassBit | PhaseSpaceBit;
    if (hasvalue("headline")) set_headline(getparam("headline"));
    put_history(outstr);
    put_snap(outstr,&btab, &n, &time, &bits);
    free(btab);
}

/*-----------------------------------------------------------------------------
 *  mkhomsph  --
 *           note:
 *             the system should be converted to the true c.o.m. coordinates 
 *             at the end of the calculations.
 *-----------------------------------------------------------------------------
 */
mkhomsph(nbody, btab, i_flag, r_flag, t_flag, b_flag, l_flag)
int   nbody;
Body *btab;
bool  i_flag;
bool  t_flag;
bool  r_flag;
bool  b_flag;
bool  l_flag;
{
    int  i;
    real  mass;
    Body *pi;

    mass = 1.0 / (real) nbody;
    
    for (pi = btab; pi < btab + nbody; pi++) {
	Mass(pi) = mass;
        if (i_flag)
	    mkiso_sph(pi, b_flag, l_flag);
	else if (r_flag)
	    mkrad_sph(pi, b_flag, l_flag);
	else if (t_flag)
	    mktang_sph(pi, b_flag, l_flag);
    }

    /* pncom(pn); */        /* C.O.M. is skipped here */
}

/*-----------------------------------------------------------------------------
 *  mkiso_sph  --
 *-----------------------------------------------------------------------------
 */
mkiso_sph(pi, b_flag, l_flag)
Body *pi;
bool  b_flag;
bool  l_flag;
{
    real  r, theta, phi;
    real  v, vmax;
    real  pot, potmax;
    
    r = pow( xrandom(0.0, 1.0), 1.0/3.0);
    theta = acos( xrandom(-1.0, 1.0));
    phi = xrandom(0.0, TWO_PI);
    
    Pos(pi)[0] = r * sin(theta) * cos(phi);
    Pos(pi)[1] = r * sin(theta) * sin(phi);
    Pos(pi)[2] = r * cos(theta);

    pot = -1.5 + 0.5*r*r;
    if (b_flag)
        potmax = 0.0;
    else if (l_flag)
        potmax = -1.0;
    else
        error("mkiso_sph: b_flag and l_flag both FALSE");

    vmax = sqrt(2.0 * (potmax - pot));
    v = pow( xrandom(0.0, vmax*vmax*vmax), 1/3.0);
    theta = acos( xrandom(-1.0, 1.0));
    phi = xrandom(0.0, TWO_PI);
    
    Vel(pi)[0] = v * sin(theta) * cos(phi);
    Vel(pi)[1] = v * sin(theta) * sin(phi);
    Vel(pi)[2] = v * cos(theta);
}

/*-----------------------------------------------------------------------------
 *  mkrad_sph  --
 *-----------------------------------------------------------------------------
 */
mkrad_sph(pi, b_flag, l_flag)
Body  *pi;
bool  b_flag;
bool  l_flag;
{
    real  r, theta, phi;
    real  v, vmax;
    real  pot, potmax;
    bool  inwards;
    
    r = pow( xrandom(0.0, 1.0), 1/3.0);
    theta = acos( xrandom(-1.0, 1.0));
    phi = xrandom(0.0, TWO_PI);
    
    Pos(pi)[0] = r * sin(theta) * cos(phi);
    Pos(pi)[1] = r * sin(theta) * sin(phi);
    Pos(pi)[2] = r * cos(theta);

    pot = -1.5 + 0.5*r*r;
    if (b_flag)
        potmax = 0.0;
    else if (l_flag)
        potmax = -1.0;
    else
        error("mkrad_sph: b_flag and l_flag both FALSE");

    vmax = sqrt(2.0 * (potmax - pot));
    v = pow( xrandom(0.0, vmax*vmax*vmax), 1/3.0);
    if (xrandom(-1.0, 1.0) > 0)
	inwards = TRUE;
    else
	inwards = FALSE;
    
    if (inwards) {
	theta = PI - theta;
	phi += PI ;
    }

    Vel(pi)[0] = v * sin(theta) * cos(phi);
    Vel(pi)[1] = v * sin(theta) * sin(phi);
    Vel(pi)[2] = v * cos(theta);
}

/*-----------------------------------------------------------------------------
 *  mktang_sph  --
 *-----------------------------------------------------------------------------
 */
mktang_sph(pi, b_flag, l_flag)
Body  *pi;
bool  b_flag;
bool  l_flag;
{
    real  r, theta, phi;
    real  v, psi;
    real  sin_psi, cos_psi;
    real  length;
    real  ez[NDIM];
    real  e1[NDIM];
    real  e2[NDIM];
    real  ev[NDIM];
    
    r = pow( xrandom(0.0, 1.0), 1/3.0);
    theta = acos( xrandom(-1.0, 1.0));
    phi = xrandom(0.0, TWO_PI);
    
    Pos(pi)[0] = r * sin(theta) * cos(phi);
    Pos(pi)[1] = r * sin(theta) * sin(phi);
    Pos(pi)[2] = r * cos(theta);

    v = r;
    UNITV(ez, 2);
    CROSSVP(e1, ez, Pos(pi));
    CROSSVP(e2, e1, Pos(pi));
    length = absv(e1);
    DIVVS(e1, e1, length);
    length = absv(e2);
    DIVVS(e2, e2, length);
    psi = xrandom(0.0, TWO_PI);
    cos_psi = cos(psi);
    MULVS(e1, e1, cos_psi);
    sin_psi = sin(psi);
    MULVS(e2, e2, sin_psi);
    ADDV(ev, e1, e2);
    MULVS(Vel(pi), ev, v);
}

/* endof: mkhomsph.c */
