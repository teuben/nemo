/*
 * SNAPTRANS.C: coordinate transformations
 *
 *       2-nov-94   V1.0    Created             Peter Teuben
 *	27-may-95   V1.0a   fixed bug in cyl vel. conversion	pjt
 *	
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {
    "in=???\n	    Input snapshot filename",
    "out=???\n	    Output snapshot filename",
    "ctypei=\n      Input coordinate system (cart,sph,cyl)",
    "ctypeo=\n      Output coordinate system (cart,sph,cyl)",
    "projpi=\n      Input projection parameters, if needed",
    "projpo=\n      Output projection parameters, if needed",
#if  0
    /* some day, we'll do real Celestial Coordinates here too   */
    /* hence the funny FITS compatible CTYPE/PROJP names        */
    "equinox=\n     Epoch, if to override **not used**",
    "radecsys=\n    if to override **not used**",
#endif    
    "scale=phase,acc\n  Items to scale (phase,pos,vel,acc) **not used**",
    "ndim=\n        Force dimensionality, if it can be overriden",
    "times=all\n    Times to select snapshots from",
    "VERSION=1.0a\n 27-may-95 PJT",
    NULL,
};

string usage="coordinate transformations of a snapshot";

#define TIMEFUZZ	0.000001	/* tolerance in time comparisons */

extern void cartesian_spherical(real *, real *, real *);
extern void spherical_cartesian(real *, real *, real *);
extern void cartesian_cylindrical(real *, real *, real *);
extern void cylindrical_cartesian(real *, real *, real *);

void nemo_main()
{
    stream instr, outstr;
    string times, ctypei, ctypeo;
    real   tsnap;
    int i, nbody, bits, mode;
    Body *btab = NULL, *bp;
    proc transform, trans2;

    times = getparam("times");
    ctypei = getparam("ctypei");
    ctypeo = getparam("ctypeo");
    if (streq(ctypei,"cart") && streq(ctypeo,"sph"))
        transform = cartesian_spherical;
    else if (streq(ctypei,"sph") && streq(ctypeo,"cart"))
        transform = spherical_cartesian;
    else if (streq(ctypei,"cart") && streq(ctypeo,"cyl"))
        transform = cartesian_cylindrical;
    else if (streq(ctypei,"cyl") && streq(ctypeo,"cart"))
        transform = cylindrical_cartesian;
    else
        error("Unimplemented ctype i/o : %s -> %s",ctypei,ctypeo);
    dprintf(0,"converting from %s to %s\n",ctypei,ctypeo);
    
    instr = stropen(getparam("in"), "r");   
    outstr = stropen(getparam("out"), "w");

    get_history(instr);
    put_history(outstr);		
    for (;;) {
    	get_history(instr);		/* skip over stuff we can forget */
        if (!get_tag_ok(instr, SnapShotTag))
		break;			/* done with work in loop */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0) {
	    continue;       /* just skip it's probably a diagnostics */
        }

        if ((bits & TimeBit) == 0)
	    tsnap = 0.0;
        else if (!streq(times,"all") && !within(tsnap, times, TIMEFUZZ))
            continue;
        dprintf (1,"Transforming snapshot at time= %f bits=0x%x\n",tsnap,bits);
#if 0
        Qmass  = MassBit & bits;       
        Qphase = PhaseSpaceBit & bits;
        Qacc   = AccelerationBit & bits;
        Qaux   = AuxBit & bits;
        Qkey   = KeyBit & bits;
#endif

        for (bp = btab; bp < btab+nbody; bp++) {
            (transform)(Pos(bp),Vel(bp),Acc(bp));
        }
        put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    }
}

/*
 *  Name            Mnemonic    Default CTYPE(1,2,3)    PROJP
 *  ----            --------    --------------------    -----
 *
 *  cartesian:      CART        X, Y, Z
 *  spherical:      SPH         RAD, PHI, THETA
 *  cylindrical:    CYL         RAD, PHI, Z
 *
 */
static int i1=0, o1=0,       /* default first axis */
           i2=1, o2=1,       /* default second axis */
           i3=2, o3=2;       /* default third axis */
           
void cartesian_spherical(real *p, real *v, real *a)
{
    real r2, rad, phi, theta;

	error("need to fix this routine");

    r2  = sqrt(sqr(a[i1]) + sqr(a[i2]));
    rad = sqrt(sqr(r2) + sqr(a[i3]));
    if (r2>0.0)
        phi = atan2(a[i2],a[i1]);
    else
        phi = 0.0;
    if (a[i3]==0.0)
        theta=HALF_PI;
    else
        theta=atan(r2/a[i3]);
    if (theta < 0) theta += PI;
    a[o1] = rad;
    a[o2] = phi;
    a[o3] = theta;
}

void spherical_cartesian(real *p, real *v, real *a)
{
    real r2, x, y, z;

	error("need to fix this routine");

    z  = a[i1] * cos(a[i3]);
    r2 = a[i1] * sin(a[i3]);
    x = r2 * cos(a[i2]);    
    y = r2 * sin(a[i2]);
    a[o1] = x;
    a[o2] = y;
    a[o3] = z;
}

void cartesian_cylindrical(real *p, real *v, real *a)
{
    real r2, phi, z, vx, vy, vz, cosp, sinp;

    r2  = sqrt(sqr(p[i1]) + sqr(p[i2]));
    if (r2>0.0) {
        phi = atan2(p[i2],p[i1]);
        cosp = p[i1]/r2;
        sinp = p[i2]/r2;
    } else {
        phi = 0.0;
        cosp = 1.0;     /* care with vel & acc */
        sinp = 0.0;
    }
    z = p[i3];
    p[o1] = r2;
    p[o2] = phi;
    p[o3] = z;

    vx = v[i1];
    vy = v[i2];
    vz = v[i3];
    v[o1] =  cosp*vx + sinp*vy;
    v[o2] = -sinp*vx + cosp*vy;
    v[o3] = vz;
    
}

void cylindrical_cartesian(real *p, real *v, real *a)
{
    real r2, x, y, z, vr, vp, vz, cosp, sinp;

    z = p[i3];
    cosp = cos(p[i2]);    
    sinp = sin(p[i2]);
    x = p[i1] * cosp;
    y = p[i1] * sinp;
    p[o1] = x;
    p[o2] = y;
    p[o3] = z;

    vr = v[i1];
    vp = v[i2];
    vz = v[i3];
    v[o1] = vr*cosp - vp*sinp;
    v[o2] = vr*sinp + vp*cosp;
    v[o3] = vz;

    /* and do the same for acc */    
}

