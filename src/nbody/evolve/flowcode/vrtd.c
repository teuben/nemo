/* 
 * VRT: flow potential (and density weighed interpolated forces) from an image :
 *      2D version
 *      z-forces are merely returned as 0.0
 *	4 or 5 images are expected: 
 *          2D vr
 *          2D vt
 *	    2D den (optional - if absent, like old vrt.c)
 *          1D R-axis descriptor
 *          1D T-axis descriptor
 *
 *	13-apr-96  cloned off vxy.c
 *	12-sep-96  added density, but keep it backwards compatible
 *
 *	TODO: - move the vscale up into the initialization routine,
 *	        and implement the rscale option
 *            - make searching a binary one
 */

#include <stdinc.h>
#include <filestruct.h>
#include <image.h>

#define CCD_VERSION "flowcode:vrtd V1.1 12-sep-96"

local double   omega = 0.0;		/* pattern speed */
local double   vscale = 100.0;		/* rescale velocity unit */
local double   rscale = 1.0;		/* rescale radius unit */

local stream   potstr = NULL;
local imageptr vr = NULL, vt = NULL, den = NULL, r = NULL, t = NULL;
local int      nr, np;
local real     *rads, *phis;
local real     dphi;

local int binsearch(real, real *, int);

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0)  omega = par[0];
    if (n>1)  vscale = par[1];
    if (n>2)  warning("inipotential(flowrt): npar=%d only 2 parameter accepted",n);

    dprintf(1,"INIPOTENTIAL: %s: %s\n",CCD_VERSION,name);
    dprintf(1,"  Parameters:  Omega=%g Vscale=%g\n",omega,vscale);

    potstr = stropen (name, "r");            /* open the image */
    n = read_image (potstr,&vr);              /* read the 1st image */
    if (n) n = read_image (potstr,&vt);              /* read the 2nd image */
    if (n) n = read_image (potstr,&den);            /* read the 3rd image */
    if (n) n = read_image (potstr,&r);              /* read the 4th image */
    if (n==0) error("%s: Could not read 4 images from %s",CCD_VERSION,name);
    if (n) n = read_image (potstr,&t);              /* read the 5th image */
    if (n==0) {
	warning("%s: no 5th image, assuming you have a VRT and not VRTD",
			CCD_VERSION);
	r = den;
	t = r;
    }



    nr = Nx(vr);    /* set the radial and tangential dimension */
    np = Ny(vr);    /* and check if other image confirm */
    if (Nx(vt) != nr)  error("incorrect R dimension of VT image");
    if (Ny(vt) != np)  error("incorrect T dimension of VT image");
    if (Nx(den)!= nr)  error("incorrect R dimension of DEN image");
    if (Ny(den)!= np)  error("incorrect T dimension of DEN image");
    if (Nx(r)  != nr)  error("incorrect R dimension of R image");
    if (Ny(r)  != 1)   error("incorrect 1 dimension of R image");
    if (Nx(t)  != np)  error("incorrect T dimension of T image");
    if (Ny(t)  != 1)   error("incorrect 1 dimension of T image");
    if (vscale == 0.0) error("Vscale = 0.0 not allowed");

    rads = Frame(r);
    phis = Frame(t);
    dprintf(0,"Rad: %g - %g  Theta: %g - %g\n",
        rads[0], rads[nr-1], phis[0], phis[np-1]);

    par[0] = omega;
}

#define VR(ix,iy)	MapValue(vr,ix,iy)
#define VT(ix,iy)	MapValue(vt,ix,iy)
#define DEN(ix,iy)	MapValue(den,ix,iy)

void potential(int *ndim,double *pos,double *acc,double *pot,double *time)
{
    real rad, phi, phi_orig, phi1, phi2, rad1, rad2;
    real x, y, c1, c2, c3, c4, a1, a2, vrad, vtan;
    real d1, d2, d3, d4, e1, e2;
    int ix, iy, ip, jp, ir, jr, i, j;
    bool mirror;

    x = pos[0];
    y = pos[1];
    rad = sqrt(x*x + y*y);
    phi = atan2(y,x);               /* make sure phi in -pi/2 : pi/2  */

    if (rad > rads[nr-1] || rad < rads[0]) {

        vrad = vtan = 0.0;

    } else {


        mirror = x < 0;
	phi_orig = phi;
        if (mirror) {
            if (phi >  HALF_PI) phi -= PI;
	    if (phi < -HALF_PI) phi += PI;
        } 

        for (ir=0; ir<nr-1; ir++)           /* simple search in Rad */
            if (rads[ir+1] > rad) break;
        jr = ir+1;
        rad1 = rads[ir];
        rad2 = rads[jr];
            
        if(phi<phis[0]) {                   /* special search in Phi */
                if (mirror) {                   /* 2nd Quadrant */
                    ip = 0;
                    jp = np-1;
                    phi1 = phis[0] + PI;
                    phi2 = phis[np-1];
                    phi = phi_orig;
                } else {                        /* 4th Quadrant */
                    ip = np-1;
                    jp = 0;
                    phi1 = phis[np-1] - PI;
                    phi2 = phis[0];
                }
        } else if (phi>phis[np-1]) {        /* because of half-symmetry */
                if (mirror) {                   /* 3rd Quadrant */
                    ip = np-1;
                    jp = 0;
                    phi1 = phis[np-1] - PI;
                    phi2 = phis[0];
                    phi = phi_orig;
                } else {                        /* 1st Quadrant */
                    ip = 0;
                    jp = np-1;
                    phi1 = phis[0] + PI;
                    phi2 = phis[np-1];
                }
        } else {                            /* simple in the middle */
                for (ip=0; ip<np-1; ip++)
                    if (phis[ip+1] > phi) break;
                jp = ip+1;
                phi1 = phis[ip];
                phi2 = phis[jp];
                if (ip>=np-1 || ir>=nr-1) {     /* should never happen */
                    error("### POTENTIAL: Odd: %d %d: ip=%d ir=%d",i,j,ip,ir);
                }
        }
        dprintf(2,"(x,y)%d %d (r,p) %g %g [@ %d %d] LL: %g %g\n",
                    i,j,rad,phi,ir,ip,rads[ir],phis[ip]);
                
        c1 = (rad-rad2)/(rad1-rad2);
        c2 = (rad-rad1)/(rad2-rad1);
        c3 = (phi-phi2)/(phi1-phi2);
        c4 = (phi-phi1)/(phi2-phi1);

        d1 = DEN(ir,ip);
        d2 = DEN(jr,ip);
        d3 = DEN(ir,jp);
        d4 = DEN(jr,jp);

        e1 = d1*c1 + d2*c2;
        e2 = d3*c1 + d4*c2;

        a1 = (d1*VR(ir,ip)*c1 + d2*VR(jr,ip)*c2)/e1;
        a2 = (d3*VR(ir,jp)*c1 + d4*VR(jr,jp)*c2)/e2;
        vrad = (e1*a1*c3 + e2*a2*c4)/(e1*c3+e2*c4);

        a1 = (d1*VT(ir,ip)*c1 + d2*VT(jr,ip)*c2)/e1;
        a2 = (d3*VT(ir,jp)*c1 + d4*VT(jr,jp)*c2)/e2;
        vtan = (e1*a1*c3 + e2*a2*c4)/(e1*c3+e2*c4);

    }
    
    *pot = 0.0;
    if (rad > 0) {
	vtan -= omega * rad;
        acc[0] = (vrad * x - vtan * y) / (rad * vscale);
        acc[1] = (vrad * y + vtan * x) / (rad * vscale);
    } else {
        acc[0] = acc[1] = 0.0;
    }
    dprintf(1,"r,t=%g %g  vr,vt=%g %g vx,vy=%g %g\n",
            rad,phi,vrad,vtan,acc[0],acc[1]);
    acc[2] = 0.0;
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

