/*
 * rh84.c:  2D spiral & bar potential that was used by:
 *	    For counterclockwise streaming, this spiral is a trailing
 *	    spiral when the pitch angle (i0) is positive.
 *	    Within a radius r0 the potential becomes barlike, with
 *	    the bar along the X axis.
 *	    At large radii the spiral is logarithmic.
 *
 *		Roberts & Haussman (1984: ApJ 277, 744) 
 *	        Roberts, Huntley & v.Albada (1979: ApJ 233, 67)
 *
 *	11-jun-92 Created		PJT
 *	14-jun-92 Fixed the j-power term	PJT
 *      19-jul-92 added the forces              PJT
 *	25-jul-92 A=0 still done the hard way   PJT
 *	   oct-93 get_pattern
 *	19-sep-04 float/double			PJT
 */

/*CTEX
 *  {\bf potname=rh84
 *       potpars={\it $\Omega,B,a,A,r_0,i_0,j$}}
 *  
 *
 *  This 2D spiral and bar potential was used by Robert and collaborators
 *  in the 70s and 80s.
 *  For counterclockwise streaming, this spiral is a trailing
 * spiral when the pitch angle ($i_0$) is positive.
 * Within a radius $r_0$ the potential becomes barlike, with
 * the bar along the X axis.
 * At large radii the spiral is logarithmic.
 * References:
 *
 *		Roberts \& Haussman (1984: ApJ 277, 744) 
 *
 *	        Roberts, Huntley \& v.Albada (1979: ApJ 233, 67)
 */
 

#include <stdinc.h>
#include <potential_float.h>

static double omega = 0.0;         /* \omega: [Myr^-1]   Pattern speed */
static double rh_B  = 0.0576;      /* B: [Myr^-1]        Mass */
static double rh_a  = 7.0;         /* a: [kpc]           Lenght scale */
static double rh_A  = 0.067;       /* A: dimensionless   Spiral Perturbation */
static double rh_r0 = 1.0;         /* r0: [kpc]          bar -> spiral */
static double rh_i0 = 10.0;        /* i0: [degrees]      pitch angle */
static double rh_j  = 5;           /* j: dimensionless   coefficient */

static double gm, aa, jt;          /* handy constants to keep around */


void inipotential (int  *npar, double *par, char *name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) rh_B  = par[1];
    if (n>2) rh_a  = par[2];
    if (n>3) rh_A  = par[3];
    if (n>4) rh_r0 = par[4];
    if (n>5) rh_i0 = par[5];
    if (n>6) rh_j  = par[6];

    if (n>7) warning("rh84: npar=%d only 7 parameters accepted",n);

    dprintf (1,"INIPOTENTIAL Roberts-Hausman 1984: %s\n",name);
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  B= %g  a= %g \n", rh_B, rh_a);
    dprintf (1,"  A= %g  r0=%g i0=%g j=%g\n", rh_A, rh_r0, rh_i0, rh_j);

    gm = sqr(rh_B)*qbe(rh_a);
    aa = 0.2*rh_A*sqr(rh_a);
    jt = 2/(rh_j*tan(rh_i0*PI/180.0));
    par[0] = omega;
}

void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double r, r2, ra2, ra2sq;
    double theta, psi, ksi, phi, tmp1, fr, ft, sinp, cosp;
    double dpsi, dksi, dphi;
        
    r2 = sqr(pos[0]) + sqr(pos[1]);
    r = sqrt(r2);
    ra2 = r2 + sqr(rh_a);
    ra2sq = sqrt(ra2);
    if (r2>0.0)
        theta = atan2(pos[1], pos[0]);
    else
        theta = 0.0;

    psi = -gm/ra2sq;


	/* this indented section is really when aa > 0, but also does aa=0 */

    	tmp1 = 1 + pow(r/rh_r0,rh_j);
        ksi =  aa*r2/sqr(ra2);
        phi =  jt*log(tmp1);
        sinp = sin(2*theta+phi);
        cosp = cos(2*theta+phi);
        *pot = psi*(1+ksi*cosp);



    dpsi = gm * r / (ra2*ra2sq);                            /* d(psi)/dr */
    dksi = aa * 2 * r * (sqr(rh_a) - r2) / qbe(ra2);        /* d(ksi)/dr */
    dphi = jt * rh_j * (tmp1-1)/(tmp1*r);                   /* d(phi)/dr */
    

    fr = (dpsi*(1+ksi*cosp) + psi*(dksi*cosp-ksi*sinp*dphi))/r;
    ft = -2*psi*ksi*sinp/r2;

    acc[0] = -fr*pos[0] + ft*pos[1];
    acc[1] = -fr*pos[1] - ft*pos[0];
    acc[2] = 0.0;
    
}
