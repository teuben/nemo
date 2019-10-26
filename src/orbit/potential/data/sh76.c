/*
 * sh76.c:  2D bar potential used by Sanders and Huntley (1976)  1976ApJ...209...53S
 *
 * density perturbation:
 *      A r^-\alpha (1+\eps*cos(2\theta-2\Omega t))
 * eq.(2) give the radial and tangential forces
 * the rotcur of the unperturbed density is
 *      V(r) = \sqrt{2\pi G A c_1} r^{(1-\alpha)/2}
 *
 *	26-oct-2019 Created after Sanders (2019)             PJT
 */

/*CTEX
 *  {\bf potname=sh76
 *       potpars={\it $\Omega,B,a,A,r_0,i_0,j$}}
 *  
 *
 *  This bar potential was used by Sanders and Huntley (1976) and
 *  revived by Sanders (2019)
 *
 */
 

#include <stdinc.h>
#include <potential_float.h>

static double omega     = 8.1;         /* \omega: [Myr^-1]   Pattern speed */
static double sh_A      = 965.0;       /* mass */
static double sh_alpha  = 1.5;         /* exponent */
static double sh_eps    = 0.0;         /* bar strenght */
static double sh_L      = 1.0;         /* length scale kpc */

static double c1, beta;
static double tpg       = 2.70378-05;  /* 2 * pi * G  ; in SH76 galactic units  */
                                       /* G = 1/2.32385e5  in Msolar, kpc, km/s */

void inipotential (int  *npar, double *par, char *name)
{
    int n;

    n = *npar;
    if (n>0) omega     = par[0];
    if (n>1) sh_A      = par[1];
    if (n>2) sh_alpha  = par[2];
    if (n>3) sh_eps    = par[3];
    if (n>4) sh_L      = par[4];

    if (n>5) warning("sh76: npar=%d only 5 parameters accepted",n);

    dprintf (1,"INIPOTENTIAL Sanders-Huntley 1976 %s\n",name);
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  A= %g  alpha= %g eps=%g \n", sh_A, sh_alpha, sh_eps);

    par[0] = omega;
    beta = (2-sh_alpha)*sh_eps/(sh_alpha*(3-sh_alpha));
    c1 = (gamma(0.5*(2-sh_alpha)) * gamma(0.5*(sh_alpha+1))) /
         (gamma(0.5*sh_alpha) * gamma(0.5*(3-sh_alpha)));
    if (c1 < 0) {
      warning("sh76: c1=%g wrong sign, fixing it",c1);
      c1 = -c1;
    }
    dprintf (1,"  c1= %g beta= %g\n", c1, beta);
}

void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double r, r2, ra2, ra2sq;
    double theta, psi, ksi, phi, tmp1, fr, ft, sinp, cosp;
    double dpsi, dksi, dphi, vr2; 
        
    r2 = sqr(pos[0]) + sqr(pos[1]);
    r = sqrt(r2);
    if (r2>0.0)
        theta = atan2(pos[1], pos[0]);
    else
        theta = 0.0;

    /* rotcur */
    //vr2 = tp2 * sh_A * c1 * pow(1-sh_alpha);

    fr = -tpg * sh_A * c1 * pow(r, -sh_alpha-1);
    ft = 0.0;
    
#if 0
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

#endif
    
    acc[0] = -fr*pos[0] + ft*pos[1];
    acc[1] = -fr*pos[1] - ft*pos[0];
    acc[2] = 0.0;
    
}
