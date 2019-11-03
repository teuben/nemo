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
 *       potpars={\it $\Omega,A,\alpha,\epsilon$}}
 *
 *  This bar potential was used by Sanders and Huntley (1976) and
 *  also used in Sanders (2019).   The density perturbation is given
 *  by
 *
 * $$ 
 *    \sigma(r,\theta) = A r^{-\alpha} (1+\epsilon*\cos{2\theta})
 * $$
 *   and the potential
 * $$
 *    \Phi(r,\theta) = -2\pi G c_1 A r^{-\alpha+1} {1 \over {1-\alpha}} ( 1 + \beta (\alpha-1) \cos{2\theta})
 * $$
 * where
 * $$
 *  	\beta =  { {(2-\alpha)} \over { \alpha(3-\alpha)} }  \epsilon
 * $$
 * and
 * $$
 *  	c_1 = { { \Gamma{[{1\over 2}(2-\alpha)]}  \Gamma{[{1\over 2}(\alpha+1)]} }   \over
 *              { \Gamma{[{1\over 2}\alpha]}  \Gamma{[{1\over 2}(3-\alpha)]} } }
 * $$
 *  
 */
 

#include <stdinc.h>
#include <potential_float.h>

static double omega     = 8.1;         /* \omega: [Myr^-1]   Pattern speed */
static double sh_A      = 965.0;       /* mass */
static double sh_alpha  = 1.5;         /* exponent */
static double sh_eps    = 0.295;       /* oval distortion strength */
static double sh_L      = 1.0;         /* length scale kpc -- not used here */

static double c1, beta;
static double G2pi      = 2.70378-05;  /* 2 * pi * G  ; in SH76 galactic units  */
                                       /* G = 1/2.32385e5  in Msolar, kpc, km/s */

/*
rotcurves sh76 8.1,965,1.5,0 debug=1 mode=omega radii=0:28:0.1 yrange=0,16

gives that IRL is more like 2.5 kpc, not 7kpc.   Is their 2.pi.G correct? Mine is 2.70378-05

### nemo Debug Info: INIPOTENTIAL Sanders-Huntley 1976 
### nemo Debug Info:   Parameters : Pattern Speed = 8.100000
### nemo Debug Info:   A= 965  alpha= 1.5 eps=0 
### nemo Debug Info:   c1= 2.18844 beta= 0
### nemo Debug Info: oILR: 2.62013

this looks better, but why do i need to scale A by ~11.8

rotcurves sh76 '8.1,965*11.8,1.5,0' debug=1 mode=omega radii=0:28:0.01 yrange=0,16
### nemo Debug Info:   A= 11387  alpha= 1.5 eps=0 
### nemo Debug Info: OLR: 21.9891
### nemo Debug Info: CR: 15.005
### nemo Debug Info: oILR: 7.03018

*/

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
    /* note that gamma() is deprecated */
    c1 = (tgamma(0.5*(2-sh_alpha)) * tgamma(0.5*(sh_alpha+1))) /
         (tgamma(0.5*sh_alpha) * tgamma(0.5*(3-sh_alpha)));
    if (c1 < 0) {
      warning("sh76: c1=%g wrong sign for alpha=%g, fixing it",c1,sh_alpha);
      c1 = -c1;
    }
    dprintf (1,"  c1= %g beta= %g\n", c1, beta);
    if (sh_L < 0)
      warning("SH76:  returning density instead of potential");
    else
      warning("SH76:  only forces returned, potential = 0.0");
}

void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double r, r2, ra2, ra2sq;
    double theta, psi, ksi, phi, tmp1, fr, ft, sinp, cosp;
    double amr, vr2; 
        
    r2 = sqr(pos[0]) + sqr(pos[1]);
    if (r2>0.0) {
      r = sqrt(r2);
      theta = atan2(pos[1], pos[0]);
      amr = sh_A * pow(r, -sh_alpha);
    } else {
      r = 0.0;
      theta = 0.0;
      amr = 0.0;
    }

    /* rotcur */
    //vr2 = tp2 * sh_A * c1 * pow(1-sh_alpha);


    fr = -G2pi * c1 * amr;
    ft = 2 * fr * beta * sin(2*theta);
    fr *= (1 + beta * (sh_alpha-1) * cos(2*theta))/r;
    
    acc[0] = -fr*pos[0] + ft*pos[1];
    acc[1] = -fr*pos[1] - ft*pos[0];
    acc[2] = 0.0;

    *pot = -G2pi * c1 * amr * r / (1-sh_alpha) * (1 + beta * (sh_alpha-1) * cos(2*theta));

    /* unofficial feature for this potential:  pass density instead of potential */
    if (sh_L < 0) 
      *pot = amr * (1 + sh_eps * cos(2*theta));
}
