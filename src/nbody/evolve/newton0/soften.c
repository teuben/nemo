/* soften.c - softener, soft_plummer, soft_power4, soft_power8,
              soft_exponential, soft_sphere, soft_shell, soft_spline */
/*
 *  soften.c:  procedures for softening Newtonian gravity, for newton0
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
#include  "newton0.h"

static void soft_plummer(real separ[3 ], real r0_soft, real *r_inv_eff, real *r_inv_cubed_eff);
static void soft_power4(real separ[3 ], real r0_soft, real *r_inv_eff, real *r_inv_cubed_eff);
static void soft_power8(real separ[3 ], real r0_soft, real *r_inv_eff, real *r_inv_cubed_eff);
static void soft_exponential(real separ[3 ], real r0_soft, real *r_inv_eff, real *r_inv_cubed_eff);
static void soft_sphere(real separ[3 ], real r0_soft, real *r_inv_eff, real *r_inv_cubed_eff);
static void soft_shell(real separ[3 ], real r0_soft, real *r_inv_eff, real *r_inv_cubed_eff);
static void soft_spline(real separ[3 ], real r0_soft, real *r_inv_eff, real *r_inv_cubed_eff);

/*-----------------------------------------------------------------------------
 *  softener  --  returns the pointer  ptr  to the function which carries out
 *                the desired method of softening. Each of these functions
 *                computes the softened inverse distance used for the softened
 *                potential and the corresponding softened inverse cubed
 *                distance which is used in the force calculation.
 *-----------------------------------------------------------------------------
 */
proc  softener(type_of_softening)
string  type_of_softening;
    {
    if (streq("plummer", type_of_softening))
	return( soft_plummer );
    else if (streq("spline", type_of_softening))
	return( soft_spline );
    else if (streq("power4", type_of_softening))
	return( soft_power4 );
    else if (streq("power8", type_of_softening))
	return( soft_power8 );
    else if (streq("exponential", type_of_softening))
	return( soft_exponential );
    else if (streq("sphere", type_of_softening))
	return( soft_sphere );
    else if (streq("shell", type_of_softening))
	return( soft_shell );
    else
	error("softener: %s not implemented as a type of softening\n",
                                                            type_of_softening);
    }

/*-----------------------------------------------------------------------------
 *  soft_plummer  --  computes the effective inverse distance and the effective
 *                    cubed inverse distance which replaced the real inverse
 *                    distance and its cube in the calculation of softened
 *                    potential and force.
 *                    In this particular softening method, point particles
 *                    are effectively replaced by small plummer model
 *                    density distributions: 1/r -> (r^2 + r0^2)^{-1/2}.
 *                    accepts: separ: the separation vector, i.e. the
 *                                    position vector of the other particle
 *                                    starting at our particle;
 *                             r0_soft: softening length;
 *                             r_inv_eff: a pointer to the effective inverse
 *                                        distance;
 *                             r_inv_cubed_eff: a pointer to the effective
 *                                        inverse cubed distance.
 *                    effect: sets the proper values *r_inv_eff and
 *                            *r_inv_cubed_eff .
 *-----------------------------------------------------------------------------
 */
local void  soft_plummer(separ, r0_soft, r_inv_eff, r_inv_cubed_eff)
real  separ[NDIM];
real  r0_soft;
real *r_inv_eff;
real *r_inv_cubed_eff;
    {
    real  dist_sqr;           /* squared distance between the particles      */

    DOTVP(dist_sqr, separ, separ);
    *r_inv_eff = 1.0 / sqrt(dist_sqr + r0_soft*r0_soft);
    *r_inv_cubed_eff = *r_inv_eff * *r_inv_eff * *r_inv_eff;
    }

/*-----------------------------------------------------------------------------
 *  soft_power4  --  computes softened 1/R and 1/R^3 according to a type
 *                   of softening with shorter range than is standard:
 *                             1/r -> (r^4 + r0^4)^{-1/4} 
 *                   arguments: see  soft_plummer()  above.
 *-----------------------------------------------------------------------------
 */
local void  soft_power4(separ, r0_soft, r_inv_eff, r_inv_cubed_eff)
real  separ[NDIM];
real  r0_soft;
real *r_inv_eff;
real *r_inv_cubed_eff;
    {
    real  dist_sqr;           /* squared distance between the particles      */

    DOTVP(dist_sqr, separ, separ);
    *r_inv_eff = 1.0 / sqrt(sqrt(dist_sqr*dist_sqr
                                         + r0_soft*r0_soft*r0_soft*r0_soft));
    *r_inv_cubed_eff = dist_sqr * *r_inv_eff * *r_inv_eff * *r_inv_eff
                                             * *r_inv_eff * *r_inv_eff;
    }

/*-----------------------------------------------------------------------------
 *  soft_power8  --  computes softened 1/R and 1/R^3 according to a type
 *                   of softening with shorter range than is standard:
 *                             1/r -> (r^8 + r0^8)^{-1/8} 
 *                   arguments: see  soft_plummer()  above.
 *-----------------------------------------------------------------------------
 */
local void  soft_power8(separ, r0_soft, r_inv_eff, r_inv_cubed_eff)
real  separ[NDIM];
real  r0_soft;
real *r_inv_eff;
real *r_inv_cubed_eff;
    {
    real  dist_sqr;               /* squared distance between the particles  */
    real  r0_sqr;                 /* squared softening length                */
    real  eff_inv_dist_sqr;       /* squared effective inverse distance      */

    r0_sqr = r0_soft * r0_soft;
    DOTVP(dist_sqr, separ, separ);
    *r_inv_eff = 1.0 / sqrt(sqrt(sqrt(dist_sqr*dist_sqr*dist_sqr*dist_sqr
                                         + r0_sqr*r0_sqr*r0_sqr*r0_sqr)));
    eff_inv_dist_sqr = *r_inv_eff * *r_inv_eff;
    *r_inv_cubed_eff = dist_sqr * dist_sqr * dist_sqr * *r_inv_eff
                                      * eff_inv_dist_sqr * eff_inv_dist_sqr
                                      * eff_inv_dist_sqr * eff_inv_dist_sqr;
    }

/*-----------------------------------------------------------------------------
 *  soft_exponential  --  computes softened 1/R and 1/R^3 according to a type
 *                        of softening with shorter range than standard:
 *                                 1/r -> (r + exp(-r/r0))^{-1} 
 *                        arguments: see  soft_plummer()  above.
 *-----------------------------------------------------------------------------
 */
local void  soft_exponential(separ, r0_soft, r_inv_eff, r_inv_cubed_eff)
real  separ[NDIM];
real  r0_soft;
real *r_inv_eff;
real *r_inv_cubed_eff;
    {
    real  dist;               /* distance between the particles              */
    real  exp_min_ratio;      /* exp( - r/r0 )                               */

    ABSV(dist, separ);
    exp_min_ratio = exp( -dist / r0_soft );
    *r_inv_eff = 1.0 / (dist + r0_soft * exp_min_ratio);
    *r_inv_cubed_eff = *r_inv_eff * *r_inv_eff * (1.0 - exp_min_ratio) / dist;
    }

/*-----------------------------------------------------------------------------
 *  soft_sphere  --  computes softened 1/R and 1/R^3 according to a type
 *                   of softening in which effectively point particles
 *                   are replaced by small (penetrable) homogeneous spherical
 *                   density distributions.
 *                   arguments: see  soft_plummer()  above.
 *-----------------------------------------------------------------------------
 */
local void  soft_sphere(separ, r0_soft, r_inv_eff, r_inv_cubed_eff)
real  separ[NDIM];
real  r0_soft;
real *r_inv_eff;
real *r_inv_cubed_eff;
    {
    real  dist;               /* distance between the particles              */

    ABSV(dist, separ);
    if (dist > r0_soft)
	{
        *r_inv_eff = 1.0 / dist;
        *r_inv_cubed_eff = *r_inv_eff * *r_inv_eff * *r_inv_eff;
	}
    else
	{
	*r_inv_eff = 0.5 * (3.0 - (dist*dist)/(r0_soft*r0_soft) ) / r0_soft;
	*r_inv_cubed_eff = 1.0 / (r0_soft * r0_soft * r0_soft);
	}
    }

/*-----------------------------------------------------------------------------
 *  soft_shell  --  computes softened 1/R and 1/R^3 according to a type
 *                  of softening in which effectively point particles
 *                  are replaced by density distributions in the form of
 *                  small (penetrable) spherical shells.
 *                  arguments: see  soft_plummer()  above.
 *-----------------------------------------------------------------------------
 */
local void  soft_shell(separ, r0_soft, r_inv_eff, r_inv_cubed_eff)
real  separ[NDIM];
real  r0_soft;
real *r_inv_eff;
real *r_inv_cubed_eff;
    {
    real  dist;               /* distance between the particles              */

    ABSV(dist, separ);
    if (dist > r0_soft)
	{
        *r_inv_eff = 1.0 / dist;
        *r_inv_cubed_eff = *r_inv_eff * *r_inv_eff * *r_inv_eff;
	}
    else
	{
	*r_inv_eff = 1.0 / r0_soft;
	*r_inv_cubed_eff = 0.0;
	}
    }

/*-----------------------------------------------------------------------------
 *  soft_spline  --  computes softened 1/R and 1/R^3 according to a type
 *                   of softening with shorter range than is standard,
 *                   according to a cubic spline.
 *                   litt.: Monaghan, J.J. 1985, Comp. Phys. Rep. 3, 71.
 *                   arguments: see  soft_plummer()  above.
 *-----------------------------------------------------------------------------
 */
local void  soft_spline(separ, r0_soft, r_inv_eff, r_inv_cubed_eff)
real  separ[NDIM];
real  r0_soft;
real *r_inv_eff;
real *r_inv_cubed_eff;
    {
    real  x;                      /* r/h , i.e. dist / r0_sqr                */
    real  x2;                     /* x * x                                   */
    real  x3;                     /* x * x * x                               */
    real  dist;                   /* distance between the particles          */
    real  dist_sqr;               /* squared distance between the particles  */
    real  eff_inv_dist_sqr;       /* squared effective inverse distance      */

    DOTVP(dist_sqr, separ, separ);
    dist = sqrt(dist_sqr);

    *r_inv_eff = 1.0 / dist;
    *r_inv_cubed_eff = *r_inv_eff * *r_inv_eff * *r_inv_eff;

    if (dist >= 2.0 * r0_soft)
        return;

    x = dist / r0_soft;
    x2 = x * x;
    x3 = x * x2;

    if (dist >= r0_soft)
	{
	*r_inv_eff *= (x3*(x3 - 9.0*x2 + 30.0*x - 40.0) + 48*x -2.0)/30.0;
        *r_inv_cubed_eff *= (x3*(-5.0*x3 + 36.0*x2 - 90.0*x + 80.0) -2.0)/30.0;
	}
    else
	{
	*r_inv_eff *= (x3*(-3.0*x3 + 9.0*x2 - 20.0) + 42*x)/30.0;
        *r_inv_cubed_eff *= x3*(15.0*x3 - 36.0*x2 + 40.0)/30.0;
	}
    }

/* endof: soften.c */
