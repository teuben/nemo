/* create - MFRAC, create_plum, com_rest */

/*
 *  create.c: creates a Plummer model for newton0
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
#include  "newton0.h"

static void com_rest(bodyptr bodies, int npart);

extern real  xrandom(real, real);            /* random number generator                */

#define  MFRAC  0.999       /* corresponding to a radial cut off at r = 22.8 */
                            /* in the current VIRIAL units. Equivalently,    */
                            /* r = 29.7 in units of the half-mass radius, or */
                            /* r = 60.2 in units of the core radius (defined */
                            /* as the distance where the projected surface   */
                            /* has dropped to half its value). In structural */
                            /* units, in which the density is proportional   */
                            /* to (1 + r*r)^(-5), r = 38.7 .                 */

/*-----------------------------------------------------------------------------
 *  create_plum  --  builds a local many-body system according to a Plummer
 *                   model, in VIRIAL units (M=G=-4E=1, with E the total
 *                   energy), truncated to contain only particles within the
 *                   inner fraction MFRAC of the cumulative mass distribution
 *                   of the ideal Plummer model.
 *		     note: the distribution function is spherical and
 *		           isotropic, and is a polytrope of index n = 5.
 *		     litt: S.J. Aarseth, M. Henon and R. Wielen (1974),
 *		           Astron. and Astrophys. 37, p. 183.
 *                   accepts: npart: the number of particles.
 *                   returns: plum: a pointer to the newly created (using
 *                                  malloc) array of particles representing
 *			            the Plummer model realization in which
 *                                  all particles have equal masses.
 *                   the random number generator used here is xrandom();
 *                   it has to be started with a call  srandom(seed)
 *                   elsewhere (in newton0 this is done in  create_system()
 *                   in the file  newton0.c ).
 *                   NOTE: after sprinkling in particles according to a Plummer
 *                         distribution, the system is shifted in position and
 *                         velocity so as to put the center of mass at rest at
 *                         the coordinate center. This may result in the system
 *                         acquiring a somewhat larger radial extent than
 *                         anticipated from the cut-off in cumulative mass
 *                         MFRAC.
 *-----------------------------------------------------------------------------
 */
bodyptr  create_plum(npart)
int  npart;
    {
    real  partmass;		/* equal mass of each particle		  */
    real  radius;		/* absolute value of position vector      */
    real  velocity;		/* absolute value of velocity vector      */
    real  theta, phi;		/* direction angles of above vectors      */
    real  x, y;		        /* for use in rejection technique         */
    real  scalefactor;          /* for converting between different units */
    real  inv_scalefactor;      /* inverse scale factor                   */
    real  sqrt_scalefactor;     /* sqare root of scale factor             */
    bodyptr  part_i;            /* pointer to an individual particle      */
    bodyptr  plum;              /* this pointer will be returned          */

    if (NDIM < 3)
        error("create_plum: NDIM = %d but should be at least 3\n", NDIM);
                                           /* could be 4 with regularization */
    plum = mk_bodies(npart);

    partmass = 1.0 / npart;
/*
 * now we construct the individual particles:
 */
    for (part_i = plum; part_i - plum < npart; part_i++)
	{
	Mass(part_i) = partmass;
/*
 * the position coordinates are determined by inverting the cumulative
 * mass-radius relation, with the cumulative mass drawn randomly from
 * [0, MFRAC]; cf. Aarseth et al. (1974), eq. (A2).
 */
	radius = 1.0 / sqrt( pow (xrandom(0.0,MFRAC), -2.0/3.0) - 1.0);
	theta = acos(xrandom(-1.0, 1.0));
	phi = xrandom(0.0, TWO_PI);
	Pos(part_i)[0] = radius * sin( theta ) * cos( phi );
	Pos(part_i)[1] = radius * sin( theta ) * sin( phi );
        Pos(part_i)[2] = radius * cos( theta );
/*
 * the velocity coordinates are determined using von Neumann's rejection
 * technique, cf. Aarseth et al. (1974), eq. (A4,5).
 * First we take initial values for x, the ratio of velocity and escape
 * velocity (q in Aarseth et al.), and y, as a trick to enter the body of the
 * while loop.
 */
	x = 0.0;
	y = 0.1;
/*
 * Then we keep spinning the random number generator until we find a pair
 * of values (x,y), so that y < g(x) = x*x*pow( 1.0 - x*x, 3.5) . Whenever
 * an y-value lies above the g(x) curve, the (x,y) pair is discarded, and
 * a new pair is selected. The value 0.1 is chosen as a good upper limit for
 * g(x) in [0,1] : 0.1 > max g(x) = 0.092 for 0 < x < 1.
 */
	while (y > x*x*pow( 1.0 - x*x, 3.5))
	    {
	    x = xrandom(0.0,1.0);
	    y = xrandom(0.0,0.1);
	    }
/*
 * If y < g(x), proceed to calculate the velocity components:
 */
	velocity = x * sqrt(2.0) * pow( 1.0 + radius*radius, -0.25);
	theta = acos(xrandom(-1.0, 1.0));
	phi = xrandom(0.0,TWO_PI);
	Vel(part_i)[0] = velocity * sin( theta ) * cos( phi );
	Vel(part_i)[1] = velocity * sin( theta ) * sin( phi );
	Vel(part_i)[2] = velocity * cos( theta );
	}

/*
 * Calculating the coordinates was easiest in STRUCTURAL units, used above;
 * conversion to VIRIAL units will be performed below.
 * VIRIAL units are defined by the requirement that the total mass M, the
 * gravitational constant G, and the total energy E of the system satisfy
 * the relations M=G=-4E=1.
 *
 *   Recipe for scaling to the proper system of units:
 *
 * Since G = M = 1, if we switch from a coordinate system with
 * length unit  r_old  to a coordinate system with length unit  r_new ,
 * the length units simply scale by a factor  C = r_new / r_old .
 * Consequently, the coordinate values of physical quantities
 * such as positions should transform inversely to maintain the same
 * coordinate-invariant meaning. Similarly, the square of the velocities
 * should transform inversely proportional to the positions,
 * since  GM = 1  (cf. a relation such as  v*v = G*M/r ).
 * To summarize: If
 *                      r_unit(new) = C * r_unit(old)  ,
 *               then
 *                      pos(new) = (1/C) * pos(old)
 *               and
 *                      vel(new) = sqrt(C) * vel(old)  .
 */
    scalefactor = 16.0 / (3.0 * PI);
    inv_scalefactor = 1.0 / scalefactor;
    sqrt_scalefactor = sqrt( scalefactor );
/*
 * Now transform to the VIRIAL coordinates by applying
 * the scaling factors to the positions and velocities:
 */
    for (part_i = plum; part_i - plum < npart; part_i++)
        {
	Pos(part_i)[0] *= inv_scalefactor;
	Pos(part_i)[1] *= inv_scalefactor;
	Pos(part_i)[2] *= inv_scalefactor;
	Vel(part_i)[0] *= sqrt_scalefactor;
	Vel(part_i)[1] *= sqrt_scalefactor;
	Vel(part_i)[2] *= sqrt_scalefactor;
        }

#ifdef REGULARIZATION
    for (part_i = plum; part_i - plum < npart; part_i++)
        {
	Pos(part_i)[3] = 0.0;
	Vel(part_i)[3] = 0.0;
        }
#endif

    com_rest(plum, npart);      /* to shift to the center-of-mass rest frame */

    return(plum);            /* return the pointer to the new particle array */
    }

/*-----------------------------------------------------------------------------
 *  com_rest  --  transforms to a coordinate system in which the center of mass
 *                remains at rest in the origin.
 *-----------------------------------------------------------------------------
 */
local void  com_rest(bodies, npart)
bodyptr  bodies;
int  npart;
    {
    real  total_mass;
    real  pos_com[NDIM];
    real  vel_com[NDIM];
    real  mass_i;
    real  weighted_pos[NDIM];
    real  weighted_vel[NDIM];
    bodyptr  body_i;
    
    total_mass = 0.0;
    CLRV(pos_com);
    CLRV(vel_com);

    for (body_i = bodies; body_i - bodies < npart; body_i++)
        {
        mass_i = Mass(body_i);
	total_mass += mass_i;
        MULVS(weighted_pos, Pos(body_i), mass_i);
        MULVS(weighted_vel, Vel(body_i), mass_i);
	INCADDV(pos_com, weighted_pos);
	INCADDV(vel_com, weighted_vel);
	}

    DIVVS(pos_com, pos_com, total_mass);
    DIVVS(vel_com, vel_com, total_mass);

    for (body_i = bodies; body_i - bodies < npart; body_i++)
        {
	INCSUBV(Pos(body_i), pos_com);
	INCSUBV(Vel(body_i), vel_com);
        }
    }

/* endof: create.c */
