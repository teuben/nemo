/*  scat.h - CARTESIAN, SPHERICAL, TRUEANOMALY, ECCANOMALY, MEANANOMALY,
             PERIPASSAGE, SCATTERING, TSCATTERING, iskeptype */

/*
 *  scat.h: header file for gravitational scattering experiments
 *		(c) 1988  Piet Hut  Princeton, NJ, USA
 */

/*-----------------------------------------------------------------------------
 *  macros for the different types of phase space coordinates, in kepler orbits
 *-----------------------------------------------------------------------------
 */
#define   CARTESIAN     11
#define   SPHERICAL     12  /* note: Carthesian velocities along unit vectors
                                     in the radial, theta and phi directions */
#define   TRUEANOMALY   21
#define   ECCANOMALY    22
#define   MEANANOMALY   23
#define   PERIPASSAGE   24

#define   SCATTERING    31
#define   TSCATTERING   32

/*-----------------------------------------------------------------------------
 *  iskeptype  --  checks the validity of a type of kepler orbit
 *		   accepts: type, a type of kepler orbit
 *		   returns: TRUE is type is a valid type, FALSE if not.
 *-----------------------------------------------------------------------------
 */
#define   iskeptype(type)	(type == CARTESIAN || type == SPHERICAL	      \
				|| type == TRUEANOMALY || type == ECCANOMALY  \
				|| type == MEANANOMALY || type == PERIPASSAGE \
				|| type == SCATTERING || type == TSCATTERING)

/*-----------------------------------------------------------------------------
 *  macros which limit the degrees of freedom in a gravitational N-body
 *  experiment. Since the gravitational constant is already taken to be unity,
 *  (G = 1), we are left with one mass, position and velocity per particle.
 *  Therefore, in three dimensions we have:
 * 
 *     Number of degrees of freedom:  7N
 *
 *  The following considerations limit this freedom:
 *  I): We can chose physical units of length, mass, and time. However,
 *      this gives us only two constraints, the third one having been used
 *      already to set  G = 1 .
 *  II): homogeneity of spacetime gives us seven additional constraints, 
 *       one for time, three for position and three for velocity.
 *  III): isotropy of space gives us three more constraints.
 *  These twelve constraints leave us with:
 *
 *     Number of independent degrees of freedom:  7N - 12  =  2 + 7(N-2)
 *
 *  Since the target and the projectile consist of either single stars or
 *  hierarchical binary systems, the number of internal binaries is  N - 2 .
 *  Thus the number of independent degrees of freedom is two to start out with,
 *  and in addition seven per internal binary.
 *  note: the starting number two can be associated with:
 *            1) the mass ratio of the target and projectile;
 *            2) the type of incoming hyperbola (e.g. the eccentricity).
 *        although we use a different convention below, where the freedom of
 *        choice is used to fix the target internal binary orbit, instead of
 *        the relative orbit of projectile and target.
 *-----------------------------------------------------------------------------
 */

/*
 *  I) Freedom of scaling of physical units.
 *
 *     scaling of mass:
 */
#define  TARGET_MASS                 1.0
/*
 *     scaling of length:
 */
#define  TARGET_SEMIMAJOR_AXIS       1.0
/*
 *     scaling of time:
 *
 * the convention  G = 1 , with the above two, already fixes the time unit
 *
 *  II) homogeneity of spacetime:
 *     
 *     of time:
 */    
#define  TARGET_TIME_OF_PERICENTER_PASSAGE   0.0
/*
 *     of positions:
 */
#define  TARGET_COM_RX    0.0
#define  TARGET_COM_RY    0.0
#define  TARGET_COM_RZ    0.0
/*
 *     of velocities:
 */
#define  TARGET_COM_VX    0.0
#define  TARGET_COM_VY    0.0
#define  TARGET_COM_VZ    0.0
/* 
 *  III) isotropy of space
 */
#define  TARGET_INCLINATION                   0.0
#define  TARGET_LONGITUDE_OF_ASCENDING_NODE   0.0
#define  TARGET_ARGUMENT_OF_PERICENTER        0.0

/* endof: scat.h */
