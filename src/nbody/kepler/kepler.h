/*  kepler.h - ASCII, HUMAN, BINARY, CARTESIAN, ECCANOMALY, KEPLER, KEPLERPTR,
	      KArgper, KAsemi, KEcc, KEccan, KImp_angle, KImp_param,
	      KImp_separ, KIncl, KLonasc, KMass1, KMass2, KMeanan, KPerpas,
	      Kphi, KPhsp, KRmod, KRvect, KSpace, KTheta, KTime, KTruean,
	      KType, KVphi, KVinf_abs, KVinf_theta, KVinf_phi, KVrad, KVtheta,
	      KVvect, MEANANOMALY, NBODY, NBODYPTR, NCom, NList, NNumber,
	      NPart, NTime, PARTICLE, PARTPTR, PERIPASSAGE, PHSPACE, PHSPPTR,
	      PMass, PPhsp, PRvect, PSpace, PType, PVvect, SCATTERING,
	      SPHERICAL, SRvect, SSpace, SType, SVvect, TRUEANOMALY,
	      TSCATTERING, iskeptype */
/*
 *  kepler.h: header file for manipulations with nbody systems
 *		(c) 1986  Piet Hut  Princeton, NJ, USA
 *
 *           $!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$
 *	     $								 $
 *	     $  BEWARE: 3-dim. Cartesian coordinates are numbered as	 $
 *	     $								 $
 *	     $	            { z, x, y }  =  { 0, 1, 2 }  .		 $
 *	     $								 $
 *	     $		Using other conventions is likely to lead to	 $
 *	     $		errors in invoking coordinate transformations!!  $
 *	     $								 $
 *	     $   BEWARE!!       BEWARE!!       BEWARE!!       BEWARE!!   $
 *	     $								 $
 *           $!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$!@#%$
 */
#include <math.h>
#include <stdinc.h>
#include <vectmath.h>

/*-----------------------------------------------------------------------------
 *  PHSPACE, PHSPPTR  --  define a structure of type phasespace,
 *			  which contains the phasespace coordinates,
 *			  together with an indicater of the type of
 *			  coordinate system used.
 *-----------------------------------------------------------------------------
 */
typedef  struct
    {
    int  phsp_type;           /* type of coordinate representation         */
    real  phsp_coord[2*NDIM]; /* free form for the six coordinates         */
    } PHSPACE, *PHSPPTR;

/*-----------------------------------------------------------------------------
 *  PARTICLE, PARTPTR  --  define a structure of type particle,
 *			   which contains its mass and its phasespace position
 *-----------------------------------------------------------------------------
 */
typedef  struct
    {
    real  part_mass;          /* mass                                      */
    PHSPACE  part_space;      /* position in phase space                   */
    } PARTICLE, *PARTPTR;

/*-----------------------------------------------------------------------------
 *  NBODY, NBODYPTR  --  define a structure of type nbody system,
 *			 which contains: the total particle number;
 *					 the time of the snapshot;
 *					 the center-of-mass information
 *					     in the form of a pseudo-particle;
 *					 an array of particles.
 *			 note: the center-of-mass (c.o.m.) information
 *			       includes the total mass of the system as well as
 *			       the position and velocity of the c.o.m.
 *			 note: the last entry is a pointer to the first
 *			       particle; an array of particle pointers may
 *			       seem more in line with the previous definitions,
 *			       being a pointer to a particle pointer, but is
 *			       less practical when reserving storage; the
 *			       implementation allows block storage allocation
 *			       using calloc(), cf.  mknbody().
 *-----------------------------------------------------------------------------
 */
typedef  struct
    {
    int  nbody_n;             /* number of particles                       */
    real  nbody_time;         /* time                                      */
    PARTPTR  nbody_com;       /* pointer to the center of mass             */
    PARTPTR  nbody_list;      /* pointer to the first particle in an array */
    } NBODY, *NBODYPTR;

/*-----------------------------------------------------------------------------
 *  KEPLER, KEPLERPTR  --  define a structure of type kepler orbit.
 *			   which contains: the time of the snapshot;
 *					   the masses of the two bodies;
 *					   the phase space variable describing
 *					       their relative orbit.
 *			   note: many types of coordinate system are supported,
 *				 including several types of orbital elements,
 *				 which require a knowledge of the time
 *				 (to determine the time of pericenter passage)
 *				 and of the total mass (to determine the
 *				 the semimajor axis from pos. and vel.).
 *			   note: in the present implementation both masses
 *				 are included which is redundant, but is useful
 *				 when transforming to the orbits of the
 *				 individual bodies with respect to the c.o.m.
 *-----------------------------------------------------------------------------
 */
typedef  struct
    {
    real  kep_time;           /* time at which this representation applies */
    real  kep_m1, kep_m2;     /* masses                                    */
    PHSPACE  kep_space;       /* relative position in phase space          */
    } KEPLER, *KEPLERPTR;

/*-----------------------------------------------------------------------------
 *  Macros to extract individual components from a phase space structure
 *-----------------------------------------------------------------------------
 */
#define   SType(phspptr)             ((phspptr)->phsp_type)
#define   SSpace(phspptr)            (&(phspptr)->phsp_coord[0])
#define   SRvect(phspptr)            (&(phspptr)->phsp_coord[0])
#define   SVvect(phspptr)            (&(phspptr)->phsp_coord[3])

/*-----------------------------------------------------------------------------
 *  Macros to extract individual components from a particle
 *-----------------------------------------------------------------------------
 */
#define   PMass(partptr)             ((partptr)->part_mass)
#define   PPhsp(partptr)             ((partptr)->part_space.phsp_coord)
#define   PType(partptr)             ((partptr)->part_space.phsp_type)
#define   PSpace(partptr)            (&(partptr)->part_space.phsp_coord[0])
#define   PRvect(partptr)            (&(partptr)->part_space.phsp_coord[0])
#define   PVvect(partptr)            (&(partptr)->part_space.phsp_coord[3])

/*-----------------------------------------------------------------------------
 *  Macros to extract individual components from a nbody system
 *-----------------------------------------------------------------------------
 */
#define   NNumber(nbodyptr)          ((nbodyptr)->nbody_n)
#define   NTime(nbodyptr)            ((nbodyptr)->nbody_time)
#define   NCom(nbodyptr)             ((nbodyptr)->nbody_com)
#define   NList(nbodyptr)            ((nbodyptr)->nbody_list)
#define   NPart(nbodyptr, i)         (&(nbodyptr)->nbody_list[(i)])

     /*
      *  example:  to obtain the Cartesian coordinate 'vy'
      *            of the velocity vector of particle 'j' in the
      *            NBODY structure   *p  :
      *            type:
      *                   vy = PVvect(NPart(p,j))[2]
      *            since:
      *               coordinates are numbered {z,x,y} = {0,1,2}
      *               particles are numbered {0,1,2,...,j,...}
      *  note:  beware!
      *         to find 'vy' in this way, you'd better check first 
      *         whether the following two conditions are satisfied:
      *
      *         NNumber(x) > j  
      *         PType(NPart(x,j)) == CARTESIAN
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
				|| type == SCATTERING  || type == TSCATTERING)

/*-----------------------------------------------------------------------------
 *  Macros to extract individual components from a kepler orbit '*kepptr'
 *-----------------------------------------------------------------------------
 */
     /*  The two masses, time, and type of phase space coordinates:  */

#define   KTime(kepptr)              ((kepptr)->kep_time)
#define   KMass1(kepptr)   	     ((kepptr)->kep_m1)
#define   KMass2(kepptr)   	     ((kepptr)->kep_m2)
#define   KPhsp(kepptr)              ((kepptr)->kep_space.phsp_coord)
#define   KType(kepptr)    	     ((kepptr)->kep_space.phsp_type)

     /*  Alphabetic list of phase space coordinates    */
     /*  Note: the radius vector points from m1 to m2  */

#define   KArgper(kepptr)            ((kepptr)->kep_space.phsp_coord[4])
#define   KAsemi(kepptr)             ((kepptr)->kep_space.phsp_coord[0])
#define   KEcc(kepptr)               ((kepptr)->kep_space.phsp_coord[1])
#define   KEccan(kepptr)             ((kepptr)->kep_space.phsp_coord[5])
#define   KImp_angle(kepptr)         ((kepptr)->kep_space.phsp_coord[1])
#define   KImp_param(kepptr)         ((kepptr)->kep_space.phsp_coord[0])
#define   KImp_separ(kepptr)         ((kepptr)->kep_space.phsp_coord[2])
#define   KImp_time(kepptr)          ((kepptr)->kep_space.phsp_coord[2])
#define   KIncl(kepptr)              ((kepptr)->kep_space.phsp_coord[2])
#define   KLonasc(kepptr)            ((kepptr)->kep_space.phsp_coord[3])
#define   KMeanan(kepptr)            ((kepptr)->kep_space.phsp_coord[5])
#define   KPerpas(kepptr)            ((kepptr)->kep_space.phsp_coord[5])
#define   KPhi(kepptr)               ((kepptr)->kep_space.phsp_coord[2])
#define   KRmod(kepptr)              ((kepptr)->kep_space.phsp_coord[0])
#define   KRvect(kepptr)             (&(kepptr)->kep_space.phsp_coord[0])
#define   KSpace(kepptr)             (&(kepptr)->kep_space.phsp_coord[0])
#define   KTheta(kepptr)             ((kepptr)->kep_space.phsp_coord[1])
#define   KTruean(kepptr)            ((kepptr)->kep_space.phsp_coord[5])
#define   KVphi(kepptr)              ((kepptr)->kep_space.phsp_coord[5])
#define   KVinf_abs(kepptr)          ((kepptr)->kep_space.phsp_coord[3])
#define   KVinf_theta(kepptr)        ((kepptr)->kep_space.phsp_coord[4])
#define   KVinf_phi(kepptr)          ((kepptr)->kep_space.phsp_coord[5])
#define   KVrad(kepptr)              ((kepptr)->kep_space.phsp_coord[3])
#define   KVtheta(kepptr)            ((kepptr)->kep_space.phsp_coord[4])
#define   KVvect(kepptr)             (&(kepptr)->kep_space.phsp_coord[3])

/*-----------------------------------------------------------------------------
 *  List of phase space coordinates names, ordered by coordinate
 *.............................................................................
 *   array names:
 *
 *        KSpace(kepptr)             (&(kepptr)->kep_space.phsp_coord[0])
 *        KRvect(kepptr)             (&(kepptr)->kep_space.phsp_coord[0])
 *        KVvect(kepptr)             (&(kepptr)->kep_space.phsp_coord[3])
 *
 *   scalar names:
 *
 *        KRmod(kepptr)              ((kepptr)->kep_space.phsp_coord[0])
 *        KAsemi(kepptr)             ((kepptr)->kep_space.phsp_coord[0])
 *        KImp_param(kepptr)         ((kepptr)->kep_space.phsp_coord[0])
 *
 *        KTheta(kepptr)             ((kepptr)->kep_space.phsp_coord[1])
 *        KEcc(kepptr)               ((kepptr)->kep_space.phsp_coord[1])
 *        KImp_angle(kepptr)         ((kepptr)->kep_space.phsp_coord[1])
 *
 *        KPhi(kepptr)               ((kepptr)->kep_space.phsp_coord[2])
 *        KIncl(kepptr)              ((kepptr)->kep_space.phsp_coord[2])
 *        KImp_separ(kepptr)         ((kepptr)->kep_space.phsp_coord[2])
 *        KImp_time(kepptr)          ((kepptr)->kep_space.phsp_coord[2])
 *
 *        KVrad(kepptr)              ((kepptr)->kep_space.phsp_coord[3])
 *        KLonasc(kepptr)            ((kepptr)->kep_space.phsp_coord[3])
 *        KVinf_abs(kepptr)          ((kepptr)->kep_space.phsp_coord[3])
 *
 *        KVtheta(kepptr)            ((kepptr)->kep_space.phsp_coord[4])
 *        KArgper(kepptr)            ((kepptr)->kep_space.phsp_coord[4])
 *        KVinf_theta(kepptr)        ((kepptr)->kep_space.phsp_coord[4])
 *
 *        KVphi(kepptr)              ((kepptr)->kep_space.phsp_coord[5])
 *        KTruean(kepptr)            ((kepptr)->kep_space.phsp_coord[5])
 *        KEccan(kepptr)             ((kepptr)->kep_space.phsp_coord[5])
 *        KMeanan(kepptr)            ((kepptr)->kep_space.phsp_coord[5])
 *        KPerpas(kepptr)            ((kepptr)->kep_space.phsp_coord[5])
 *        KVinf_phi(kepptr)          ((kepptr)->kep_space.phsp_coord[5])
 *

 *.............................................................................
 *
 *   List of phase space coordinates names, ordered by type:
 *.............................................................................
 *   general:
 *
 *        KSpace(kepptr)             (&(kepptr)->kep_space.phsp_coord[0])
 *
 *   CARTESIAN:
 *
 *        KRvect(kepptr)             (&(kepptr)->kep_space.phsp_coord[0])
 *        KVvect(kepptr)             (&(kepptr)->kep_space.phsp_coord[3])
 *
 *   SPHERICAL:
 *
 *        KRmod(kepptr)              ((kepptr)->kep_space.phsp_coord[0])
 *        KTheta(kepptr)             ((kepptr)->kep_space.phsp_coord[1])
 *        KPhi(kepptr)               ((kepptr)->kep_space.phsp_coord[2])
 *        KVrad(kepptr)              ((kepptr)->kep_space.phsp_coord[3])
 *        KVtheta(kepptr)            ((kepptr)->kep_space.phsp_coord[4])
 *        KVphi(kepptr)              ((kepptr)->kep_space.phsp_coord[5])
 *
 *   TRUEANOMALY:
 *
 *        KAsemi(kepptr)             ((kepptr)->kep_space.phsp_coord[0])
 *        KEcc(kepptr)               ((kepptr)->kep_space.phsp_coord[1])
 *        KIncl(kepptr)              ((kepptr)->kep_space.phsp_coord[2])
 *        KLonasc(kepptr)            ((kepptr)->kep_space.phsp_coord[3])
 *        KArgper(kepptr)            ((kepptr)->kep_space.phsp_coord[4])
 *        KTruean(kepptr)            ((kepptr)->kep_space.phsp_coord[5])
 *
 *   ECCANOMALY:
 *
 *        KAsemi(kepptr)             ((kepptr)->kep_space.phsp_coord[0])
 *        KEcc(kepptr)               ((kepptr)->kep_space.phsp_coord[1])
 *        KIncl(kepptr)              ((kepptr)->kep_space.phsp_coord[2])
 *        KLonasc(kepptr)            ((kepptr)->kep_space.phsp_coord[3])
 *        KArgper(kepptr)            ((kepptr)->kep_space.phsp_coord[4])
 *        KEccan(kepptr)             ((kepptr)->kep_space.phsp_coord[5])
 *
 *   MEANANOMALY:
 *
 *        KAsemi(kepptr)             ((kepptr)->kep_space.phsp_coord[0])
 *        KEcc(kepptr)               ((kepptr)->kep_space.phsp_coord[1])
 *        KIncl(kepptr)              ((kepptr)->kep_space.phsp_coord[2])
 *        KLonasc(kepptr)            ((kepptr)->kep_space.phsp_coord[3])
 *        KArgper(kepptr)            ((kepptr)->kep_space.phsp_coord[4])
 *        KMeanan(kepptr)            ((kepptr)->kep_space.phsp_coord[5])
 *
 *   PERIPASSAGE:
 *
 *        KAsemi(kepptr)             ((kepptr)->kep_space.phsp_coord[0])
 *        KEcc(kepptr)               ((kepptr)->kep_space.phsp_coord[1])
 *        KIncl(kepptr)              ((kepptr)->kep_space.phsp_coord[2])
 *        KLonasc(kepptr)            ((kepptr)->kep_space.phsp_coord[3])
 *        KArgper(kepptr)            ((kepptr)->kep_space.phsp_coord[4])
 *        KPerpas(kepptr)            ((kepptr)->kep_space.phsp_coord[5])
 *

 *   SCATTERING:
 *
 *        KImp_param(kepptr)         ((kepptr)->kep_space.phsp_coord[0])
 *
 *            the impact parameter, which is the offset at infinity from
 *            a head-on collision, i.e. the closest distance to the origin
 *            if the incoming particle would move in a straight line,
 *            starting from an infinite separation, instead of following
 *            a hyperbolic orbit; in other words, the impact parameter is
 *            the closest distance to the origin of the asymptote to the
 *            hyperbolic orbit.
 *
 *        KImp_angle(kepptr)         ((kepptr)->kep_space.phsp_coord[1])
 *
 *            for an incoming particle, the impact angle is defined as
 *            the angle from a fiducial vector to the impact vector,
 *            in counter-clockwise direction as seen from the origin.
 * 
 *            the fiducial vector is a unit vector in the plane spanned
 *            by the z-axis and the velocity vector at infinity,
 *            perpendicular to the latter and with a positive z-component,
 *            i.e. pointing upwards.
 *
 *            the impact vector is a vector perpendicular to the velocity
 *            at infinity, with a length equal to the impact parameter,
 *            pointing in the direction of the offset with respect to 
 *            a head-on collision.
 *              a constructive definition for the impact vector:
 *              introduce an extra incoming particle with the same
 *              asymptotic velocity (velocity at infinity) but on a
 *              head-on collision orbit. the impact vector is the vector
 *              from the position of this particle to the position of the
 *              real incoming particle positioned at such a point on its
 *              orbit that this vector is perpendicular to their common
 *              velocity vector, and very far away from the origin.
 *            
 *            in the case of an outgoing orbit, the impact vector is similarly
 *            defined as the offset with respect to a collision orbit at very
 *            late times after the scattering event. also in this case the
 *            fiducial vector points upwards, and the impact angle is counted
 *            counter-clockwise as seen from the origin.
 * 
 *        KImp_separ(kepptr)         ((kepptr)->kep_space.phsp_coord[2])
 *
 *            the distance between the two particles.
 *            note that for any finite distance, the relative velocity
 *            between the two particles is larger than the asymptotic value
 *            of their velocities,  KVinf_abs  .
 *

 *        KVinf_abs(kepptr)          ((kepptr)->kep_space.phsp_coord[3])
 *
 *            for ingoing orbits:
 *            the absolute value of the relative velocity at infinity,
 *            i.e. in the limit that the two particles are infinitely
 *            far removed and therefore approach the asymptote of their
 *            relative kepler orbit.
 *            NOTE: a trick is used to discriminate between
 *                  ingoing and outgoing orbits: outgoing orbits
 *                  are indicated by adding a minus sign to  KVinf_abs() ,
 *                  while ingoing orbits have the usual plus sign.
 *            BEWARE !!
 *
 *        KVinf_theta(kepptr)        ((kepptr)->kep_space.phsp_coord[4])
 *
 *            theta and phi are the usual spherical angles as used in
 *            spherical coordinates, and are used here to point in the
 *            direction of the ingoing or outgoing particle 
 *            (as the case may be) as seen from the origin at asymptotically
 *            early or late times.
 *            thus for an incoming orbit, the angles point in the direction
 *            exactly opposite to the velocity at time t --> -infinity,
 *            while for an outgoing orbit, the angles point in the same
 *            direction as the velocity at time t --> +infinity.
 *
 *        KVinf_phi(kepptr)          ((kepptr)->kep_space.phsp_coord[5])
 *
 *            see above.
 *

 *   TSCATTERING:
 *  
 *        KImp_param(kepptr)         ((kepptr)->kep_space.phsp_coord[0])
 *
 *            as in SCATTERING.
 *
 *        KImp_angle(kepptr)         ((kepptr)->kep_space.phsp_coord[1])
 *
 *            as in SCATTERING.
 *
 *        KImp_time(kepptr)          ((kepptr)->kep_space.phsp_coord[2])
 *
 *            the time before pericenter passage (negative
 *            for an outgoing orbit)
 *
 *        KVinf_abs(kepptr)          ((kepptr)->kep_space.phsp_coord[3])
 *
 *            as in SCATTERING, with the fortunate difference that there is
 *            no need for the hack which used a negative value to indicate
 *            outgoing orbits. In TSCATTERING, we have only positive values
 *            for the absolute value of the relative velocity at infinity.
 *
 *        KVinf_theta(kepptr)        ((kepptr)->kep_space.phsp_coord[4])
 *
 *            as in SCATTERING.
 *
 *        KVinf_phi(kepptr)          ((kepptr)->kep_space.phsp_coord[5])
 *
 *            as in SCATTERING.
 *
 */

/* end of: kepler.h */
