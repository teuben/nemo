/* transformext.c - */

/*
 *  transformext.c:  for transforming to & from polynomial orbit extrapolation
 *
 *      Feb 1988  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
#include  "newton0.h"

static orbitsegmentptr mk_orbitsegments(int npart);
static void nonext_to_ext(bodyptr bodies, orbitsegmentptr extbodies, int npart);
static void lay_the_tracks(orbitsegmentptr extbodies, int npart, real t_now, real eta);
static void firstpass(orbitsegmentptr pi, orbitsegmentptr pj);
static void secondpass(orbitsegmentptr pi, orbitsegmentptr pj);
static void thirdpass(orbitsegmentptr pi, real t_now, real eta);
static void synchronize_particles(orbitsegmentptr extbodies, int npart, real sync_time);
static void get_pos_vel(bodyptr bodies, orbitsegmentptr extbodies, int npart);

/*-----------------------------------------------------------------------------
 *  setup_ext  --  sets up an orbit EXTrapolation representation for the
 *                 system  initial_sys , in three steps:
 *                 1) memory allocation for polynomial orbit approximation;
 *                 2) reading in of masses, positions, velocities;
 *                 3) initialization of polynomial orbit approximation.
 *                 Litt.: Chap. 2 of S.J. Aarseth, 
 *                          in Multiple Time Scales, eds. J.U. Brackbill
 *                          and B.I. Cohen (Academic Press, 1985), p. 377.
 *                        Appendix 4.B of J. Binney and S. Tremaine,
 *	                    Galactic Dynamics (Princeton Univ. Press, 1987),
 *	                    which contains a Fortran version of a simple
 *                          Aarseth code, from which the present version grew.
 *-----------------------------------------------------------------------------
 */
void  setup_ext(initial_sys, initial_specs)
systptr  initial_sys;
specptr  initial_specs;
    {
    int  npart;
    bodyptr  bodies;
    orbitsegmentptr  extbodies;
    real  t_now;           /* current time */
    real  eta;             /* dimensionless integration accuracy parameter,  */
                           /* used in determining the size of the time steps */

    t_now = Tnow(initial_sys);
    eta = Stepparam(initial_specs);

    npart = Nbody(initial_sys);
    bodies = Bodies(initial_sys);
    extbodies = Orbitsegments(initial_sys) = mk_orbitsegments(npart);

    nonext_to_ext(bodies, extbodies, npart);

    lay_the_tracks(extbodies, npart, t_now, eta);
    }

/*-----------------------------------------------------------------------------
 *  read_ext  --  reads the relevant contents of an orbit EXTrapolation
 *                representation and copies those back into the original
 *                 newton0  representation.
 *                This requires three steps:
 *                1) synchronize all particles in an orbitsegment
 *                   representation;
 *                2) copy the positions and velocities of all particles
 *                   back into the original  newton0 representation.
 *                3) compute the accelerations and potentials for all
 *                   particles at the shared new time, by direct summation
 *                   over all particle pair interactions.
 *                Litt.: see above, in  setup_ext() .
 *-----------------------------------------------------------------------------
 */
void  read_ext(sys, specs, readout_time)
systptr  sys;
specptr  specs;
real  readout_time;
    {
    int  npart;
    bodyptr  bodies;
    orbitsegmentptr  extbodies;

    npart = Nbody(sys);
    bodies = Bodies(sys);
    extbodies = Orbitsegments(sys);

    synchronize_particles(extbodies, npart, readout_time);
    get_pos_vel(bodies, extbodies, npart);
    deriv_system(sys, specs);   /* in file  differentiate.c  */
    }

/*-----------------------------------------------------------------------------
 *  mk_orbitsegments  --  constructs an array of new particles in the
 *                        representation of orbitsegments.
 *                        accepts: npart: the number of particles.
 *                        returns: pointer to a contiguus block of memory of
 *                                 the right size to contain an array of npart
 *                                 particles in orbitsegment nrepresentation.
 *		          initialization: all entrees are initialized to zero,
 *			                  using calloc() rather than malloc().
 *-----------------------------------------------------------------------------
 */
local orbitsegmentptr  mk_orbitsegments(npart)
int  npart;
    {
    orbitsegmentptr  extbodies;

    extbodies = (orbitsegmentptr) calloc((unsigned)npart, sizeof(orbitsegment));
    if (extbodies == NULL)
	error("mk_orbitsegments: too little memory left for %d particles\n",
                                                                       npart);
    return(extbodies);
    }

/*-----------------------------------------------------------------------------
 *  nonext_to_ext  --  initializes masses, positions and velocities
 *                     for  extbodies[npart]  from  bodies[npart] .
 *-----------------------------------------------------------------------------
 */
local void  nonext_to_ext(bodies, extbodies, npart)
bodyptr  bodies;
orbitsegmentptr  extbodies;
int  npart;
    {
    bodyptr  body_i;      /* points to an individual particle of bodies[] */
    orbitsegmentptr  extbody_i;       /* points to an extbodies[] element */

    for (body_i = bodies, extbody_i = extbodies; body_i - bodies < npart;
	 body_i++, extbody_i++)
        {
        AMass(extbody_i) = Mass(body_i);
	SETV(APos(extbody_i), Pos(body_i));
	SETV(AVel(extbody_i), Vel(body_i));
	}
    }

/*-----------------------------------------------------------------------------
 *  lay_the_tracks  --  starts up an Aarseth type polynomial orbit
 *                      extrapolation representation for the array of
 *                      particles contained in  extbodies[npart] ,
 *                      in three separate passes, each of which loops
 *                      over all particles.
 *-----------------------------------------------------------------------------
 */
local void  lay_the_tracks(extbodies, npart, t_now, eta)
orbitsegmentptr  extbodies;
int  npart;
real  t_now;
real  eta;
    {
    orbitsegmentptr  pi;
    int  i, j, k;

   /* initialize force and first derivative between i and j */

    for (i = 0, pi = &extbodies[0]; i < npart; i++, pi++)
        for (k = 0; k < NDIM; k++)
	    AForce(pi)[k] = AJerk(pi)[k] = 0;
    for (i = 0, pi = &extbodies[0]; i < npart; i++, pi++)
        for (j = 0; j < npart; j++)
            if (i != j)
                firstpass(&extbodies[i], &extbodies[j]);

   /* calculate second and third derivatives of forces between i, j */

    for (i = 0, pi = &extbodies[0]; i < npart; i++, pi++)
        for (k = 0; k < NDIM; k++)
	    D2(pi)[k] = D3(pi)[k] = 0;
    for (i = 0; i < npart; i++)
        for (j = 0; j < npart; j++)
            if (i != j)
                secondpass(&extbodies[i], &extbodies[j]);

   /* initialize integration steps and convert to force differences */

    for (i = 0; i < npart; i++)
        thirdpass(&extbodies[i], t_now, eta);
    }

/*-----------------------------------------------------------------------------
 *  firstpass  --  determine the forces and jerks on all particles,
 *                 on a first double pass over the whole system
 *-----------------------------------------------------------------------------
 */
local void  firstpass(pi, pj)
orbitsegmentptr  pi, pj;	/* pointers to the i'th and j'th particle */
    {
    int  k;			/* counter */
    real  rij[NDIM];		/* difference of position of i and j */
    real  vij[NDIM];		/* difference of velocities of i and j */
    real  pfij[NDIM];		/* force between the pair of particles i,j */
    real  rinv2;		/* inverse distance between i and j squared */
    real  mrinv3;		/* mass divided by |rij| cubed */
    real  a;			/* as defined by Aarseth (1985) p. 382 */

   /* find differences in position and velocity */

    for (k = 0; k < NDIM; k++)
	{
        rij[k] = APos(pi)[k] - APos(pj)[k];
        vij[k] = AVel(pi)[k] - AVel(pj)[k];
	}
    rinv2 = 1.0 / dotvp(rij, rij);
    mrinv3 = AMass(pj) * rinv2 * sqrt(rinv2);
    a = dotvp(rij, vij) * rinv2;

   /* calculate force and jerk using eq.(5), Aarseth (1985) */

    for (k = 0; k < NDIM; k++)
	{
        pfij[k] = -rij[k] * mrinv3;
        AForce(pi)[k] += pfij[k];
        AJerk(pi)[k]  -= vij[k] * mrinv3 + 3 * a * pfij[k];
	}
    }

/*-----------------------------------------------------------------------------
 *  secondpass  --  determine the higher time derivatives of the forces
 *                  on a second double pass over the whole system
 *-----------------------------------------------------------------------------
 */
local void  secondpass(pi, pj)
orbitsegmentptr  pi, pj;	/* pointers to the i'th and j'th particle */
    {
    int  k;			/* counter */
    real  rij[NDIM];		/* difference in position */
    real  vij[NDIM];		/* difference in velocity */
    real  fij[NDIM];		/* difference in force */
    real  f1ij[NDIM];		/* difference in jerk (d/dt of force) */

   /* NOTE: the previous two definitions give the difference between total
    *       force and jerk of each particle, as exerted by all other particles.
    *       The following definitions pertain to the single mutual force
    *       between only one pair of particles, and the time derivatives of
    *       that single quantity.
    * BEWARE: In Aarseth's (1985) notation:
    *         our  "fij"  is his  "F sub i - F sub j" ;	   ===>	GLOBAL !
    *         our  "pfij"  is his  "F sub ij".		   ===>	LOCAL !
    */
    real  pfij[NDIM];		/* force between the pair of particles i,j */
    real  pf1ij[NDIM];		/* first time derivative of pfij[] */
    real  pf2ij[NDIM];		/* second time derivative of pfij[] */
    real  pf3ij[NDIM];		/* third time derivative of pfij[] */
    real  a;			/* as defined by Aarseth (1985) p. 382 */
    real  b;			/* as defined by Aarseth (1985) p. 383 */
    real  c;			/* as defined by Aarseth (1985) p. 383 */
    real  rinv2;		/* inverse distance between i and j squared */
    real  mrinv3;		/* mass divided by |rij| cubed */

   /* calculate differences of global quantities between particles i and j */

    for (k = 0; k < NDIM; k++)
	{
        rij[k] = APos(pi)[k] - APos(pj)[k];
        vij[k] = AVel(pi)[k] - AVel(pj)[k];
        fij[k] = AForce(pi)[k] - AForce(pj)[k];
        f1ij[k] = AJerk(pi)[k] - AJerk(pj)[k];
	}
    rinv2 = 1 / dotvp(rij, rij);
    mrinv3 = AMass(pj) * rinv2 * sqrt(rinv2);
    a = dotvp(rij, vij) * rinv2;
    b = (dotvp(vij, vij) + dotvp(rij, fij)) * rinv2 + a * a;
    c = (3 * dotvp(vij, fij) + dotvp(rij, f1ij)) * rinv2 
        + a * (3 * b - 4 * a * a);

   /*  calculate mutual force between a pair of particles,
    *  and its first three time derivatives using eq.(6), Aarseth (1985)
    */
    for (k = 0; k < NDIM; k++)
	{
	pfij[k]  = -mrinv3 * rij[k];
	pf1ij[k] = -mrinv3 * vij[k] - 3 * a * pfij[k];
	pf2ij[k] = -mrinv3 * fij[k] - 6 * a * pf1ij[k] - 3 * b * pfij[k];
	pf3ij[k] = -mrinv3 * f1ij[k] - 9 * a * pf2ij[k] - 9 * b * pf1ij[k]
							- 3 * c * pfij[k];
       /* NOTE: alas, there are no memory allocations in which to accumulate
        *       and store the pf2 and pf3 components (because this would
	*       significantly increase the total amount of memory required).
	*       Therefore, a hack is introduced:
	*       we use the D2 and D3 arrays instead, but only from here until
	*       the beginning the next function (thirdpass).
	*/
        D2(pi)[k] += pf2ij[k];  /* NOTE: secondary usage of memory location! */
        D3(pi)[k] += pf3ij[k];  /* NOTE: secondary usage of memory location! */
	}
    }

/*-----------------------------------------------------------------------------
 *  thirdpass  --  determine individual time steps and
 *                 initialize the force polynomials
 *                 on a final single pass over the whole system
 *-----------------------------------------------------------------------------
 */
local void  thirdpass(pi, t_now, eta)
orbitsegmentptr  pi;		  /* pointer to the i'th particle */
real  t_now;
real  eta;
    {
    int  k;
    real  f2i[NDIM];      /* second time derivative of force on ith particle */
    real  f3i[NDIM];      /* third time derivative of force on ith particle  */
    real  fi_square;	  /* absolute value squared of force on ith particle */
    real  f1i_square;	  /* abs.v.sq. of 1st time derivative of force on i  */
    real  f2i_square;	  /* abs.v.sq. of 2nd time derivative of force on i  */
    real  f3i_square;	  /* abs.v.sq. of 3rd time derivative of force on i  */
    real  tstep_sq;	  /* square of initial time step scale               */

   /* NOTE: we first have to undo the secondary usage of the D2 and D3 arrays
    *       introduced at the end of the previous function (secondpass),
    *       so that they can be used to store the polynomial coefficients
    *       for which they were designed in the first place,
    *       at the end of this function.
    */
    for (k = 0; k < NDIM; k++)
	{
	f2i[k] = D2(pi)[k];     /* NOTE: secondary usage of memory location! */
	f3i[k] = D3(pi)[k];     /* NOTE: secondary usage of memory location! */
	}

   /* initialize time steps using eq.(9), Aarseth (1985) */

    fi_square = dotvp(AForce(pi), AForce(pi));
    f1i_square = dotvp(AJerk(pi), AJerk(pi));
    f2i_square = dotvp(f2i, f2i);
    f3i_square = dotvp(f3i, f3i);
   /*
    * note: do not multiply the ..._square's below before taking the sqrt;
    *       they are likely to cause overflow on a VAX
    */
    tstep_sq = (sqrt(fi_square) * sqrt(f2i_square) + f1i_square)
	       / (sqrt(f1i_square) * sqrt(f3i_square) + f2i_square);
    ATimestep(pi) = eta * sqrt(tstep_sq);      /* Aarseth has eta under sqrt */

   /*  initialize the roots of the force polynomial,
    *  i.e. the backward times needed to construct the polynomial
    */
    T0(pi) = t_now;
    T1(pi) = t_now - ATimestep(pi);
    T2(pi) = t_now - 2*ATimestep(pi);
    T3(pi) = t_now - 3*ATimestep(pi);

   /*  initialize coefficients of the force polynomial
    *  using eq.(7), Aarseth (1985)
    */
    for (k = 0; k < NDIM; k++)
	{
        D1(pi)[k] = (f3i[k] * ATimestep(pi)/6 - f2i[k]/2) * ATimestep(pi)
		   + AJerk(pi)[k];
        D2(pi)[k] = (f2i[k] - f3i[k] * ATimestep(pi)) / 2;
        D3(pi)[k] = f3i[k] / 6;
	}
    }

/*-----------------------------------------------------------------------------
 *  synchronize_particles  --  push all particles forward in time along their
 *                             orbit segments, until a shared time  sync_time .
 *-----------------------------------------------------------------------------
 */
local void  synchronize_particles(extbodies, npart, sync_time)
orbitsegmentptr  extbodies;
int  npart;
real  sync_time;
    {
    int  i, k;
    real  dt;
    real  f2i[NDIM];
    real  f3i[NDIM];
    orbitsegmentptr  pi;

    for (pi = extbodies, i = 0; i < npart; pi++, i++)
	{
        dt = sync_time - T0(pi);
       /*
        * predict the position and velocity of a particle
        * to higher order accuracy
	* (up to the third time derivative of the force)
	*/
        for (k = 0; k < NDIM; k++)
	    {
           /*
            *  find higher order force derivatives using eq.(4), Aarseth (1985)
            */
            f2i[k] = 2 * ( D3(pi)[k] * ((T0(pi) - T1(pi)) + (T0(pi) - T2(pi)))
		           + D2(pi)[k] );
	    f3i[k] = 6 * D3(pi)[k];
           /*
	    *  integrate the equation of motion
            *  using the higher order force approximation
	    */
            APos_now(pi)[k] = (((((f3i[k]/120) * dt + f2i[k]/24) * dt
			      + AJerk(pi)[k]/6) * dt + AForce(pi)[k]/2) * dt
			      + AVel(pi)[k])*dt + APos(pi)[k];
            AVel_now(pi)[k] = ((((f3i[k]/24) * dt + f2i[k]/6)*dt
			      + AJerk(pi)[k]/2) * dt + AForce(pi)[k]) * dt
			      + AVel(pi)[k];
            }
	}
    }

/*-----------------------------------------------------------------------------
 *  get_pos_vel  --  updates the positions and velocities in the newton0
 *                   representation, starting from the Aarseth-type
 *                   orbit extrapolation representation.
 *-----------------------------------------------------------------------------
 */
local void  get_pos_vel(bodies, extbodies, npart)
bodyptr  bodies;
orbitsegmentptr  extbodies;
int  npart;
    {
    bodyptr  body_i;      /* points to an individual particle of bodies[] */
    orbitsegmentptr  extbody_i;       /* points to an extbodies[] element */

    for (body_i = bodies, extbody_i = extbodies; body_i - bodies < npart;
	 body_i++, extbody_i++)
        {
	SETV(Pos(body_i), APos_now(extbody_i));
	SETV(Vel(body_i), AVel_now(extbody_i));
	}
    }

/* endof: transformext.c */
