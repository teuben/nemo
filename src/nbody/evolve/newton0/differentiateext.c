/* differentiateext.c - base_update, get_force, integrate, 
	                next_body, poly_update, pred_low_order,
	                pred_high_order, timestep_update, trak_update */
/*
 *  differentiateext.c: contains procedures for nbody integration using
 *	                polynomials, based on Aarseth's NBODY0.
 *
 *      Feb 1988  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 *
 *  NOTE: this is a PRELIMINARY version: I have not yet included softening,
 *        or proper diagnostics such as closest encounter distances, etc.
 *        Also, the code could be made quite a bit more efficient. In due
 *        time I plan to improve these points.
 */

#include  "newton0.h"

local orbitsegmentptr  abody;     /* pointer to the first particle           */
                                  /*   (the name abody instead of body is    */
                                  /*    used here to avoid name conflict     */
                                  /*    with "body" as defined in standard   */
                                  /*    newton0 terms (see file  body.h ).   */
local int  nobj;                  /* number of particles                     */
local real  t_now;                /* present time in the simulation          */
local real  eta;           /* dimensionless integration accuracy parameter,  */
                           /* used in determining the size of the time steps */

static orbitsegmentptr next_body(void);
static void pred_low_order(orbitsegmentptr pj);
static void pred_high_order(orbitsegmentptr pi);
static void get_force(orbitsegmentptr pi, real force[], realptr min_pair_ptr);
static void poly_update(orbitsegmentptr pi, timesegments *dt, real force[], real d4_now[]);
static void base_update(orbitsegmentptr pi, timesegments *dt, real d4_now[]);
static void trak_update(orbitsegmentptr pi, timesegments *dt, real force[]);
static void timestep_update(orbitsegmentptr pi, timesegments *dt);


/*-----------------------------------------------------------------------------
 *  integrate_ext_invidual_timestep  --  the main integrating loop.
 *                                       before invoking this procedure, make
 *                                       sure to have converted the standard 
 *                                       newton0 representation into an
 *                                       Aarseth representation, using the
 *                                       procedures provided in the file
 *                                        transformext.c .
 *-----------------------------------------------------------------------------
 */
void  integrate_ext_invidual_timestep(the_state, final_time)
stateptr  the_state;
real  final_time;
    {
    int  j;
    orbitsegmentptr  pi, pj;	/* pointers to the i'th and j'th particle    */
    timesegments  time_int;	/* a structure containing all time intervals */
				/* used in the force polynomial; see         */
				/* "nbody0.h" for the detailed information   */
    real  force[NDIM];		/* total force on particle i                 */
				/* from all other particles j                */
    real  d4_now[NDIM];	   /* d4 coefficients in the force polynomial        */
			   /* needed for the semi-iteration in base_update() */
    real  t_integration;   /* total integration time before interruption     */
    realptr  min_pair_ptr;
    systptr  sys;
    specptr  specs;
    diagptr  diags;

    sys = System(the_state);
    specs = Specs(the_state);
    diags = Diags(the_state);
/*
 * initialize the global variables in this file
 *   (in a future version these variables should be passed explicitly)
 */
    abody = Orbitsegments(sys);
    nobj = Nbody(sys);
    t_now = Tnow(sys);
    eta = Stepparam(specs);
    min_pair_ptr = &SPECmin_pair(specs);
/*
 * keep pushing the first particle for which the extrpolation prediction
 * validity runs out first, as detected by  next_body() :
 */
    while (t_now < final_time)
        {
    	pi = next_body();
    	for (pj = abody, j = 0; j < nobj; pj++, j++)
            pred_low_order(pj);
        pred_high_order(pi);
        get_force(pi, force, min_pair_ptr);
	poly_update(pi, &time_int, force, d4_now);
	base_update(pi, &time_int, d4_now);
	trak_update(pi, &time_int, force);
	timestep_update(pi, &time_int);

	DIAGnsteps(diags) += 1;          /* counts number of steps performed */
        } 
    }

/*-----------------------------------------------------------------------------
 *  next_body --  return a pointer to the particle that ought to be moved next
 *-----------------------------------------------------------------------------
 */
local orbitsegmentptr  next_body()
    {
    int  j;
    orbitsegmentptr  pi, pj;

    pi = &abody[0];
    t_now = T0(pi) + ATimestep(pi);

    for (pj = &abody[1], j = 1; j < nobj; pj++, j++)
        if (t_now > T0(pj) + ATimestep(pj))
            {
            pi = pj;
            t_now = T0(pj) + ATimestep(pj);
            }
    return(pi);
    }

/*-----------------------------------------------------------------------------
 *  pred_low_order  --  predict the position of a particle to low order
 *		       (up to the jerk, the first time derivative of the force)
 *-----------------------------------------------------------------------------
 */
local void  pred_low_order(pj)
orbitsegmentptr  pj;	/* jth particle */
    {
    int  k;
    real  dt;		/*  advance jth particle over a time increase dt
			 *  using eq.(8), Aarseth (1985)
			 */
    real  dt2;		/*  dt/2  */
    real  dt3;		/*  dt/3  */
    
    dt = t_now - T0(pj);
    dt2 = dt/2;
    dt3 = dt/3;
    for (k = 0; k < NDIM; k++)
	APos_now(pj)[k] = ((AJerk(pj)[k]*dt3+AForce(pj)[k])*dt2+AVel(pj)[k])*dt
			  + APos(pj)[k];
    }
    

/*-----------------------------------------------------------------------------
 *  pred_high_order  --  predict the position and velocity of a particle
 *                       to higher order accuracy
 *			 (up to the third time derivative of the force)
 *-----------------------------------------------------------------------------
 */
local void  pred_high_order(pi)
orbitsegmentptr  pi;
    {
    int  k;
    real  f2i[NDIM];
    real  f3i[NDIM];
    real  dt;

    dt = t_now - T0(pi);

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
        APos_now(pi)[k] += ((f3i[k]/120) * dt + f2i[k]/24)*dt*dt*dt*dt;
        AVel(pi)[k] += ((((f3i[k]/24) * dt + f2i[k]/6) * dt +
    			AJerk(pi)[k]/2) * dt + AForce(pi)[k]) * dt;
        }
    }

#define  BIGNUMBER  1.0e20

/*-----------------------------------------------------------------------------
 *  get_force  --  compute the force on particle i exerted by all other
 *		   particles j and return it in the vector 'force'
 *-----------------------------------------------------------------------------
 */
local void  get_force(pi, force, min_pair_ptr)
orbitsegmentptr  pi;
real  force[];
realptr  min_pair_ptr;
    {
    int  j, k;
    orbitsegmentptr  pj;	/* pointer to the j'th particle */
    real  rij[NDIM];		/* vector from j to i */
    real  r2;			/* distance between i and j squared */
    real  mrinv3;		/* mass divided by |rij| cubed */
    real  r2_min;               /* minimum pair separation squared */

    for (k = 0; k < NDIM; k++)
        force[k] = 0;

    r2_min = BIGNUMBER;

    for (pj = abody, j = 0; j < nobj; pj++, j++)
        if (pj != pi)
            {
	    for (k = 0; k < NDIM; k++)
		rij[k] = APos_now(pi)[k] - APos_now(pj)[k];
            r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
            mrinv3 = AMass(pj)/(r2 * sqrt(r2));
	    for (k = 0; k < NDIM; k++)
		force[k] -= mrinv3 * rij[k];

	    r2_min = MIN(r2_min, r2);
            }

    *min_pair_ptr = MIN(*min_pair_ptr, sqrt(r2_min));
    }

/*-----------------------------------------------------------------------------
 *  poly_update  --  update the polynomial representation part of particle i
 *-----------------------------------------------------------------------------
 */
local void  poly_update(pi, dt, force, d4_now)
orbitsegmentptr  pi;
timesegments  *dt;
real  force[];
real  d4_now[];		/*  d4 coefficients in the force polynomial,       */
    {			/*  needed for the semi-iteration in base_update() */
    int  k;
    real  d1_now[NDIM], d2_now[NDIM], d3_now[NDIM];  /* updates D1(pi), etc. */

    ANow(dt)[0] = t_now - T0(pi);
    ANow(dt)[1] = t_now - T1(pi);
    ANow(dt)[2] = t_now - T2(pi);
    ANow(dt)[3] = t_now - T3(pi);
    AOld(dt)[0] = T0(pi) - T0(pi);	/* for symmetry's sake; not used */
    AOld(dt)[1] = T0(pi) - T1(pi);
    AOld(dt)[2] = T0(pi) - T2(pi);
    AOld(dt)[3] = T0(pi) - T3(pi);
    T3(pi) = T2(pi);
    T2(pi) = T1(pi);
    T1(pi) = T0(pi);
    T0(pi) = t_now;

    for (k = 0; k < NDIM; k++)
        {
	d1_now[k] = (force[k] - AForce(pi)[k])/ANow(dt)[0];
        d2_now[k] = (d1_now[k] - D1(pi)[k])/ANow(dt)[1];
        d3_now[k] = (d2_now[k] - D2(pi)[k])/ANow(dt)[2];
        d4_now[k] = (d3_now[k] - D3(pi)[k])/ANow(dt)[3];
        D1(pi)[k] = d1_now[k];
        D2(pi)[k] = d2_now[k];
        D3(pi)[k] = d3_now[k];
	}
    }

/*-----------------------------------------------------------------------------
 *  base_update  --  update the position and velocity of particle i
 *		     using semi-iteration to obtain even higher accuracy
 *		     (up to the fourth time derivative of the force)
 *-----------------------------------------------------------------------------
 */
local void  base_update(pi, dt, d4_now)
orbitsegmentptr  pi;
timesegments  *dt;
real  d4_now[];		/*  d4 coefficients in the force polynomial        */
    {			/*  needed for the semi-iteration in base_update() */
    int  k;
    real  df1i[NDIM];
    real  mdf2i[NDIM], mdf3i[NDIM], mdf4i[NDIM];
	/* the notation is a twofold extension of that used in the functions
	 * pred_high_order() and thirdpass() above where, e.g., f2i denoted the
	 * second time derivative of the total force on the i'th particle:
	 *     df2i  (differential of f2i)
	 * denotes the change added in the semi-iteration to the
	 * value of f2i calculated already before in pred_high_order();
	 *    mdf2i  (modified differential of f2i)
	 * denotes the incorporation of the factorial coefficients as
	 * given in Aarseth (1985), eq. (4), and therefore
	 *    mdf2i = df2i/2;    mdf3i = df3i/6;    mdf4i = df4i/24 .
	 * this is all done to avoid unnecessary multiplications leading
	 * to unnecessary loss of accuracy (although at the expense of clarity)
	 */
    real  tstep;		/* most recent time step for i'th particle */

    tstep = ANow(dt)[0];

    for (k = 0; k < NDIM; k++)
        {
        APos(pi)[k] = APos_now(pi)[k];

        df1i[k]  = AOld(dt)[1] * AOld(dt)[2] * AOld(dt)[3] * d4_now[k];
        mdf2i[k] = (AOld(dt)[1] * AOld(dt)[2]
		   + AOld(dt)[3] * (AOld(dt)[1] + AOld(dt)[2])) * d4_now[k];
        mdf3i[k] = (AOld(dt)[1] + AOld(dt)[2] + AOld(dt)[3]) * d4_now[k];
	mdf4i[k] = d4_now[k];

        APos(pi)[k] += (((mdf4i[k]*tstep/30 + 0.05*mdf3i[k])*tstep
			+ mdf2i[k]/12)*tstep + df1i[k]/6)*tstep*tstep*tstep;
        AVel(pi)[k] += (((0.2*mdf4i[k]*tstep + 0.25*mdf3i[k])*tstep
			+ mdf2i[k]/3)*tstep + 0.5*df1i[k])*tstep*tstep;
        }
    }

/*-----------------------------------------------------------------------------
 *  trak_update  --  update the force and jerk
 *		     in the 'trak' representation of particle i
 *-----------------------------------------------------------------------------
 */
local void  trak_update(pi, dt, force)
orbitsegmentptr  pi;
timesegments  *dt;
real  force[];
    {
    int  k;

    for (k = 0; k < NDIM; k++)
        {
        AForce(pi)[k] = force[k];
        AJerk(pi)[k] = (D3(pi)[k] * ANow(dt)[1] + D2(pi)[k]) * ANow(dt)[0]
		       + D1(pi)[k];
        }
    }

/*-----------------------------------------------------------------------------
 *  timestep_update  -- update the timestep for particle i
 *-----------------------------------------------------------------------------
 */
local void  timestep_update(pi, dt)
orbitsegmentptr  pi;		  /* pointer to the i'th particle */
timesegments  *dt;
    {
    int  k;
    real  f2i[NDIM];      /* second time derivative of force on ith particle */
    real  f3i[NDIM];      /* third time derivative of force on ith particle  */
    real  fi_square;	  /* absolute value squared of force on ith particle */
    real  f1i_square;	  /* abs.v.sq. of 1st time derivative of force on i  */
    real  f2i_square;	  /* abs.v.sq. of 2nd time derivative of force on i  */
    real  f3i_square;	  /* abs.v.sq. of 3rd time derivative of force on i  */
    real  tstep_sq;	  /* square of time step scale                       */

    for (k = 0; k < NDIM; k++)
        {
        f2i[k] = 2 * (D3(pi)[k]*(ANow(dt)[0] + ANow(dt)[1]) + D2(pi)[k]);
	f3i[k] = 6 * D3(pi)[k];
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
    }

/* end of: differentiateext.c */
