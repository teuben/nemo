/* timestep.c - BIGNUMBER, constant_timestep, collision_time_timestep,
                find_timestep_method, proj_free_time */

/*
 *  timestep.c: different time step criteria for newton0
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
#include  "newton0.h"

static real proj_free_time(systptr sys, real r0_soft);


/*-----------------------------------------------------------------------------
 *  find_timestep_method()  --  returns the pointer  ptr  to the function which
 *                              determines the new common timestep, 
 *                              to be shared by all particles.
 *-----------------------------------------------------------------------------
 */
rproc  find_timestep_method(type_of_timestep)
string  type_of_timestep;
    {
    if (streq("constant", type_of_timestep))
	return( constant_timestep );
    else if (streq("collision_time", type_of_timestep))
	return( collision_time_timestep );
    else
	error("find_timestep_method: %s not implemented as timestep method",
                                                             type_of_timestep);
    }

/*-----------------------------------------------------------------------------
 *  constant_timestep  --  determines the new common timestep, to be shared by
 *                         all particles. The absolute value of the timestep
 *                         has the constant value  Eta_acc .
 *                         The sign of the new timestep is determined
 *                         depending on whether the integration is
 *                         forwards or backwards in time.
 *                         accepts: some_state: pointer to a nbody state in 
 *                                              standard newton0 form.
 *                         returns: the timestep.
 *-----------------------------------------------------------------------------
 */
real  constant_timestep(some_state)
stateptr  some_state;
    {
    real  dt;

    dt = Stepparam(Specs(some_state));
    if (Forwards(Ctrls(some_state)))
	return(MIN(dt, DTmax(Ctrls(some_state))));
    else                                 /* dt_max_step < 0 set in newton0.c */
	return(MAX(-dt, DTmax(Ctrls(some_state))));
    }

/*-----------------------------------------------------------------------------
 *  collision_time_timestep  --  determines the new common timestep, to be 
 *                               shared by all particles. The absolute value of
 *                               the timestep is a small fraction  Stepparam of
 *                               the time scale "proj_free_time()" on which
 *                               close encounters are predicted to occur by
 *                               linear extrapolation, as long as it does not
 *                               exceed  DTmax .
 *                               The sign of the new timestep is determined
 *                               depending on whether the integration is
 *                               forwards or backwards in time.
 *                               accepts: some_state: pointer to a nbody state
 *                                                    in standard newton0 form.
 *-----------------------------------------------------------------------------
 */
real  collision_time_timestep(some_state)
stateptr  some_state;
    {
    real  dt;

    dt = Stepparam(Specs(some_state))
            * proj_free_time(System(some_state), Softparam(Specs(some_state)));
    if (Forwards(Ctrls(some_state)))
	return(MIN(dt, DTmax(Ctrls(some_state))));
    else                                 /* dt_max_step < 0 set in newton0.c */
	return(MAX(-dt, DTmax(Ctrls(some_state))));
    }

#define  BIGNUMBER  1.0e20

/*-----------------------------------------------------------------------------
 *  proj_free_time  --  estimates the time at which a collision might occur
 *                      somewhere in the system. Expressed more accurately:
 *                      the relative distance and speed in each particle pair
 *                      is determined, and from these the time scale of change;
 *                      accepts: some_state: pointer to a nbody state in 
 *                                              standard newton0 form.
 *                      returns: min_coll_time: the minimum of all these
 *                                              pair time scales.
 *                      NOTE: when softening is used, an effective minimum
 *                            distance is used equal to the softening length.
 *-----------------------------------------------------------------------------
 */
local real  proj_free_time(sys, r0_soft)
systptr  sys;
real r0_soft;
    {
    int  npart;
    real  del_r_sqr;             /* square of effective position separ. vect.*/
    real  del_v_sqr;             /* square of velocity separation vector     */
    real  coll_time_sqr;         /* projected duration until a possible      */
                                 /* collision between a particle pair (i,j)  */
    real  min_coll_time_sqr;     /* minimum of all such collision times      */
    bodyptr  bodies;             /* point to the first particle              */
    bodyptr  body_i, body_j;     /* point to individual particles            */
    real  del_r_ji[NDIM];        /* separation vector from i to j            */
    real  del_v_ji[NDIM];        /* velocity separation vector from i to j   */
    real  r0_sqr;                /* square of softening parameter            */

    npart = Nbody(sys);
    bodies = Bodies(sys);

    r0_sqr = r0_soft * r0_soft;
    min_coll_time_sqr = BIGNUMBER;
    for (body_i = bodies + 1; body_i - bodies < npart; body_i++)
        for (body_j = bodies; body_j < body_i; body_j++)
	    {
	    SUBV(del_r_ji, Pos(body_j), Pos(body_i));
            DOTVP(del_r_sqr, del_r_ji, del_r_ji);
            del_r_sqr += r0_sqr;         /* from real to effective distance */
	    SUBV(del_v_ji, Vel(body_j), Vel(body_i));
            DOTVP(del_v_sqr, del_v_ji, del_v_ji);
            coll_time_sqr = del_r_sqr / del_v_sqr;
	    if (coll_time_sqr < min_coll_time_sqr)
	        min_coll_time_sqr = coll_time_sqr;
	    }
    return (sqrt(min_coll_time_sqr));       /* from square root to real time */
    }

/* endof: timestep.c */
