/* orbit.c - propagate, evolve, rev_engine, time_out, too_far */

/*
 *  orbit.c: orbit integration module for newton0
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
#include  "newton0.h"

static int propagate(stateptr some_state, real t_out);
static int time_out(stateptr some_state);
static bool too_far(stateptr some_state, real tnow);
/* REGULARIZATION */
static void goto_hyperspace(stateptr new_state);
/* EXTRAPOLATION */
static void goto_orbitsegments(stateptr new_state);
static int ext_propagate(stateptr some_state, real t_out);


extern real  cputime();


#define  HALT  -1         /* error condition returned by  propagate()        */
#define  OKAY   0         /* indicates succesful completion of  propagate()  */

/*-----------------------------------------------------------------------------
 *  rev_engine  --  invoke the integrator without performing an actual
 *                  integration step, with the purpose of updating the
 *                  the potential and the acceleration of the system (or
 *                  dQ/ds and dP/ds and total energy for the case of a 
 *                  regularized system).
 *-----------------------------------------------------------------------------
 */
void  rev_engine(new_state)
stateptr  new_state;
    {
    init_system_rev(System(new_state), Specs(new_state));
                                          /* resides in the file integrate.c */
#ifdef REGULARIZATION
    goto_hyperspace(new_state);
#endif
#ifdef EXTRAPOLATION               /* first statement could be skipped; but  */
    goto_orbitsegments(new_state); /* perhaps useful for initial diagnostics */
#endif
    }

/*-----------------------------------------------------------------------------
 *  evolve  --  orchestrates a full orbit integration.
 *              In the following endless loop, each cycle performs
 *              orbit integration until the next output time. There
 *              are three ways to break out of the endless loop,
 *              the first two most gracefully, the last one more abrupt:
 *              (1) the simulated time passes the desired ending time  Tend ;
 *              (2) a soft time limit is reached, for example when
 *                  the CPU time exceeds a soft time limit, given by
 *                   CPUlastcall , after which integration is allowed
 *                  to proceed until the next major output time;
 *                  this in indicated by the value HALT returned by
 *                  the procedure "time_out()" which oversees output.
 *              (3) a hard time limit is reached, by exceeding the
 *                  maximal number of integration steps  Nmaxstep ,
 *                  or by exceeding the hard CPU limit  CPUmax ,
 *                  or by any other condition which may have triggered
 *                  the procedure "propagate()" to return the value HALT,
 *                  indicating that some bound has been exceeded.
 *-----------------------------------------------------------------------------
 */

void  evolve(some_state)
stateptr  some_state;
    {
    int  halt_flag;
    bool  forwards;
    real  tnow;
    real  till_then;
    ctrlptr  ctr;

    ctr = Ctrls(some_state);
    forwards = Forwards(ctr);             /* determine the direction of time */

    for(;;)          /* forever, i.e. until a  break  command is encountered */
        {                                     /* one orbit integration cycle */
        till_then = earliest(Tminor(ctr), Tmajor(ctr));
	till_then = earliest(till_then, Tsave(ctr));
	till_then = earliest(till_then, Tend(ctr));

#ifndef EXTRAPOLATION
        halt_flag = propagate(some_state, till_then);
#else
        halt_flag = ext_propagate(some_state, till_then);
#endif
        tnow = Tnow(System(some_state));

        if (halt_flag == HALT)            /* if CPU or other limits exceeded */
	    break;
        if (later(tnow, till_then))     /* oops, slippery path thru time ... */
	    announce("\n evolve: overshooting propagate cycle by dt = %g\n",
                     ABS(tnow - till_then));

        if (! earlier(tnow, Tend(ctr)))   /* normal succesful completion     */
	    break;                        /* leave last output for finish()  */
        else
	    halt_flag = time_out(some_state);  /* take a break, and return   */
                                               /* after output has been done */
        if (halt_flag == HALT)
	    break;
        }
    }

/*-----------------------------------------------------------------------------
 *  propagate  --  pushes all particles along their orbits until time  t_out .
 *                 accepts: some_state: pointer to a state in standard newton0
 *                               form;
 *                        t_out: when  Tnow  reaches this time, control is
 *                               returned to the calling procedure.
 *                 returns: halt_flag: OKAY upon successful completion;
 *                                     HALT otherwise.
 *-----------------------------------------------------------------------------
 */
local int  propagate(some_state, t_out)
stateptr  some_state;
real  t_out;
    {
    int  halt_flag;
    bool  forwards;
    real  tnow;
    real  ds;                                        /* time parameter */
    systptr  sys;
    specptr  specs;
    ctrlptr  ctr;
    diagptr  diags;
    rproc  timestep_method;        /* determines the time step length.       */
    int  niter;                    /* number of iterations in the last       */
                                   /* time step before arriving at t_out.    */

#ifndef REGULARIZATION
    sys = System(some_state);
#else
    sys = Regsystem(some_state);
#endif
    specs = Specs(some_state);
    ctr = Ctrls(some_state);
    diags = Diags(some_state);

#ifndef REGULARIZATION
    niter = 0;
#else
    niter = Ntimingiter(ctr);
#endif
/*
 * determine the direction of time 
 */
    forwards = Forwards(ctr);
/*
 * select the desired procedure for determining the time step length:
 */
    timestep_method = find_timestep_method(Timestepmethod(specs));

/*
 * main integration loop:
 */
    halt_flag = OKAY;
    while ( (earlier(Tnow(sys), t_out) || niter > 0 )
	    && (DIAGnsteps(diags) < Nmajor(ctr))      )
	{
        tnow = Tnow(sys);

        if (! earlier(tnow, t_out))
	    niter--;

        if (too_far(some_state, tnow))
            {
	    halt_flag = HALT;
            break;
	    }

	ds = (*timestep_method)(some_state);

#ifndef REGULARIZATION
	if (later(tnow + ds, t_out))
	    ds = t_out - tnow;
#else                                            /* ds is d(fictitious time) */
	if (later(tnow + ds/Reglagrangian(sys), t_out))
	    ds = (t_out - tnow) * Reglagrangian(sys);
#endif

	system_step(sys, specs, ds);           /* advances positions,        */
                                               /* velocities & accelerations */
	DIAGnsteps(diags) += 1;          /* counts number of steps performed */
	DIAGde_nsteps(diags) += 1;       /* steps after last energy check    */
	}

#ifdef REGULARIZATION
    Tnow(System(some_state)) = Tnow(Regsystem(some_state));
#endif

    return(halt_flag);
    }

/*-----------------------------------------------------------------------------
 *  time_out  --  takes a break to oversee the output activities, after which
 *                control is returned to the  evolve()  procedure.
 *                returns: halt_flag: OKAY upon successful completion;
 *                                    HALT when evolution should be halted.
 *-----------------------------------------------------------------------------
 */
local int  time_out(some_state)
stateptr  some_state;
    {
    int  halt_flag;
    bool  forwards;
    real  tnow;
    ctrlptr  ctr;
    diagptr  diags;

    halt_flag = OKAY;
    tnow = Tnow(System(some_state));
    ctr = Ctrls(some_state);
    diags = Diags(some_state);
    forwards = Forwards(ctr);

    if ( (! earlier(tnow, Tmajor(ctr)))
	|| (DIAGnsteps(diags) >= Nmajor(ctr)) )
        {                                     /* attempt a major output,     */
	if (cputime() > CPUlastcall(ctr))     /* but if past "last call"     */
	    halt_flag = HALT;                 /* no new orders are accepted, */
	else
	    maj_out(some_state, FALSE);       /* maj_out() updates t_maj_out */
        }
    else if (! earlier(tnow, Tminor(ctr)))    /* minor only if no major      */
        min_out(some_state, FALSE);           /* min_out() updates t_min_out */
    if (! earlier(tnow, Tsave(ctr)))          /* after all output done       */
        save_state(some_state);               /* save_state() updates t_save */

    return(halt_flag);
    }

/*-----------------------------------------------------------------------------
 *  too_far  --  checks whether any hard limits has been exceeded which should
 *               cause integration to be halted.
 *-----------------------------------------------------------------------------
 */
local bool  too_far(some_state, tnow)
stateptr  some_state;
real  tnow;
    {
    systptr  sys;
    specptr  specs;
    ctrlptr  ctr;
    diagptr  diags;
    real  total_energy;
    real  delta_total_energy;    /* energy drift, compared to original value */
    bool  gone_too_far = FALSE;

#ifndef REGULARIZATION
    sys = System(some_state);
#else
    sys = Regsystem(some_state);
#endif
    specs = Specs(some_state);
    ctr = Ctrls(some_state);
    diags = Diags(some_state);

    if (DIAGnsteps(diags) >= Nmaxstep(ctr))
        {
	announce("\n\tmax. number of integration steps exceeded at time %g\n",
                 tnow);
        gone_too_far = TRUE;
	}

    if (DIAGde_nsteps(diags) >= Nstep_de(specs))
	{
#ifndef REGULARIZATION                      /* in file  diagnose.c  */
        total_energy = compute_total_energy_only(sys, diags);
#else
        total_energy = Reghamiltonian(Regsystem(some_state));
#endif
	delta_total_energy = total_energy - DIAGetot(diags)[INITIAL];

	if (ABS(delta_total_energy) > DEmax(ctr))
	    {
	    announce("\n\ttoo much energy drift detected at time %g\n",
                      tnow);
	    gone_too_far = TRUE;
	    }
	else
	    DIAGde_nsteps(diags) = 0;   /* start a new round */
	}

    if (SPECmin_pair(specs) < Min_pairdist(ctr))
        {
	announce("\n\tminimum pair distance passed at time %g\n", tnow);
        gone_too_far = TRUE;
	}

    if (cputime() >= CPUmax(ctr))
	{
	announce("\n\tCPU time limit exceeded at time %g\n", tnow);
        gone_too_far = TRUE;
	}

    return(gone_too_far);
    }

#ifdef REGULARIZATION

/*-----------------------------------------------------------------------------
 *  goto_hyperspace  --  i.e. go to the final frontier: the 4th dimension.
 *                       effect: creates a regularized system (making use of
 *                               the total energy, computed in the previous
 *                               revving of the nonregularized system),
 *                               and revs this regularized system (to set it
 *                               up for the first integration time step).
 *-----------------------------------------------------------------------------
 */
local void  goto_hyperspace(new_state)
stateptr  new_state;
    {
    systptr  sys, regsys;
    specptr  specs;
    diagptr  diags;

    sys = System(new_state);
    Regsystem(new_state) = setup_reg(sys);              /* in transformreg.c */
    regsys = Regsystem(new_state);
    specs = Specs(new_state);
    diags = Diags(new_state);
                              
    compute_energy(sys, diags);                     /* needed for Regenergy  */
    Regenergy(regsys) = Reghamiltonian(regsys) =    /* needed for system_rev */
                        DIAGekin(diags)[CURRENT] + DIAGepot(diags)[CURRENT];
    Reglagrangian(regsys) = DIAGekin(diags)[CURRENT]-DIAGepot(diags)[CURRENT];
    system_rev(regsys, specs);
    }

#endif

#ifdef EXTRAPOLATION

/*-----------------------------------------------------------------------------
 *  goto_orbitsegments  --  make a polynomial approximation to the orbits of
 *                          all particles, as in Aarseth-type codes.
 *-----------------------------------------------------------------------------
 */
local void  goto_orbitsegments(new_state)
stateptr  new_state;
    {
    systptr  sys;
    specptr  specs;

    sys = System(new_state);
    specs = Specs(new_state);
                              
    setup_ext(sys, specs);                             /* in transformext.c */
    }

/*-----------------------------------------------------------------------------
 *  ext_propagate  --  polynomial-orbit-approximation version of  propagate() .
 *                     accepts: some_state: pointer to a state in standard
 *                              newton0 form;
 *                              t_out: when  Tnow  reaches this time, control
 *                                     is returned to the calling procedure.
 *                     returns: halt_flag: OKAY upon successful completion;
 *                                         HALT otherwise.
 *                     NOTE: this is a working version, but should soon be
 *                           augmented with an implementation of  too_far()
 *                           and softening choices, time reversal, etc; also
 *                           diagnostics such as number of steps, close
 *                           encounter information, etc.
 *-----------------------------------------------------------------------------
 */
local int  ext_propagate(some_state, t_out)
stateptr  some_state;
real  t_out;
    {
    systptr  sys;
    specptr  specs;
    int  halt_flag;

    sys = System(some_state);
    specs = Specs(some_state);

    halt_flag = OKAY;         /* in the future, breaks should be implemented */
/*
 * perform an integration using polynomial orbit integration,
 *   (the procedure resides in file  differentiateext.c)
 */
    integrate_ext_invidual_timestep(some_state, t_out);

    Tnow(sys) = t_out;
/*
 * translate the results back into standard newton0 representation,
 * while keeping the polynomial orbit approximation as well:
 *   (the procedure resides in file  transformext.c)
 */
    read_ext(sys, specs, t_out);

    return(halt_flag);
    }

#endif

/* endof: orbit.c */
