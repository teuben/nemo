/* integrate.c - backward_euler, forward_euler, init_system_rev, leapfrog, 
                 runge_kutta_2, runge_kutta_4, runge_kutta_6, rk_fehlberg_45,
                 system_rev, system_step */

/*
 *  integrate.c: contains a choice of integration schemes for newton0
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
#include  "newton0.h"

static void forward_euler(systptr sys, specptr specs, real ds);
static void backward_euler(systptr sys, specptr specs, real ds);
static void runge_kutta_2(systptr sys, specptr specs, real ds);
static void runge_kutta_4(systptr sys, specptr specs, real ds);
static void runge_kutta_6(systptr sys, specptr specs, real ds);
static void rk_fehlberg_45(systptr sys, specptr specs, real ds);
static void leapfrog(systptr sys, specptr specs, real ds);


/*****************************************************************************/
/*                                                                           */
/*  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#  */
/*  @                                                                     @  */
/*  #                                                                     #  */
/*  @   integrate.c:                                                      @  */
/*  #                                                                     #  */
/*  @           a collection of integration schemes is offered here,      @  */
/*  #           in the form of recipes which use the lower-level          #  */
/*  @           procedures provided in  systemalgebra.c  and              @  */
/*  #           systemconversion.c  as basic ingredients.                 #  */
/*  @                                                                     @  */
/*  #           all integration schemes implemented in newton0            #  */
/*  @           reside in this file, and all are in the class of          @  */
/*  #           SELF-STARTING, COMMON TIME-STEP  methods, i.e.            #  */
/*  @           the time steps are equal for all particles at any         @  */
/*  #           given time, but the time length can vary through          #  */
/*  @           time.                                                     @  */
/*  #                                                                     #  */
/*  @   Part I: contains the administrative procedure  system_step()      @  */
/*  #           which selects the desired integration scheme.             #  */
/*  @                                                                     @  */
/*  #   Part II: integration schemes for FIRST-ORDER differential         #  */
/*  @            equations, which are applied to the system as a whole:   @  */
/*  #            the positions and velocities are considered as one       #  */
/*  @            2*NDIM*nbody - dimensional vector generalizations        @  */
/*  #            of the scalar prototypes, such as Runge-Kutta schemes.   #  */
/*  @                                                                     @  */
/*  #            A basic concept here is the differential of a system in  #  */
/*  @            the form of a psystem (containing pbodies with the full  @  */
/*  #            phase space information) in which the differentials      #  */
/*  @            of positions and velocities are treated simultaneously.  @  */
/*  #                                                                     #  */
/*  @   Part III: integration schemes for SECOND-ORDER differential       @  */
/*  #             equations, which treat the positions part and the       #  */
/*  @             velocities part of a system not necessarily in the      @  */
/*  #             same way. A prototype is the leap-frog scheme.          #  */
/*  @                                                                     @  */
/*  #             A basic concept is the use of separate differentials    #  */
/*  @             of positions and velocities in the form of csystems,    @  */
/*  #             which are selected and manipulated in different ways,   #  */
/*  @             and only later combined to form the full phase space    @  */
/*  #             information.                                            #  */
/*  @                                                                     @  */
/*  #                                                                     #  */
/*  @                                                                     @  */
/*  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#  */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#  */
/*  @                                                                     @  */
/*  #                                                                     #  */
/*  @                               PART I                                @  */
/*  #                                                                     #  */
/*  @           contains the administrative procedure  system_step()      @  */
/*  #           which selects the desired integration scheme.             #  */
/*  @                                                                     @  */
/*  #                                                                     #  */
/*  @                                                                     @  */
/*  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  system_step  --  performs one integration step for a whole nbody system,
 *                   with the choice of integration scheme determined by the
 *                   specs variable  Integrationscheme , and the choice of
 *                   softening method determined by the specs variable
 *                   Softfocus.
 *                   accepts: sys: pointer to a nbody system in standard
 *                                 newton0 form;
 *                          specs: pointer to the specifications block,
 *                                 which contains that part of the control and
 *                                 diagnostics information which is
 *                                 intertwined with the force calculation.
 *                             ds: the pseudo-time length of a system step;
 *                                 for non-regularized system, ds is simply
 *                                 the time differential dt, but for a
 *                                 regularized system ds is the time-like
 *                                 parameter introduced to regularize the
 *                                 equations of motion.
 *                   effect: when  system_step()  returns, the system to which
 *                           "sys" points has been evolved in (pseudo-) time
 *                           by an amount  ds.
 *              side effect: when  system_step()  returns, the potential and
 *                           acceleration have been updated to the current
 *                           time, through a call to  system_rev()  just
 *                           before returning, so that diagnostics can be 
 *                           obtained directly from the system.
 *                           note: this effort is not wasted: the remembered
 *                                 potential and acceleration which are stored
 *                                 in the system are available through the
 *                                 handles Pot() and Acc() and are used
 *                                 directly at the start of the execution of
 *                                 each  integration scheme -- by skipping a
 *                                 call to  deriv_system() .
 *              remark: deriv_system()  is encapsulated in the procedure
 *                      system_rev() , in order to keep an extra conceptual
 *                      barrier between the low-level procedures which directly
 *                      act on the systems and the higher-level entry procedure
 *                      system_step()  which deals only with intermediate-level
 *                      procedures such as  (*integration_scheme)()  and
 *                      system_rev() .
 *-----------------------------------------------------------------------------
 */

void  system_step(sys, specs, ds)
systptr  sys;
specptr  specs;
real  ds;
    {
    proc  integration_scheme;   /* procedure name, short for: void (*i_s)(); */
    bool  rev_needed;
    string  scheme;

/*
 * most schemes need to rev the engines at the end of an integration step,
 * to update the potential in order to allow an evaluation of the total
 * energy of the system (accelerations are produced as well as a useful
 * byproduct); those schemes that end up with updated potentials and
 * accelerations have no need for extra engine revving, and will set the
 * following flag to FALSE:
 */
    rev_needed = TRUE;
/*
 * select the desired integration scheme:
 */
    scheme = Integrationscheme(specs);

    if (streq("forward_euler", scheme))              /*  streq()  is a macro */
        integration_scheme = forward_euler;          /* defined in  stdinc.h */
    else if (streq("backward_euler", scheme))
        integration_scheme = backward_euler;
    else if (streq("runge_kutta_2", scheme))
        integration_scheme = runge_kutta_2;
    else if (streq("runge_kutta_4", scheme))
        integration_scheme = runge_kutta_4;
    else if (streq("runge_kutta_6", scheme))
        integration_scheme = runge_kutta_6;
    else if (streq("rk_fehlberg_45", scheme))
        integration_scheme = rk_fehlberg_45;
    else if (streq("leapfrog", scheme))
        {
        integration_scheme = leapfrog;
	rev_needed = FALSE;
#ifdef REGULARIZATION
	error("system_step: leapfrog not compatable with regularization\n");
#endif
        }
    else
	error("system_step: %s not implemented as an integration scheme\n",
              scheme);

    (*integration_scheme)(sys, specs, ds);

    if (rev_needed)
        system_rev(sys, specs);   /* updates the accelerations and potential */
    }

/*-----------------------------------------------------------------------------
 *  init_system_rev  --  revs the engines in the system "sys" without actually
 *                       advancing the system state; this only updates the
 *                       accelerations and potentials (or for a regularized 
 *                       system dQ/ds and dP/ds), so that the total energy
 *                       of the system can be computed.
 *                       note: cf. system_rev() below; the need for a separate
 *                             init_system_rev arises from the difference in 
 *                             treatment of regularized systems, for which 
 *                             the initial revving is 3D, and further ones 4D.
 *-----------------------------------------------------------------------------
 */
void  init_system_rev(sys, specs)
systptr  sys;
specptr  specs;
    {
    init_deriv_system(sys, specs);            /* resides in differentiate.c  */
    }    

/*-----------------------------------------------------------------------------
 *  system_rev  --  revs the engines in the system "sys" without actually
 *                  advancing the system state; this only updates the
 *                  accelerations and potentials (or for a regularized 
 *                  system dQ/ds and dP/ds), so that the total energy
 *                  of the system can be computed.
 *                  note: future versions may feature more than just one line;
 *                        but even now, system_rev() acts as a conceptual
 *                        barrier between the entry procedure  rev_engine()
 *                        in the file  orbit.c  and the procedures such as
 *                        deriv_system() in the file differentiate.c  which
 *                        directly act on components of the system.
 *-----------------------------------------------------------------------------
 */
void  system_rev(sys, specs)
systptr  sys;
specptr  specs;
    {
    deriv_system(sys, specs);               /* resides in differentiate.c    */
    }    

/*===========================================================================*/

/*****************************************************************************/
/*                                                                           */
/*  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#  */
/*  @                                                                     @  */
/*  #                                                                     #  */
/*  @                              PART II                                @  */
/*  #                                                                     #  */
/*  @            integration schemes for FIRST-ORDER differential         @  */
/*  #            equations, which are applied to the system as a whole:   #  */
/*  @            the positions and velocities are considered as one       @  */
/*  #            2*NDIM*nbody - dimensional vector generalizations        #  */
/*  @            of the scalar prototypes, such as Runge-Kutta schemes.   @  */
/*  #                                                                     #  */
/*  @            A basic concept here is the differential of a system     @  */
/*  #            state in the form of a psystem (containing the full      #  */
/*  @            phase space information) in which the differentials      @  */
/*  #            of positions and velocities are treated simultaneously.  #  */
/*  @                                                                     @  */
/*  #                                                                     #  */
/*  @                                                                     @  */
/*  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  forward_euler  --  integrates the orbits of all particles in a
 *                     gravitational many-body system over one timestep.
 *                     The method of integration is forward Euler.
 *-----------------------------------------------------------------------------
 */
local void  forward_euler(sys, specs, ds)
systptr  sys;
specptr  specs;
real  ds;
    {
    psystptr  k1;

    k1 = d_system(sys, specs, ds);   /* deriv_system(sys, npart) has already */
                                     /* been invoked by  system_rev() .      */
    push_system(sys, k1);

    rm_psystem(k1);
    }

/*-----------------------------------------------------------------------------
 *  backward_euler  --  integrates the orbits of all particles in a
 *                      gravitational many-body system over one timestep.
 *                      The method of integration is backward Euler.
 *-----------------------------------------------------------------------------
 */
local void  backward_euler(sys, specs, ds)
systptr  sys;
specptr  specs;
real  ds;
    {
    systptr  p_aux;
    psystptr  k1, k2;

    k1 = d_system(sys, specs, ds);   /* deriv_system(sys, npart) has already */
                                     /* been invoked by  system_rev() .      */
    p_aux = cp_system(sys);
    push_system(p_aux, k1);
    deriv_system(p_aux, specs);
    k2 = d_system(p_aux, specs, ds);
    rm_system(p_aux);

    push_system(sys, k2);

    rm_psystem(k1);
    rm_psystem(k2);
    }

/*-----------------------------------------------------------------------------
 *  runge_kutta_2  --  integrates the orbits of all particles in a
 *                     gravitational many-body system over one timestep.
 *                     The method of integration is a 2nd-order Runge-Kutta.
 *-----------------------------------------------------------------------------
 */
local void  runge_kutta_2(sys, specs, ds)
systptr  sys;
specptr  specs;
real  ds;
    {
    psystptr  k1, k2;     /* Runge-Kutta variables, see Abramowitz & Stegun, */
    psystptr  k_tot;      /*   Handbook of Mathematical Functions, sec. 25.5 */
    systptr  p_aux;

    k1 = d_system(sys, specs, ds);   /* deriv_system(sys, npart) has already */
                                     /* been invoked by  system_rev() .      */
    p_aux = cp_system(sys);
    push_system(p_aux, k1);
    deriv_system(p_aux, specs);
    k2 = d_system(p_aux, specs, ds);
    rm_system(p_aux);

    k_tot = add_psystem(k1, k2);
    scale_psystem(k_tot, 0.5);

    push_system(sys, k_tot);

    rm_psystem(k1);
    rm_psystem(k2);
    rm_psystem(k_tot);
    }

/*-----------------------------------------------------------------------------
 *  runge_kutta_4  --  integrates the orbits of all particles in a
 *                     gravitational many-body system over one timestep.
 *                     The method of integration is a 4th-order Runge-Kutta.
 *-----------------------------------------------------------------------------
 */
local void  runge_kutta_4(sys, specs, ds)
systptr  sys;
specptr  specs;
real  ds;
    {
    psystptr  k1, k2, k3, k4;     /* Runge-Kutta variables, see Abramowitz & */
    psystptr  k_tot;              /*   Stegun, Handbook of Mathematical      */
    psystptr  k_aux1, k_aux2;     /*   Functions, sec. 25.5                  */
    systptr  p_aux;

    k1 = d_system(sys, specs, ds);   /* deriv_system(sys, npart) has already */
                                     /* been invoked by  system_rev() .      */
    k_aux1 = mul_psystem(k1, 0.5);
    p_aux = cp_system(sys);
    push_system(p_aux, k_aux1);
    deriv_system(p_aux, specs);
    k2 = d_system(p_aux, specs, ds);
    rm_system(p_aux);
    rm_psystem(k_aux1);

    k_aux1 = mul_psystem(k2, 0.5);
    p_aux = cp_system(sys);
    push_system(p_aux, k_aux1);
    deriv_system(p_aux, specs);
    k3 = d_system(p_aux, specs, ds);
    rm_system(p_aux);
    rm_psystem(k_aux1);

    p_aux = cp_system(sys);
    push_system(p_aux, k3);
    deriv_system(p_aux, specs);
    k4 = d_system(p_aux, specs, ds);
    rm_system(p_aux);

    scale_psystem(k1, 1.0/6.0);
    scale_psystem(k2, 1.0/3.0);
    scale_psystem(k3, 1.0/3.0);
    scale_psystem(k4, 1.0/6.0);

    k_aux1 = add_psystem(k1, k2);
    k_aux2 = add_psystem(k3, k4);
    k_tot = add_psystem(k_aux1, k_aux2);
    rm_psystem(k_aux1);
    rm_psystem(k_aux2);
    push_system(sys, k_tot);

    rm_psystem(k1);
    rm_psystem(k2);
    rm_psystem(k3);
    rm_psystem(k4);
    rm_psystem(k_tot);
    }

/*-----------------------------------------------------------------------------
 *  runge_kutta_6  --  integrates the orbits of all particles in a
 *                     gravitational many-body system over one timestep.
 *                     The method of integration is a 6th-order Runge-Kutta,
 *                     with a coefficient set of type Butcher.
 *                     Litt.: M.K. Jain: Numerical Solution of Differential
 *                            Equations (1984, second edition, Wiley),
 *                            section (2.3.4).
 *-----------------------------------------------------------------------------
 */
local void  runge_kutta_6(sys, specs, ds)
systptr  sys;
specptr  specs;
real  ds;
    {
    psystptr  k1, k2, k3, k4, k5, k6, k7;
    psystptr  k_tot;
    psystptr  k_aux0, k_aux00, k_aux1, k_aux2, k_aux3, k_aux4, k_aux5, k_aux6;
    systptr  p_aux;

    k1 = d_system(sys, specs, ds);   /* deriv_system(sys, npart) has already */
                                     /* been invoked by  system_rev() .      */
                                                                /* k1 is set */
    k_aux1 = mul_psystem(k1, 1.0/3.0);
    p_aux = cp_system(sys);
    push_system(p_aux, k_aux1);
    rm_psystem(k_aux1);
    deriv_system(p_aux, specs);
    k2 = d_system(p_aux, specs, ds);                            /* k2 is set */
    rm_system(p_aux);

    k_aux2 = mul_psystem(k2, 2.0/3.0);     /* k_aux1 = 0 */
    p_aux = cp_system(sys);
    push_system(p_aux, k_aux2);
    rm_psystem(k_aux2);
    deriv_system(p_aux, specs);
    k3 = d_system(p_aux, specs, ds);                            /* k3 is set */
    rm_system(p_aux);

    k_aux1 = mul_psystem(k1, 1.0/12.0);
    k_aux2 = mul_psystem(k2, 1.0/3.0);
    k_aux3 = mul_psystem(k3, -1.0/12.0);
    k_aux00 = add_psystem(k_aux1, k_aux2);
    rm_psystem(k_aux1);
    rm_psystem(k_aux2);
    k_aux0 = add_psystem(k_aux00, k_aux3);
    rm_psystem(k_aux00);
    rm_psystem(k_aux3);
    p_aux = cp_system(sys);
    push_system(p_aux, k_aux0);
    rm_psystem(k_aux0);
    deriv_system(p_aux, specs);
    k4 = d_system(p_aux, specs, ds);                            /* k4 is set */
    rm_system(p_aux);

    k_aux1 = mul_psystem(k1, -1.0/16.0);
    k_aux2 = mul_psystem(k2, 9.0/8.0);
    k_aux3 = mul_psystem(k3, -3.0/16.0);
    k_aux4 = mul_psystem(k4, -3.0/8.0);
    k_aux0 = add_psystem(k_aux1, k_aux2);
    rm_psystem(k_aux1);
    rm_psystem(k_aux2);
    k_aux00 = add_psystem(k_aux0, k_aux3);
    rm_psystem(k_aux0);
    rm_psystem(k_aux3);
    k_aux0 = add_psystem(k_aux00, k_aux4);
    rm_psystem(k_aux00);
    rm_psystem(k_aux4);
    p_aux = cp_system(sys);
    push_system(p_aux, k_aux0);
    rm_psystem(k_aux0);
    deriv_system(p_aux, specs);
    k5 = d_system(p_aux, specs, ds);                            /* k5 is set */
    rm_system(p_aux);

    k_aux2 = mul_psystem(k2, 9.0/8.0);     /* k_aux1 = 0 */
    k_aux3 = mul_psystem(k3, -3.0/8.0);
    k_aux4 = mul_psystem(k4, -3.0/4.0);
    k_aux5 = mul_psystem(k5, 1.0/2.0);
    k_aux0 = add_psystem(k_aux2, k_aux3);
    rm_psystem(k_aux2);
    rm_psystem(k_aux3);
    k_aux00 = add_psystem(k_aux0, k_aux4);
    rm_psystem(k_aux0);
    rm_psystem(k_aux4);
    k_aux0 = add_psystem(k_aux00, k_aux5);
    rm_psystem(k_aux00);
    rm_psystem(k_aux5);
    p_aux = cp_system(sys);
    push_system(p_aux, k_aux0);
    rm_psystem(k_aux0);
    deriv_system(p_aux, specs);
    k6 = d_system(p_aux, specs, ds);                            /* k6 is set */
    rm_system(p_aux);

    k_aux1 = mul_psystem(k1, 9.0/44.0);
    k_aux2 = mul_psystem(k2, -9.0/11.0);
    k_aux3 = mul_psystem(k3, 63.0/44.0);
    k_aux4 = mul_psystem(k4, 18.0/11.0);
    k_aux6 = mul_psystem(k6,-16.0/11.0);    /* k_aux5 = 0 */
    k_aux00 = add_psystem(k_aux1, k_aux2);
    rm_psystem(k_aux1);
    rm_psystem(k_aux2);
    k_aux0 = add_psystem(k_aux00, k_aux3);
    rm_psystem(k_aux00);
    rm_psystem(k_aux3);
    k_aux00 = add_psystem(k_aux0, k_aux4);
    rm_psystem(k_aux0);
    rm_psystem(k_aux4);
    k_aux0 = add_psystem(k_aux00, k_aux6);
    rm_psystem(k_aux00);
    rm_psystem(k_aux6);
    p_aux = cp_system(sys);
    push_system(p_aux, k_aux0);
    rm_psystem(k_aux0);
    deriv_system(p_aux, specs);
    k7 = d_system(p_aux, specs, ds);                            /* k7 is set */
    rm_system(p_aux);

    scale_psystem(k1, 11.0/120.0);
    scale_psystem(k3, 27.0/40.0);        /* yes, there is no k2 term */
    scale_psystem(k4, 27.0/40.0);
    scale_psystem(k5, -4.0/15.0);
    scale_psystem(k6, -4.0/15.0);
    scale_psystem(k7, 11.0/120.0);

    k_aux00 = add_psystem(k1, k3);
    rm_psystem(k1);
    rm_psystem(k2);      /* possible earlier, but more systematic to do here */
    rm_psystem(k3);
    k_aux0 = add_psystem(k_aux00, k4);
    rm_psystem(k_aux00);
    rm_psystem(k4);
    k_aux00 = add_psystem(k_aux0, k5);
    rm_psystem(k_aux0);
    rm_psystem(k5);
    k_aux0 = add_psystem(k_aux00, k6);
    rm_psystem(k_aux00);
    rm_psystem(k6);
    k_tot = add_psystem(k_aux0, k7);                 /* k_tot is set */
    rm_psystem(k_aux0);
    rm_psystem(k7);
    push_system(sys, k_tot);              /* integration step is done */
    rm_psystem(k_tot);
    }

/*-----------------------------------------------------------------------------
 *  rk_fehlberg_45  --  integrates the orbits of all particles in a
 *                      gravitational many-body system over one timestep.
 *                      The method of integration is Runge-Kutta-Fehlberg,
 *                      with 5th order accuracy; it also computes an
 *                      intermediate 4th order result (the reason for the "45"
 *                      in the name). Comparing the two results gives a good
 *                      estimate for the accuracy obtained; however, in the 
 *                      present implementation I do not yet use this option.
 *                      Litt.: J.R. Rice: Numerical Methods, Software, and
 *                             analysis (1983, McGraw-Hill), eq. 9.2.16.
 *                             Note: in that equation, f(x,y)=f(y) since there
 *                                   is no explicit time dependence in Newton's
 *                                   law of gravity if the system is modeled by
 *                                   a single point in a space with
 *                                   2*NDIM*npart dimensions.
 *-----------------------------------------------------------------------------
 */
local void  rk_fehlberg_45(sys, specs, ds)
systptr  sys;
specptr  specs;
real  ds;
    {
    psystptr  k1, k2, k3, k4, k5, k6;
    psystptr  k_tot;
    psystptr  k_aux0, k_aux00, k_aux1, k_aux2, k_aux3, k_aux4, k_aux5;
    systptr  p_aux;

    k1 = d_system(sys, specs, ds);   /* deriv_system(sys, npart) has already */
                                     /* been invoked by  system_rev() .      */
                                                                /* k1 is set */
    k_aux1 = mul_psystem(k1, 0.25);
    p_aux = cp_system(sys);
    push_system(p_aux, k_aux1);
    rm_psystem(k_aux1);
    deriv_system(p_aux, specs);
    k2 = d_system(p_aux, specs, ds);                            /* k2 is set */
    rm_system(p_aux);

    k_aux1 = mul_psystem(k1, 3.0/32.0);
    k_aux2 = mul_psystem(k2, 9.0/32.0);
    k_aux0 = add_psystem(k_aux1, k_aux2);
    rm_psystem(k_aux1);
    rm_psystem(k_aux2);
    p_aux = cp_system(sys);
    push_system(p_aux, k_aux0);
    rm_psystem(k_aux0);
    deriv_system(p_aux, specs);
    k3 = d_system(p_aux, specs, ds);                            /* k3 is set */
    rm_system(p_aux);

    k_aux1 = mul_psystem(k1, 1932.0/2197.0);
    k_aux2 = mul_psystem(k2, -7200.0/2197.0);
    k_aux3 = mul_psystem(k3, 7296.0/2197.0);
    k_aux00 = add_psystem(k_aux1, k_aux2);
    rm_psystem(k_aux1);
    rm_psystem(k_aux2);
    k_aux0 = add_psystem(k_aux00, k_aux3);
    rm_psystem(k_aux00);
    rm_psystem(k_aux3);
    p_aux = cp_system(sys);
    push_system(p_aux, k_aux0);
    rm_psystem(k_aux0);
    deriv_system(p_aux, specs);
    k4 = d_system(p_aux, specs, ds);                            /* k4 is set */
    rm_system(p_aux);

    k_aux1 = mul_psystem(k1, 439.0/216.0);
    k_aux2 = mul_psystem(k2, -8.0);
    k_aux3 = mul_psystem(k3, 3680.0/513.0);
    k_aux4 = mul_psystem(k4, -845.0/4104.0);
    k_aux0 = add_psystem(k_aux1, k_aux2);
    rm_psystem(k_aux1);
    rm_psystem(k_aux2);
    k_aux00 = add_psystem(k_aux0, k_aux3);
    rm_psystem(k_aux0);
    rm_psystem(k_aux3);
    k_aux0 = add_psystem(k_aux00, k_aux4);
    rm_psystem(k_aux00);
    rm_psystem(k_aux4);
    p_aux = cp_system(sys);
    push_system(p_aux, k_aux0);
    rm_psystem(k_aux0);
    deriv_system(p_aux, specs);
    k5 = d_system(p_aux, specs, ds);                            /* k5 is set */
    rm_system(p_aux);

    k_aux1 = mul_psystem(k1, -8.0/27.0);
    k_aux2 = mul_psystem(k2, 2.0);
    k_aux3 = mul_psystem(k3, -3544.0/2565.0);
    k_aux4 = mul_psystem(k4, 1859.0/4104.0);
    k_aux5 = mul_psystem(k5, -11.0/40.0);
    k_aux00 = add_psystem(k_aux1, k_aux2);
    rm_psystem(k_aux1);
    rm_psystem(k_aux2);
    k_aux0 = add_psystem(k_aux00, k_aux3);
    rm_psystem(k_aux00);
    rm_psystem(k_aux3);
    k_aux00 = add_psystem(k_aux0, k_aux4);
    rm_psystem(k_aux0);
    rm_psystem(k_aux4);
    k_aux0 = add_psystem(k_aux00, k_aux5);
    rm_psystem(k_aux00);
    rm_psystem(k_aux5);
    p_aux = cp_system(sys);
    push_system(p_aux, k_aux0);
    rm_psystem(k_aux0);
    deriv_system(p_aux, specs);
    k6 = d_system(p_aux, specs, ds);                            /* k6 is set */
    rm_system(p_aux);

    scale_psystem(k1, 16.0/135.0);
    scale_psystem(k3, 6656.0/12825.0);   /* yes, there is no k2 term */
    scale_psystem(k4, 28561.0/56430.0);
    scale_psystem(k5, -9.0/50.0);
    scale_psystem(k6, 2.0/55.0);

    k_aux0 = add_psystem(k1, k3);
    rm_psystem(k1);
    rm_psystem(k2);      /* possible earlier, but more systematic to do here */
    rm_psystem(k3);
    k_aux00 = add_psystem(k_aux0, k4);
    rm_psystem(k_aux0);
    rm_psystem(k4);
    k_aux0 = add_psystem(k_aux00, k5);
    rm_psystem(k_aux00);
    rm_psystem(k5);
    k_tot = add_psystem(k_aux0, k6);                 /* k_tot is set */
    rm_psystem(k_aux0);
    rm_psystem(k6);
    push_system(sys, k_tot);              /* integration step is done */
    rm_psystem(k_tot);
    }

/*===========================================================================*/

/*****************************************************************************/
/*                                                                           */
/*  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#  */
/*  @                                                                     @  */
/*  #                                                                     #  */
/*  @                              PART III                               @  */
/*  #                                                                     #  */
/*  @             integration scheme for SECOND-ORDER differential        @  */
/*  #             equations, which treat the positions part and the       #  */
/*  @             velocities part of a system not necessarily in the      @  */
/*  #             same way. A prototype is the leap-frog scheme.          #  */
/*  @                                                                     @  */
/*  #             A basic concept is the use of separate differentials    #  */
/*  @             of positions and velocities in the form of csystems,    @  */
/*  #             which are selected and manipulated in different ways,   #  */
/*  @             and only later combined to form the full phase space    @  */
/*  #             information.                                            #  */
/*  @                                                                     @  */
/*  #                                                                     #  */
/*  @                                                                     @  */
/*  #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  leapfrog  --  integrates the orbits of all particles in a
 *                gravitational many-body system over one timestep.
 *                The method of integration is a leapfrog scheme.
 *                note: the present leapfrog variant is selfstarting, and can
 *                      be viewed as a mapping of system state at  t  to
 *                      system state at  t + ds . This is accomplished by
 *                      doubling the number of (implicit) velocity evaluations,
 *                      half of which are (explicitly) in sinc with the
 *                      positions and accelerations, while the other half are
 *                      (implicitly) at the half-interval-points in time
 *                      between the time at which pos. and acc. are evaluated.
 *-----------------------------------------------------------------------------
 */
local void  leapfrog(sys, specs, ds)
systptr  sys;
specptr  specs;
real  ds;
    {
    csystptr  v_old, a_old, a_new;
    csystptr  dr1, dr2, dr, dv, dv_helper;

    v_old = sel_vel(sys);                /* v_old                     */
    a_old = sel_acc(sys);                /* a_old                     */
    dr1 = mul_csystem(v_old, ds);
    rm_csystem(v_old);
    dr2 = mul_csystem(a_old, 0.5*ds*ds);
    dr = add_csystem(dr1, dr2);          /* dr = v_old*dt             */
    rm_csystem(dr1);                     /*      + a_old*dt*dt/2      */
    rm_csystem(dr2);
    annex_pos(sys, dr);                  /* r_new = r_old + dr        */
    rm_csystem(dr);
    deriv_system(sys, specs);            /* force calculation: a_new  */
    a_new = sel_acc(sys);
    dv_helper = add_csystem(a_old, a_new);
    rm_csystem(a_old);
    rm_csystem(a_new);
    dv = mul_csystem(dv_helper, 0.5*ds); /* dv = (a_old + a_new)*dt/2 */
    rm_csystem(dv_helper);
    annex_vel(sys, dv);                  /* v_new = v_old + dr        */
    rm_csystem(dr);
    rm_csystem(dv);
    Tnow(sys) += ds;                     /* t_new = t_old + dt        */
    }

/* endof: integrate.c */
