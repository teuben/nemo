/* diagnose.c - check_energy, compute_energy, compute_total_energy_only,
                diagnostics, final_diagnostics, grinding_halt,
                micro_diagnostics */

/*
 *  diagnose.c: diagnostics module for newton0.c : equal time steps
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
#include  "newton0.h"

static void check_energy(stateptr the_state, int init_flag);
static bool isbinary(bodyptr body1, bodyptr body2);
static bool is_ancestor_or_self(int i, int j, int nbody, int membership_list[200 ][3]);
static void init_binary(bodyptr body0, bodyptr body1, bodyptr body2, realptr aptr, realptr eptr);
static void scatter_diagnostics(stateptr the_state, int init_flag);
static void micro_diagnostics(stateptr the_state, int init_flag);
static void print_object(int i, int nbody, int membership_list[200 ][3]);


/*-----------------------------------------------------------------------------
 *  diagnostics  --  performs a generic set of diagnostics
 *-----------------------------------------------------------------------------
 */
void  diagnostics(the_state, init_flag)
stateptr  the_state;
bool  init_flag;
    {
    if (init_flag)                       /* should perhaps go in  newton0.c  */
        DIAGnsteps(Diags(the_state)) = DIAGde_nsteps(Diags(the_state)) = 0;

    check_energy(the_state, init_flag);

    if (streq("microscopic", Diagnostics(Specs(the_state))))
	micro_diagnostics(the_state, init_flag);

    if (streq("scatter", Diagnostics(Specs(the_state))))
	scatter_diagnostics(the_state, init_flag);

    /* .... */                   /* more to follow in future implementations */
    }

/*-----------------------------------------------------------------------------
 *  final_diagnostics  --  to be performed at closing time
 *-----------------------------------------------------------------------------
 */
void  final_diagnostics(the_state)
stateptr  the_state;
    {
    ;
    /* .... */                   /* more to follow in future implementations */
    }

/*-----------------------------------------------------------------------------
 *  check_energy  --  compute global kinetic, potential and total energy,
 *                    and compare with previous values, since nowadays it has
 *                    become progressive to conserve energy.
 *-----------------------------------------------------------------------------
 */
local void  check_energy(the_state, init_flag)
stateptr  the_state;
bool  init_flag;
    {
    real  kinetic_energy;
    real  velocity_squared;
    bodyptr  body_i;
    diagptr  diags;
    systptr  sys;
#ifdef REGULARIZATION
    systptr  regsys;
#endif

    sys = System(the_state);
    diags = Diags(the_state);

    if (! init_flag)
	{
        DIAGtime(diags)[PREVIOUS] = DIAGtime(diags)[CURRENT];
        DIAGetot(diags)[PREVIOUS] = DIAGetot(diags)[CURRENT];
        DIAGekin(diags)[PREVIOUS] = DIAGekin(diags)[CURRENT];
        DIAGepot(diags)[PREVIOUS] = DIAGepot(diags)[CURRENT];
        }

#ifndef REGULARIZATION
    compute_energy(sys, diags);
#else
    regsys = Regsystem(the_state);
    DIAGekin(diags)[CURRENT] = (Reghamiltonian(Regsystem(the_state))
                                + Reglagrangian(Regsystem(the_state))) / 2.0;
    DIAGepot(diags)[CURRENT] = (Reghamiltonian(Regsystem(the_state))
                                - Reglagrangian(Regsystem(the_state))) / 2.0;
#endif

    if (init_flag)
	{
        DIAGtime(diags)[INITIAL] = DIAGtime(diags)[PREVIOUS] =
                                          DIAGtime(diags)[CURRENT] = Tnow(sys);
        DIAGetot(diags)[INITIAL] = DIAGetot(diags)[PREVIOUS] =
                           DIAGetot(diags)[CURRENT] =
                           DIAGekin(diags)[CURRENT] + DIAGepot(diags)[CURRENT];
        }
    else
	{
	DIAGtime(diags)[CURRENT] = Tnow(sys);
	DIAGetot(diags)[CURRENT] =
                           DIAGekin(diags)[CURRENT] + DIAGepot(diags)[CURRENT];
	}
    }

/*-----------------------------------------------------------------------------
 *  compute_energy  --  compute global kinetic and potential energy
 *                      note: only the relative kinetic energy is computed,
 *                            with respect to the center of mass of the system.
 *-----------------------------------------------------------------------------
 */
void  compute_energy(sys, diags)
systptr  sys;
diagptr  diags;
    {
    int  nbody;
    real  kinetic_energy;
    real  velocity_squared;
    real  sum_kin_energy;
    real  sum_pot_energy;
    real  com_velocity[NDIM];
    real  mi_vi[NDIM];
    real  total_mass;
    real  com_kin_energy;
    real  com_vel_squared;
    bodyptr  bodies;
    bodyptr  body_i;

    sum_kin_energy = sum_pot_energy = 0.0;
    bodies = Bodies(sys);
    nbody = Nbody(sys);
    for (body_i = bodies; body_i - bodies < nbody; body_i++)
	{
        DOTVP(velocity_squared, Vel(body_i), Vel(body_i));
	kinetic_energy = Mass(body_i) * velocity_squared;    /* 2 x too high */
	sum_kin_energy += kinetic_energy;
	sum_pot_energy += Mass(body_i) * Pot(body_i);  /* counts pairs twice */
	}
    DIAGekin(diags)[CURRENT] = 0.5 * sum_kin_energy;
                                       /* corrects for the factor 2 left out */
    DIAGepot(diags)[CURRENT] = 0.5 * sum_pot_energy;
                                       /* corrects for double counting pairs */

    total_mass = 0.0;
    CLRV(com_velocity);
    for (body_i = bodies; body_i - bodies < nbody; body_i++)
	{
	total_mass += Mass(body_i);
	MULVS(mi_vi, Vel(body_i), Mass(body_i));
	INCADDV(com_velocity, mi_vi);                 /* here: mass-weighted */
	}
    INCDIVVS(com_velocity, total_mass);     /* here: the real c.o.m.velocity */
    DOTVP(com_vel_squared, com_velocity, com_velocity);
    com_kin_energy = 0.5 * total_mass * com_vel_squared;
    DIAGepot(diags)[CURRENT] -= com_kin_energy;
    }

/*-----------------------------------------------------------------------------
 *  compute_total_energy_only  --  compute total energy for checking purposes
 *                                 note: only the relative kinetic energy
 *                                       is computed, with respect to the 
 *                                       center of mass of the system.
 *-----------------------------------------------------------------------------
 */
real  compute_total_energy_only(sys, diags)
systptr  sys;
diagptr  diags;
    {
    int  nbody;
    real  kinetic_energy;
    real  velocity_squared;
    real  sum_kin_energy;
    real  sum_pot_energy;
    real  com_velocity[NDIM];
    real  mi_vi[NDIM];
    real  total_mass;
    real  com_kin_energy;
    real  com_vel_squared;
    bodyptr  bodies;
    bodyptr  body_i;

    sum_kin_energy = sum_pot_energy = 0.0;
    bodies = Bodies(sys);
    nbody = Nbody(sys);
    for (body_i = bodies; body_i - bodies < nbody; body_i++)
	{
        DOTVP(velocity_squared, Vel(body_i), Vel(body_i));
	kinetic_energy = Mass(body_i) * velocity_squared;    /* 2 x too high */
	sum_kin_energy += kinetic_energy;
	sum_pot_energy += Mass(body_i) * Pot(body_i);  /* counts pairs twice */
	}
    sum_kin_energy *= 0.5;             /* corrects for the factor 2 left out */
    sum_pot_energy *= 0.5;             /* corrects for double counting pairs */

    total_mass = 0.0;
    CLRV(com_velocity);
    for (body_i = bodies; body_i - bodies < nbody; body_i++)
	{
	total_mass += Mass(body_i);
	MULVS(mi_vi, Vel(body_i), Mass(body_i));
	INCADDV(com_velocity, mi_vi);                 /* here: mass-weighted */
	}
    INCDIVVS(com_velocity, total_mass);     /* here: the real c.o.m.velocity */
    DOTVP(com_vel_squared, com_velocity, com_velocity);
    com_kin_energy = 0.5 * total_mass * com_vel_squared;

    return( sum_pot_energy + sum_kin_energy - com_kin_energy );
    }

/*-----------------------------------------------------------------------------
 *  grinding_halt  --  halt the present integration, supply information to the
 *                     proper output channels about the error condition which
 *                     prompted this premature halt, and exit as graceful as
 *                     possible under the circumstances.
 *-----------------------------------------------------------------------------
 */
void  grinding_halt(the_state)
stateptr  the_state;
    {
    announce("\ngrinding_halt: not yet implemented\n",0.0);
    }

#define  NMAX_TMP  200   /* temporary limit, to avoid using  malloc()  */
#define  VALID   0       /* see below, under  membership_list[][2]     */
#define  TANGLE  1       /* see below, under  membership_list[][2]     */
#define  N_ITER  4
#define  LARGE_NUMBER  1000000000.0   /* to start search for smallest number */

/*-----------------------------------------------------------------------------
 *  cp_NMAX_bodies  --  makes a copy of a bodies in standard newton0 form
 *                 and padd memory allocation to contain room for NMAX bodies.
 *                 accepts: old_bod: pointer to bodies;
 *                            npart: the number of particles in that bodies.
 *                 returns: new_bod: pointer to the new copy
 *-----------------------------------------------------------------------------
 */
bodyptr  cp_NMAX_bodies(old_bod, npart, n_alloc)
bodyptr  old_bod;
int  npart, n_alloc;
    {
    bodyptr  new_bod;
    bodyptr  new_i;      /* points to an individual particle of new_bod */
    bodyptr  old_i;      /* points to an individual particle of old_bod */
    
    new_bod = mk_bodies(n_alloc);
    
    for (old_i = old_bod, new_i = new_bod; old_i - old_bod < npart;
                                                              old_i++, new_i++)
        {
        SETV(Pos(new_i), Pos(old_i));
        SETV(Vel(new_i), Vel(old_i));
	Mass(new_i) = Mass(old_i);
	Pot(new_i) = Pot(old_i);
        SETV(Acc(new_i), Acc(old_i));
        }
    return(new_bod);
    }

/*-----------------------------------------------------------------------------
 *  isbinary  --  returns TRUE is the two bodies are bound (E_kin < |E_pot|)
 *-----------------------------------------------------------------------------
 */
local bool  isbinary(body1, body2)
bodyptr  body1, body2;
    {
    real delta_r, delta_v, m_sum;
    
    delta_r = distv(Pos(body1), Pos(body2));
    delta_v = distv(Vel(body1), Vel(body2));
    m_sum = Mass(body1) + Mass(body2);  /* sum because of reduced mass below */

    if (delta_r * delta_v * delta_v < 2.0 * m_sum)        /* bound ? */
	return(TRUE);
    else
	return(FALSE);
    }

/*-----------------------------------------------------------------------------
 *  is_ancestor_or_self  --  returns TRUE if body1 is a composite object
 *                             with body2 as one of its components or (sub...)
 *                             components  (body1 is #i ; body2 is #j).
 *-----------------------------------------------------------------------------
 */
local bool  is_ancestor_or_self(i, j, nbody, membership_list)
int  i, j, nbody;
int  membership_list[NMAX_TMP][3];
    {
    if (i == j)
	return(TRUE);
    else if (i < nbody)
	return(FALSE);
    else if (is_ancestor_or_self(membership_list[i][0], j, nbody,
				   membership_list))
	return(TRUE);
    else if (is_ancestor_or_self(membership_list[i][1], j, nbody,
				   membership_list))
	return(TRUE);
    else
	return(FALSE);
    }

/*-----------------------------------------------------------------------------
 *  init_binary  --  initializes the new binary  body0  which has as components
 *                   body0  and  body1  (which may themselves be either single
 *                   stars or binaries or even more composite objects).
 *                   Also computes a and e; stores those values where desired.
 *-----------------------------------------------------------------------------
 */
local void  init_binary(body0, body1, body2, aptr, eptr)
bodyptr  body0, body1, body2;
realptr  aptr, eptr;              /* semi major axis and eccentricity, resp. */
    {
    real  delta_r, delta_v, m_sum;
    real  r_rel[NDIM], v_rel[NDIM], r_out_v[NDIM];
    real  lever1[NDIM], lever2[NDIM], lever_sum[NDIM];
    real  r_out_v_squared;

    delta_r = distv(Pos(body1), Pos(body2));
    delta_v = distv(Vel(body1), Vel(body2));
    m_sum = Mass(body1) + Mass(body2);

    *aptr = 1.0 / (2.0/delta_r - delta_v*delta_v / m_sum);

    SUBV(r_rel, Pos(body1), Pos(body2));
    SUBV(v_rel, Vel(body1), Vel(body2));
    r_out_v[0] = r_rel[1] * v_rel[2] - r_rel[2] * v_rel[1];
    r_out_v[1] = r_rel[2] * v_rel[0] - r_rel[0] * v_rel[2];
    r_out_v[2] = r_rel[0] * v_rel[1] - r_rel[1] * v_rel[0];
    DOTVP(r_out_v_squared, r_out_v, r_out_v);
    *eptr = sqrt(1.0 - (r_out_v_squared / (m_sum * *aptr)));

    MULVS(lever1, Pos(body1), Mass(body1));
    MULVS(lever2, Pos(body2), Mass(body2));
    ADDV(lever_sum, lever1, lever2);
    DIVVS(Pos(body0), lever_sum, m_sum);
    MULVS(lever1, Vel(body1), Mass(body1));
    MULVS(lever2, Vel(body2), Mass(body2));
    ADDV(lever_sum, lever1, lever2);
    DIVVS(Vel(body0), lever_sum, m_sum);
    Mass(body0) = m_sum;
    }

/*-----------------------------------------------------------------------------
 *  scatter_diagnostics  --  analyse the state of a scattering experiment, by
 *                           decomposing the group of all particles in mutually
 *                           unbound subgroups, and classifying each subgroup
 *                           as a single star, binary, or k-body clump, k>2.
 *                 init_flag: not used right now.
 *-----------------------------------------------------------------------------
 */
local void  scatter_diagnostics(the_state, init_flag)
stateptr  the_state;
bool  init_flag;
    {
    int  i, j, j_nearest;
    int  nbody;
    int  nobject;
    int  present_nobject, previous_nobject;
    int  iter_count;
    real  nearest_neighbor_distance;
    bodyptr  body_i;
    bodyptr  body_j;
    bodyptr  bodies;
    bodyptr  objects;
    bodyptr  new_binary;
    int  membership_list[NMAX_TMP][3];
    real  delta_r;
    real  semi_major_axis[NMAX_TMP];
    real  eccentricity[NMAX_TMP];
    systptr  sys;

    if (NDIM != 3)
        error("micro_diagnostics: NDIM = %d != 3 not implemented\n", NDIM);
    for (i = 0; i < NMAX_TMP; i++)
	membership_list[i][2] = VALID;

    sys = System(the_state);
    bodies = Bodies(sys);
    nbody = Nbody(sys);
    if (nbody > NMAX_TMP)
        error("micro_diagnostics: nbody = %d > NMAX_TMP = %d , sorry!\n",
               nbody, NMAX_TMP);
    objects = cp_NMAX_bodies(bodies, nbody, NMAX_TMP);
    nobject = nbody;
    previous_nobject = 0;

/*
 * first pass: find elementary binaries, made up of single stars:
 */

    iter_count=0;
    while (iter_count++ < N_ITER && nobject != previous_nobject)
                                       /* hack to get started; then to check */
        {
        present_nobject = nobject;

        for (body_i = objects + previous_nobject, i = previous_nobject;
             i < present_nobject; body_i++, i++)
	    {
	    nearest_neighbor_distance = LARGE_NUMBER;
	    for (body_j = objects, j = 0; j < present_nobject; body_j++, j++)
	        if (!is_ancestor_or_self(i, j, nbody, membership_list))
		    {
                    delta_r = distv(Pos(body_i), Pos(body_j));
	            if (delta_r < nearest_neighbor_distance)
		        {
		        nearest_neighbor_distance = delta_r;
		        j_nearest = j;
		        }
		    }
	    j = j_nearest;
            body_j = objects + j;
	    if ( (i > j) && (isbinary(body_i, body_j) == TRUE) )
		{
                if (nobject > NMAX_TMP - 1)
		    error("micro_diagnostics: too many binaries, sorry!\n");
		membership_list[nobject][0] = body_i - objects;
		membership_list[nobject][1] = body_j - objects;
		new_binary = objects + nobject;
		init_binary(new_binary, body_i, body_j,
			    &semi_major_axis[nobject], &eccentricity[nobject]);
		nobject += 1;
		}
	    }
        previous_nobject = present_nobject;    
        }


/*
 * second pass: weed out overlapping binaries:
 */

/*  for (body_i = objects + nobject - 1, i = nobject - 1; i >= nbody;
 *	 body_i--, i--)
 *      if (membership_list[i][2] == VALID)
 *	    for (body_j = body_i - 1, j = i - 1; j >= nbody; body_j--, j--)
 *	        if (membership_list[i][0] == membership_list[j][0] ||
 *	            membership_list[i][0] == membership_list[j][1] ||
 *	            membership_list[i][1] == membership_list[j][0] ||
 *	            membership_list[i][1] == membership_list[j][1])
 *		    membership_list[j][2] = TANGLE;
 */

/*
 * (temporary) spot where output is printed:
 */
    printf("\n");    

    for (body_i = objects + nbody, i = nbody; i < nobject; body_i++, i++)
	if (membership_list[i][2] == VALID)
	    {
	    printf("t = %lf  ", Tnow(sys));
	    print_object(i, nbody, membership_list);
	    printf(" : a = %lf  e = %lf\n",
		   semi_major_axis[i], eccentricity[i]);
	    }
    free(objects);
    }

/*-----------------------------------------------------------------------------
 *  micro_diagnostics  --  to obtain detailed information concerning binaries,
 *                         escapers, etc.
 *                         note: it is not very efficient to recompute
 *                               potential and kinetic energies here;
 *                               we will later worry about efficiency, though.
 *                         init_flag: not used right now.
 *                         note: the members of binaries[i] are
 *                               membership_list[i][0] & membership_list[i][1].
 *                               stars are numbered from 0 to nbody-1;
 *                               binaries are numbered from 0 to nbinary-1;
 *                               in a hierarchical triple, the inner binary #j
 *                               is classified as a "star" with number
 *                               nbody + j .
 *                         note: membership_list[i][2] = VALID for OK binaries,
 *                               membership_list[i][2] = TANGLE for invalid
 *                               hierarchical triples, where the inner 
 *                               apocenter exceeds the outer pericenter.
 *-----------------------------------------------------------------------------
 */
local void  micro_diagnostics(the_state, init_flag)
stateptr  the_state;
bool  init_flag;
    {
    int  i, j, j_nearest;
    int  nbody;
    int  nobject;
    int  present_nobject, previous_nobject;
    int  iter_count;
    real  nearest_neighbor_distance;
    bodyptr  body_i;
    bodyptr  body_j;
    bodyptr  bodies;
    bodyptr  objects;
    bodyptr  new_binary;
    int  membership_list[NMAX_TMP][3];
    real  delta_r;
    real  semi_major_axis[NMAX_TMP];
    real  eccentricity[NMAX_TMP];
    systptr  sys;

    if (NDIM != 3)
        error("micro_diagnostics: NDIM = %d != 3 not implemented\n", NDIM);
    for (i = 0; i < NMAX_TMP; i++)
	membership_list[i][2] = VALID;

    sys = System(the_state);
    bodies = Bodies(sys);
    nbody = Nbody(sys);
    if (nbody > NMAX_TMP)
        error("micro_diagnostics: nbody = %d > NMAX_TMP = %d , sorry!\n",
               nbody, NMAX_TMP);
    objects = cp_NMAX_bodies(bodies, nbody, NMAX_TMP);
    nobject = nbody;
    previous_nobject = 0;

/*
 * first pass: find elementary binaries, made up of single stars:
 */

    iter_count=0;
    while (iter_count++ < N_ITER && nobject != previous_nobject)
                                       /* hack to get started; then to check */
        {
        present_nobject = nobject;

        for (body_i = objects + previous_nobject, i = previous_nobject;
             i < present_nobject; body_i++, i++)
	    {
	    nearest_neighbor_distance = LARGE_NUMBER;
	    for (body_j = objects, j = 0; j < present_nobject; body_j++, j++)
	        if (!is_ancestor_or_self(i, j, nbody, membership_list))
		    {
                    delta_r = distv(Pos(body_i), Pos(body_j));
	            if (delta_r < nearest_neighbor_distance)
		        {
		        nearest_neighbor_distance = delta_r;
		        j_nearest = j;
		        }
		    }
	    j = j_nearest;
            body_j = objects + j;
	    if ( (i > j) && (isbinary(body_i, body_j) == TRUE) )
		{
                if (nobject > NMAX_TMP - 1)
		    error("micro_diagnostics: too many binaries, sorry!\n");
		membership_list[nobject][0] = body_i - objects;
		membership_list[nobject][1] = body_j - objects;
		new_binary = objects + nobject;
		init_binary(new_binary, body_i, body_j,
			    &semi_major_axis[nobject], &eccentricity[nobject]);
		nobject += 1;
		}
	    }
        previous_nobject = present_nobject;    
        }


/*
 * second pass: weed out overlapping binaries:
 */

/*  for (body_i = objects + nobject - 1, i = nobject - 1; i >= nbody;
 *	 body_i--, i--)
 *      if (membership_list[i][2] == VALID)
 *	    for (body_j = body_i - 1, j = i - 1; j >= nbody; body_j--, j--)
 *	        if (membership_list[i][0] == membership_list[j][0] ||
 *	            membership_list[i][0] == membership_list[j][1] ||
 *	            membership_list[i][1] == membership_list[j][0] ||
 *	            membership_list[i][1] == membership_list[j][1])
 *		    membership_list[j][2] = TANGLE;
 */

/*
 * (temporary) spot where output is printed:
 */
    printf("\n");    

    for (body_i = objects + nbody, i = nbody; i < nobject; body_i++, i++)
	if (membership_list[i][2] == VALID)
	    {
	    printf("t = %lf  ", Tnow(sys));
	    print_object(i, nbody, membership_list);
	    printf(" : a = %lf  e = %lf\n",
		   semi_major_axis[i], eccentricity[i]);
	    }
    free(objects);
    }

/*-----------------------------------------------------------------------------
 *  print_object  --  prints a (possibly recursively composite) object
 *-----------------------------------------------------------------------------
 */
local void  print_object(i, nbody, membership_list)
int  i, nbody;
int  membership_list[NMAX_TMP][3];
    {
    if (i < nbody)
	printf("%d", i);
    else
	{
	printf("(");
	print_object(membership_list[i][0], nbody, membership_list);
	printf(" ");
	print_object(membership_list[i][1], nbody, membership_list);
	printf(")");
	}
    }

/* endof: diagnose.c */
