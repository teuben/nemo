/* scatdiagnose.c - check_energy, compute_energy, compute_total_energy_only,
                    diagnostics, final_diagnostics, grinding_halt,
                    scatter_diagnostics */
/*
 *  scatdiagnose.c: diagnostics module for  scatter.c
 *
 *      Oct 1988  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
#include  "newton0.h"

extern bool  debug;       /* flag for output control, defined in  scatter.c  */

/*-----------------------------------------------------------------------------
 *  diagnostics  --  performs a generic set of diagnostics, and initialize
 *                   the diagnostics substructure (remember: this could not
 *                   be set earlier, together with the specs and ctrls in the
 *                   file containing the main() procedure, since that would
 *                   have required  energy diagnostics, which would have
 *                   necessitated revving the engines, which could not have
 *                   been done before the other three substructures were set).
 *-----------------------------------------------------------------------------
 */
int  diagnostics(the_state, init_flag)
stateptr  the_state;
bool  init_flag;
    {
    int  diagnosis;
    int  scatter_diagnostics();
    diagptr  diags;

    diags = Diags(the_state);

    if (init_flag == TRUE)               /* should perhaps go in  newton0.c  */
        DIAGnsteps(diags) = DIAGde_nsteps(diags) = 0;

    Hier_string(diags)[0] = NULL;

    check_energy(the_state, init_flag);
/*
 *  if (streq("microscopic", Diagnostics(Specs(the_state))))
 *	micro_diagnostics(the_state, init_flag);
 */
    if (streq("scatter", Diagnostics(Specs(the_state))))
	diagnosis = scatter_diagnostics(the_state, init_flag);

    /* .... */                   /* more to follow in future implementations */

    return(diagnosis);
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
 *  grinding_halt  --  halt the present integration, supply information to the
 *                     proper output channels about the error condition which
 *                     prompted this premature halt, and exit as graceful as
 *                     possible under the circumstances.
 *-----------------------------------------------------------------------------
 */
void  grinding_halt(the_state)
stateptr  the_state;
    {
    error("\ngrinding_halt: not yet implemented\n",0.0);
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

    sys = System(the_state);
    diags = Diags(the_state);

    if (init_flag == FALSE)
	{
        DIAGtime(diags)[PREVIOUS] = DIAGtime(diags)[CURRENT];
        DIAGetot(diags)[PREVIOUS] = DIAGetot(diags)[CURRENT];
        DIAGekin(diags)[PREVIOUS] = DIAGekin(diags)[CURRENT];
        DIAGepot(diags)[PREVIOUS] = DIAGepot(diags)[CURRENT];
        }

#ifndef REGULARIZATION
    compute_energy(&(DIAGekin(diags)[CURRENT]),
		   &(DIAGepot(diags)[CURRENT]), sys, diags);
#else
    DIAGekin(diags)[CURRENT] = (Reghamiltonian(Regsystem(the_state))
                                + Reglagrangian(Regsystem(the_state))) / 2.0;
    DIAGepot(diags)[CURRENT] = (Reghamiltonian(Regsystem(the_state))
                                - Reglagrangian(Regsystem(the_state))) / 2.0;
#endif

    if (init_flag == TRUE)
	{
        DIAGtime(diags)[INITIAL] = DIAGtime(diags)[PREVIOUS] =
                                          DIAGtime(diags)[CURRENT] = Tnow(sys);
        DIAGekin(diags)[INITIAL] = DIAGekin(diags)[CURRENT];
        DIAGepot(diags)[INITIAL] = DIAGepot(diags)[CURRENT];
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
 *-----------------------------------------------------------------------------
 */
void  compute_energy(ekinptr, epotptr, sys, diags)
realptr  ekinptr;
realptr  epotptr;
systptr  sys;
diagptr  diags;
    {
    int  nbody;
    real  kinetic_energy;
    real  velocity_squared;
    real  sum_kin_energy;
    real  sum_pot_energy;
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

    *ekinptr = 0.5 * sum_kin_energy;   /* corrects for the factor 2 left out */
    *epotptr = 0.5 * sum_pot_energy;   /* corrects for double counting pairs */
    }

/*-----------------------------------------------------------------------------
 *  scatter_diagnostics  --  provides diagnostic information to describe the
 *                           hierarchical structure of a few-body system,
 *                           to be used in a few-body scattering experiment.
 *-----------------------------------------------------------------------------
 */
local int  scatter_diagnostics(the_state, init_flag)
stateptr  the_state;
bool  init_flag;
    {
    int  nbody;
    int  diagnosis;
    int  find_hierarchical_structure();
    bodyptr  bodies;
    systptr  sys;
    diagptr  diags;

    if (NDIM != 3)
        error("scatter_diagnostics: NDIM = %d != 3 not implemented\n", NDIM);

    sys = System(the_state);
    bodies = Bodies(sys);
    nbody = Nbody(sys);
    diags = Diags(the_state);

    diagnosis = find_hierarchical_structure(bodies, nbody, diags, init_flag,
					    Tnow(sys));
    return(diagnosis);
    }

/*-----------------------------------------------------------------------------
 *  macros for the hierarchical decomposition of a few-body system:
 *      the following macro definitions govern the bookkeeping;
 *      they should migrate to their own header file in due time.
 *             NOTE: only one-digit star numbering implemented, i.e.
 *                   only valid for an N-body sub-system with N < 11
 *-----------------------------------------------------------------------------
 */
#define  NBODY_MAX    10                     /* maximum number of particles  */
#define  NNODE_MAX    (2 * NBODY_MAX - 1)    /* maximum number of nodes      */

/*-----------------------------------------------------------------------------
 *  the following macros are used to access
 *                int  admin_table[NNODE_MAX][N_ADMIN];
 *  the entrees for the second variable are:
 *                                           number of daughter nodes   
 *                                           covered or visible         
 *                                           stablve, unstable, barely unbound
 *                                              or strongly unbound
 *                                           index of first daughter  
 *                                           index of second daughter 
 *                                           ...
 *                                           index of last daughter   
 *-----------------------------------------------------------------------------
 */
#define  N_ADMIN      (3 + NBODY_MAX)

#define  Ndaughter(ptr, offset)          ((ptr)[(offset)][0])

#define  Visibility(ptr, offset)         ((ptr)[(offset)][1])
#define  COVERED      1
#define  VISIBLE      2

#define  Stability(ptr, offset)          ((ptr)[(offset)][2])
#define  STABLE            1                          /* bound, and stable   */
#define  UNSTABLE          2                          /* bound, and unstable */
#define  BARELY_UNBOUND    3                          /* barely unbound      */
#define  STRONGLY_UNBOUND  4                          /* barely unbound      */
#define  UNKNOWN           5

#define  Daughter(ptr, offset, i)        ((ptr)[(offset)][3 + (i)])

/*-----------------------------------------------------------------------------
 *  the following macros are used to access
 *                real  struct_table[NNODE_MAX][N_STRUCT];
 *  the entrees of the second variable are:
 *                                          radius
 *                                          a: semimajor axis
 *                                          e: eccentricity
 *-----------------------------------------------------------------------------
 */
#define  N_STRUCT     3

#define  Radius(ptr, offset)    ((ptr)[(offset)][0])
#define  Apair(ptr, offset)     ((ptr)[(offset)][1])
#define  Epair(ptr, offset)     ((ptr)[(offset)][2])

#define  INVERSE_FREQUENCY_OF_DEBUG_OUTPUTS  5

/*-----------------------------------------------------------------------------
 *  find_hierarchical_structure  --  
 *-----------------------------------------------------------------------------
 */
local int  find_hierarchical_structure(bodies, nbody, diags, init_flag, tnow)
bodyptr  bodies;
int  nbody;
diagptr  diags;
bool  init_flag;
real  tnow;
    {
    int  nnode;
    int *nnodeptr;
    body  nodes[NNODE_MAX];
    int  admin_table[NNODE_MAX][N_ADMIN];
    permanent int  debug_output_counter;
    real  struct_table[NNODE_MAX][N_STRUCT];
    real  dist_table[NNODE_MAX][NNODE_MAX];
    bool  news;
    bool  anything_new_to_report();
    bool  stop_the_show();

    if (nbody > NBODY_MAX)
	error("find_hierarchical_structure: nbody = %d > NBODY_MAX = %d\n",
	      nbody, NBODY_MAX);

    if (init_flag == TRUE && debug == TRUE)
	debug_output_counter = INVERSE_FREQUENCY_OF_DEBUG_OUTPUTS;

    nnodeptr = &nnode;

    init_nodes(nodes, nnodeptr, bodies, nbody);

    init_dist_table(nodes, nnode, dist_table);
    init_admin_table(admin_table);
    init_struct_table(struct_table);

    iterate_find_new_node(nodes, nnodeptr, admin_table, struct_table,
			  dist_table, 2);

    iterate_dissolve_visible_strongly_unbound_nodes(nodes, nnodeptr,
						    admin_table, struct_table,
						    dist_table);

    news = anything_new_to_report(admin_table, nnode, init_flag);

    if (news == TRUE)
	{
/*
 * the following is meant for human reading,
 * and can be commented out for machine use.
 */
	print_hierarchical_structure(nnode, admin_table, struct_table,
				     dist_table, tnow);
/*
 * the following is meant for machine reading,
 * and can be commented out for human use.
 */
	unique_decomposition(diags, nnode, admin_table, struct_table);

	return(REPORT_RUN);
	}
    else
	{
	if (stop_the_show(nodes, nnode, admin_table, struct_table, dist_table,
			  tnow, init_flag) == TRUE)
	    {
/*
 * the following is meant for human reading,
 * and can be commented out for machine use.
 */
	    print_hierarchical_structure(nnode, admin_table, struct_table,
					 dist_table, tnow);
/*
 * the following is meant for machine reading,
 * and can be commented out for human use.
 */
	    unique_decomposition(diags, nnode, admin_table, struct_table);

	    return(HALT_RUN);
	    }
	else
	    {
	    if (debug == TRUE && --debug_output_counter == 0)
		{
	        print_hierarchical_structure(nnode, admin_table, struct_table,
					     dist_table, tnow);
		debug_output_counter = INVERSE_FREQUENCY_OF_DEBUG_OUTPUTS;
		}
	    return(SILENT_RUN);
	    }
	}
    }

/*-----------------------------------------------------------------------------
 *  iterate_find_new_node  --  
 *-----------------------------------------------------------------------------
 */
local void  iterate_find_new_node(nodes, nnodeptr, admin_table, struct_table,
				  dist_table, k)
bodyptr  nodes;
int  *nnodeptr;
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
real  dist_table[NNODE_MAX][NNODE_MAX];
int  k;
    {
    int  n_visible_nodes;
    int  count_visible_nodes();
    bool  find_new_node();

    n_visible_nodes = count_visible_nodes(admin_table, *nnodeptr);

    if (k > n_visible_nodes)
	return;

    if (find_new_node(nodes, nnodeptr, admin_table, struct_table, dist_table,
		      k) == TRUE)
	iterate_find_new_node(nodes, nnodeptr, admin_table, struct_table,
			      dist_table, 2);
    else
	iterate_find_new_node(nodes, nnodeptr, admin_table, struct_table,
			      dist_table, k+1);
    }

/*-----------------------------------------------------------------------------
 *  find_new_node  --  
 *-----------------------------------------------------------------------------
 */
local bool  find_new_node(nodes, nnodeptr, admin_table, struct_table,
			  dist_table, k)
bodyptr  nodes;
int  *nnodeptr;
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
real  dist_table[NNODE_MAX][NNODE_MAX];
int  k;
    {
    int  i;
    int  n;                                       /* number of visible nodes */
    int  tuple[NBODY_MAX];
    int  ordered_tuple[NBODY_MAX];
    int  count_visible_nodes();
    bool  next_tuple();
    bool  is_new_node();

    n = count_visible_nodes(admin_table, *nnodeptr);
    if (k > n)
	return(FALSE);

    init_tuple(tuple, admin_table, n);
    for (i = 0; i < n; i++)
        ordered_tuple[i] = tuple[i];

    do
        if (is_new_node(tuple, nodes, nnodeptr, admin_table, struct_table,
			dist_table, k) == TRUE)
	    return(TRUE);
    while
	(next_tuple(tuple, k, n, ordered_tuple) == TRUE);
	
    return(FALSE);
    }

/*-----------------------------------------------------------------------------
 *  is_new_node  --  checks wether a node is ilsolated
 *              note: old version:
 *                   if this particular tuple is not isolated, or does not form
 *                     a bound subsystem, then the value  FALSE  is returned;
 *                   if the tuple is isolated and bound,  TRUE  is returned
 *                     after the following actions have been carried out:
 *                        1) the members of this tuple are covered,
 *                        2) a new visible node is assigned to this tuple,
 *                        3) both tables are updated.
 *-----------------------------------------------------------------------------
 */
local bool  is_new_node(tuple, nodes, nnodeptr, admin_table, struct_table,
			dist_table, k)
int  tuple[NBODY_MAX];
bodyptr  nodes;
int  *nnodeptr;
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
real  dist_table[NNODE_MAX][NNODE_MAX];
int  k;
    {
    real  radius;
    bool  is_approximately_isolated_tuple();
    bool  is_isolated_tuple();

    if (k < 2)
	error("is_new_node: k = %d < 2\n", k);

    if (is_approximately_isolated_tuple(tuple, *nnodeptr, admin_table,
					dist_table, k) == FALSE)
	return(FALSE);

    if (is_isolated_tuple(tuple, nodes, *nnodeptr, admin_table, struct_table,
			  dist_table, &radius, k) == FALSE)
	return(FALSE);

    install_node(tuple, nodes, nnodeptr, admin_table, struct_table, dist_table,
		 radius, k);

    return(TRUE);
    }

/*-----------------------------------------------------------------------------
 *  is_approximately_isolated_tuple  --  a tuple is defined to be approximately
 *                                       isolated if the minimum distance
 *                                       between a member of the tuple and a
 *                                       non-member exceeds the maximum
 *                                       distance between a pair of members.
 *-----------------------------------------------------------------------------
 */
local bool  is_approximately_isolated_tuple(tuple, nnode, admin_table,
					    dist_table, k)
int  tuple[NBODY_MAX];
int  nnode;
int  admin_table[NNODE_MAX][N_ADMIN];
real  dist_table[NNODE_MAX][NNODE_MAX];
int  k;
    {
    int  i, j, jt;
    real  max_internal_pair_distance;
    int  visible_non_tuple[NBODY_MAX];
    int  n_visible_non_tuple;

    if (k < 2)
	error("is_approximately_isolated_tuple: k = %d < 2\n", k);

    max_internal_pair_distance = 0.0;

    for (i = 0; i < k - 1; i++)
	for (j = i + 1; j < k; j++)
	    if (dist_table[tuple[i]][tuple[j]] > max_internal_pair_distance)
		max_internal_pair_distance = dist_table[tuple[i]][tuple[j]];

    jt = 0;
    j = 0;
    for (i = 0; i < nnode; i++)
	if (Visibility(admin_table, i) == VISIBLE)
	    if (jt >= k)                            /* i.e.:  i > tuple[k-1] */
		visible_non_tuple[j++] = i;
	    else if (i < tuple[jt])
		visible_non_tuple[j++] = i;
	    else if (i == tuple[jt])
		jt++;
            else
		error("is_approximately_isolated_tuple: non-visible tuple?\n");

    n_visible_non_tuple = j;
    
    for (i = 0; i < n_visible_non_tuple; i++)
	for (j = 0; j < k; j++)
	    if (dist_table[visible_non_tuple[i]][tuple[j]]
		< max_internal_pair_distance)
		return(FALSE);

    return(TRUE);    
    }

/*-----------------------------------------------------------------------------
 *  is_isolated_tuple  -- two objects, with radius Ri and Rj, are considered
 *                        isolated, if each point inside a sphere belonging
 *                        to an object is closer to any point in its own sphere
 *                        than to any point in the other sphere. This implies
 *                        that the distance between the sphere surfaces should
 *                        be at least as large as the largest of the two 
 *                        diameters, i.e. 2.0*max(Ri, Rj), and the distance
 *                        between the two centers of mass of the objects should
 *                        be at least  
 *                                   Ri + Rj + 2.0*max(Ri, Rj)
 *-----------------------------------------------------------------------------
 */
local bool  is_isolated_tuple(tuple, nodes, nnode, admin_table, struct_table,
			      dist_table, radiusptr, k)
int  tuple[NBODY_MAX];
bodyptr  nodes;
int  nnode;
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
real  dist_table[NNODE_MAX][NNODE_MAX];
realptr  radiusptr;
int  k;
    {
    int  i, j, jt;
    real  radius;
    real  radius_i;
    real  tuple_radius();
    real  com_pos[NDIM];
    real  distance_to_com;
    int  visible_non_tuple[NBODY_MAX];
    int  n_visible_non_tuple;

    radius = tuple_radius(tuple, nodes, struct_table, dist_table, com_pos, k);
    *radiusptr = radius;

    jt = 0;
    j = 0;

    for (i = 0; i < nnode; i++)
	if (Visibility(admin_table, i) == VISIBLE)
	    if (jt >= k)                            /* i.e.:  i > tuple[k-1] */
                visible_non_tuple[j++] = i;
	    else if (i < tuple[jt])
		visible_non_tuple[j++] = i;
	    else if (i == tuple[jt])
		jt++;
            else
		error("is_isolated_tuple: non-visible tuple?\n");

    n_visible_non_tuple = j;
    for (i = 0; i < n_visible_non_tuple; i++)
	{
	DISTV(distance_to_com, Pos(nodes + visible_non_tuple[i]), com_pos);
	radius_i = Radius(struct_table, visible_non_tuple[i]);
	if (distance_to_com < radius + radius_i + 2.0 * MAX(radius, radius_i))
	    return(FALSE);
	}

    return(TRUE);
    }

/*-----------------------------------------------------------------------------
 *  how_bound_tuple  --  
 *-----------------------------------------------------------------------------
 */
local void  how_bound_tuple(tuple, nodes, dist_table, k, isbound,
			    isbarelyunbound)
int  tuple[NBODY_MAX];
bodyptr  nodes;
real  dist_table[NNODE_MAX][NNODE_MAX];
int  k;
bool *isbound;
bool *isbarelyunbound;
    {
    real  kinetic_energy;
    real  abs_potential_energy;
    real  tuple_kin_energy();
    real  tuple_pot_energy();

    kinetic_energy = tuple_kin_energy(tuple, nodes, k);
    abs_potential_energy = -tuple_pot_energy(tuple, nodes, dist_table, k);

    *isbound = *isbarelyunbound = FALSE;

    if (kinetic_energy < abs_potential_energy)
	*isbound = TRUE;
    else if (kinetic_energy < 1.5 * abs_potential_energy)
	*isbarelyunbound = TRUE;
    }

/*-----------------------------------------------------------------------------
 *  tuple_pot_energy  --  return the total potential energy of all the nodes in
 *                        the tuple, due to their mutual gravitational
 *                        attraction.
 *-----------------------------------------------------------------------------
 */
local real  tuple_pot_energy(tuple, nodes, dist_table, k)
int  tuple[NBODY_MAX];
bodyptr  nodes;
real  dist_table[NNODE_MAX][NNODE_MAX];
int  k;
    {
    int  i, j;
    real  pot_energy;

    pot_energy = 0.0;

    for (i = 0; i < k - 1; i++)
	for (j = i + 1; j < k; j++)
	    pot_energy -= ( Mass(nodes + tuple[i]) * Mass(nodes + tuple[j]) )
		          / dist_table[ tuple[i] ][ tuple[j] ];

    return(pot_energy);
    }

/*-----------------------------------------------------------------------------
 *  tuple_kin_energy  --  return the total kinetic energy of all the nodes in
 *                        the tuple, in the center of mass frame of the tuple.
 *-----------------------------------------------------------------------------
 */
local real  tuple_kin_energy(tuple, nodes, k)
int  tuple[NBODY_MAX];
bodyptr  nodes;
int  k;
    {
    int  i;
    real  velocity_squared;
    real  kin_energy;
    real  total_mass;
    real  mass_times_vel[NDIM];
    real  rel_vel[NDIM];                  /* velocity with respect to c.o.m. */
    real  com_vel[NDIM];
    bodyptr  body_i;

    total_mass = 0.0;
    CLRV(com_vel);
    for (i = 0; i < k; i++)
	{
        body_i = nodes + tuple[i];
	total_mass += Mass(body_i);
	MULVS(mass_times_vel, Vel(body_i), Mass(body_i));
	INCADDV(com_vel, mass_times_vel);
	}
    INCDIVVS(com_vel, total_mass);               /* the real c.o.m. velocity */

    kin_energy = 0.0;

    for (i = 0; i < k; i++)
	{
	body_i = nodes + tuple[i];
	SUBV(rel_vel, Vel(body_i), com_vel);
	DOTVP(velocity_squared, rel_vel, rel_vel);
	kin_energy += Mass(body_i) * velocity_squared;       /* 2 x too high */
	}

    kin_energy *= 0.5;                 /* corrects for the factor 2 left out */

    return(kin_energy);
    }

/*-----------------------------------------------------------------------------
 *  anything_new_to_report  --  returns TRUE if anything has changed;
 *                              returns FALSE if all has remained the same.
 *-----------------------------------------------------------------------------
 */
local bool  anything_new_to_report(new_admin_table, new_nnode, init_flag)
int  new_admin_table[NNODE_MAX][N_ADMIN];
int  new_nnode;
bool  init_flag;
    {
    int  i, j;
    permanent int  old_admin_table[NNODE_MAX][N_ADMIN];            /* static */
    permanent int  old_nnode;                                      /* static */
    bool  same_admin_table();

    if (init_flag == FALSE && new_nnode == old_nnode &&
	same_admin_table(new_admin_table, old_admin_table, new_nnode) == TRUE)
	return(FALSE);

    old_nnode = new_nnode;
    for (i = 0; i < old_nnode; i++)
	for (j = 0; j < N_ADMIN; j++)
	    old_admin_table[i][j] = new_admin_table[i][j];

    return(TRUE);
    }

/*-----------------------------------------------------------------------------
 *  same_admin_table  --  returns TRUE if both admin_table's are identical,
 *                        returns FALSE otherwise
 *-----------------------------------------------------------------------------
 */
local bool  same_admin_table(new_admin_table, old_admin_table, nnode)
int  new_admin_table[NNODE_MAX][N_ADMIN];
int  old_admin_table[NNODE_MAX][N_ADMIN];
int  nnode;
    { 
    int  i, j, n_daughter;

    for (i = 0; i < nnode; i++)
	{
	if ( (n_daughter = Ndaughter(new_admin_table, i)) !=
	     Ndaughter(old_admin_table, i) )
	    return(FALSE);
	if (Visibility(new_admin_table, i) != Visibility(old_admin_table, i))
	    return(FALSE);
	if (Stability(new_admin_table, i) != Stability(old_admin_table, i))
	    return(FALSE);
	for (j = 0; j < n_daughter; j++)
	    if (Daughter(new_admin_table,i,j) != Daughter(old_admin_table,i,j))
		return(FALSE);
	}

    return(TRUE);
    }    

#define   KINETIC_ENERGY_FACTOR           10.0
#define   NODE_SEPARATION_FACTOR          10.0
#define   SEVERAL_TIMES_LARGER_THAN_ONE    5.0
#define   BIGNUMBER                        1.e9
#define   DEBUG_HALTING_SPEEDUP_FRACTION   0.25

/*-----------------------------------------------------------------------------
 *  stop_the_show  --  returns TRUE if the scattering experiment has reached
 *                     the asymptotic future state, in which all fragments
 *                     move away from the interaction region, in such a way
 *                     that their cumulative projected future interaction with
 *                     each other is estimated to remain below a certain
 *                     threshold;
 *                     returns FALSE otherwise.
 *                note:
 *                     if debug == TRUE, the whole procedure has to be executed
 *                     in order to determine the fraction to which the halting
 *                     criterion has been satisfied; if debug == FALSE,
 *                     execution can be halted earlier, thereby saving some
 *                     time.
 *                     when debug == TRUE, the fractional to which the halt 
 *                     criterion has been satisfied is printed whenever it has
 *                     changed by more than 5%.
 *                NOTE:
 *                     The hierarchical check is still very 
 *                     preliminary, and is only implemented for
 *                     checking just one hierarchical triple.
 *-----------------------------------------------------------------------------
 */
local bool  stop_the_show(nodes, nnode, admin_table, struct_table, dist_table,
			  tnow, init_flag)
bodyptr  nodes;
int  nnode;
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
real  dist_table[NNODE_MAX][NNODE_MAX];
real  tnow;
bool  init_flag;
    {
    int  i, j;
    real  pos_dot_vel;
    real  kin_fact, pot_fact;    /* energy factors, note use of reduced mass */
    real  rel_pos[NDIM];
    real  rel_vel[NDIM];
    real  how_freely_escaping();
    real  halt_criterion;
    real  min_halt_criterion;
    real  local_separation_criterion;
    real  local_energy_criterion;
    real  global_energy_criterion;
    permanent real  last_min_halt_criterion;
    bodyptr  node_i;
    permanent bool  triple_quarantine;  /* checking out hierarchical triples */
    permanent real  trip_quar_t_end;    /*   up till  trip_quar_t_end        */
    real  quarantine_duration;
    bool  find_hier_triple();
    permanent int  old_outer_rep_node;
    permanent int  old_inner_rep_node;
    int  outer_rep_node;                /* representative body of outer node */
    int  inner_rep_node;                /* representative body of inner node */

    if (init_flag == TRUE)
        triple_quarantine = FALSE;

    if (init_flag == TRUE && debug == TRUE)
	last_min_halt_criterion = 0.0;
/*
 * check whether each node is either a single particle or a stable binary:
 */
    for (i = 0; i < nnode; i++)
	if (Stability(admin_table, i) != STABLE)
	    {
            triple_quarantine = FALSE;	
	    return(FALSE);
	    }
/*
 * make a redundant check, to see whether all particles are bound together.
 * This is redundant in the sense that a scattering experiment has an 
 * infinitesimal probability of giving rise to a stable bound system.
 */
    if (nnode < 2)
	error("stop_the_show: nnode = %d in a stable system?!?!?\n", nnode);

/*
 * check whether each visible node is moving outwards:
 */
    for (node_i = nodes; node_i - nodes < nnode; node_i++)
	if (Visibility(admin_table, node_i - nodes) == VISIBLE)
	    {
	    DOTVP(pos_dot_vel, Pos(node_i), Vel(node_i));
	    if (pos_dot_vel < 0.0)
	        {
                triple_quarantine = FALSE;	
	        return(FALSE);
	        }
	    }
/*
 * check whether in each pair of visible nodes the two nodes are moving
 * outwards from each other:
 */
    for (i = 0; i < nnode - 1; i++)
	if (Visibility(admin_table, i) == VISIBLE)
	    for (j = i+1; j < nnode; j++)
	        if (Visibility(admin_table, j) == VISIBLE)
		    {
		    SUBV(rel_pos, Pos(nodes + i), Pos(nodes + j));
		    SUBV(rel_vel, Vel(nodes + i), Vel(nodes + j));
		    DOTVP(pos_dot_vel, rel_pos, rel_vel);
		    if (pos_dot_vel < 0.0)
	                {
                        triple_quarantine = FALSE;	
	                return(FALSE);
	                }
		    }

    local_separation_criterion = BIGNUMBER;
    local_energy_criterion = BIGNUMBER;
    global_energy_criterion = BIGNUMBER;

/*
 * check that all visible nodes are very well separated: there distances 
 * being at least  NODE_SEPARATION_FACTOR  times as large as their sizes:
 */
    for (i = 0; i < nnode - 1; i++)
	if (Visibility(admin_table, i) == VISIBLE)
	    for (j = i+1; j < nnode; j++)
	        if (Visibility(admin_table, j) == VISIBLE)
		    {
		    halt_criterion = dist_table[i][j] / (NODE_SEPARATION_FACTOR
		       * (Radius(struct_table, i) + Radius(struct_table, j)));
		    if (local_separation_criterion > halt_criterion)
			local_separation_criterion = halt_criterion;
		    }

    if (debug == FALSE && local_separation_criterion < 1.0)
	{
        triple_quarantine = FALSE;	
	return(FALSE);
	}
/*
 * check that each visible node is moving out with a kinetic energy in
 * its absolute motion which is at least a factor  KINETIC_ENERGY_FACTOR 
 * larger than their potential energy:
 * NOTE: this criterion is sometimes much too strict: consider ionization of
 *       three particles, the middle of which remains near the center of mass;
 *       See below for a temporary solution, involving diluting the global
 *       energy criterion with the local one, taking the geometric mean,
 *       as a temporary stopgap.
 */
    for (i = 0; i < nnode; i++)
	if (Visibility(admin_table, i) == VISIBLE)
	    {
	    halt_criterion =
		(how_freely_escaping(nodes, nnode, admin_table, dist_table, i)
		/ KINETIC_ENERGY_FACTOR);
	    if (global_energy_criterion > halt_criterion)
		global_energy_criterion = halt_criterion;
	    }
/*
 * check that each pair of visible nodes is separating with a kinetic energy in
 * their relative motion which is at least a factor  KINETIC_ENERGY_FACTOR 
 * larger than their potential energy:
 */
    for (i = 0; i < nnode - 1; i++)
	if (Visibility(admin_table, i) == VISIBLE)
	    for (j = i+1; j < nnode; j++)
	        if (Visibility(admin_table, j) == VISIBLE)
		    {
		    SUBV(rel_vel, Vel(nodes + i), Vel(nodes + j));
		    DOTVP(kin_fact, rel_vel, rel_vel);
                    pot_fact = (Mass(nodes + i) + Mass(nodes + j))
			/ dist_table[i][j];  /* reduced mass, 2 body problem */
		    halt_criterion = (0.5*kin_fact / pot_fact) /
			KINETIC_ENERGY_FACTOR;
		    if (local_energy_criterion > halt_criterion)
			local_energy_criterion = halt_criterion;
		    }

/*
 * trick -- see above, under 'NOTE'
 */
    if (local_energy_criterion > SEVERAL_TIMES_LARGER_THAN_ONE)
	min_halt_criterion = MIN(local_energy_criterion,
				 local_separation_criterion);
    else if (local_energy_criterion > global_energy_criterion)
	{
	min_halt_criterion = sqrt(global_energy_criterion
				  * local_energy_criterion);
	if (min_halt_criterion > local_separation_criterion)
	    min_halt_criterion = local_separation_criterion;
	}
    else
	min_halt_criterion = MIN(local_energy_criterion,
				 local_separation_criterion);

    if (debug == TRUE && triple_quarantine == FALSE
	&& ABS(last_min_halt_criterion - min_halt_criterion) > 0.05)
	{
	printf("\nhalt_criterion satisfied for %2.0f%% at t = %lf\n",
	       100 * min_halt_criterion, tnow);
	last_min_halt_criterion = min_halt_criterion;

	printf("  local_separation_criterion satisfied for %2.0lf%%\n",
	       100 * local_separation_criterion);
	printf("  global_energy_criterion satisfied for %2.0lf%%\n",
	       100 * global_energy_criterion);
	printf("  local_energy_criterion satisfied for %2.0lf%%\n",
	       100 * local_energy_criterion);
	}

    if (debug == TRUE && triple_quarantine == FALSE
	&& min_halt_criterion > DEBUG_HALTING_SPEEDUP_FRACTION)
	{
	printf("\nhalt_criterion satisfied for %2.0lf%% at t = %lf ,\n",
	       100 * min_halt_criterion, tnow);
	printf("which exceeds the DEBUG_HALTING_SPEEDUP_FRACTION = %2.0lf%%\n",
	       100 * DEBUG_HALTING_SPEEDUP_FRACTION);

	printf("  local_separation_criterion satisfied for %2.0lf%%\n",
	       100 * local_separation_criterion);
	printf("  global_energy_criterion satisfied for %2.0lf%%\n",
	       100 * global_energy_criterion);
	printf("  local_energy_criterion satisfied for %2.0lf%%\n",
	       100 * local_energy_criterion);
	}

    if (debug == TRUE)
	{
	if  (min_halt_criterion < DEBUG_HALTING_SPEEDUP_FRACTION)
	    {
            triple_quarantine = FALSE;	
	    return(FALSE);
	    }
	}
    else if (min_halt_criterion < 1.0)
	    {
            triple_quarantine = FALSE;	
	    return(FALSE);
	    }

/*
 * check that all higher-order binaries have completed at least ... revolutions
 */
    if (find_hier_triple(nodes, nnode, &quarantine_duration, &outer_rep_node,
			 &inner_rep_node, admin_table, struct_table) == FALSE)
	return(TRUE);                                     /* no triple found */
    else if (triple_quarantine == FALSE)                  /* triple found    */
	{
	trip_quar_t_end = tnow + quarantine_duration;
	triple_quarantine = TRUE;
        old_outer_rep_node = outer_rep_node;
        old_inner_rep_node = inner_rep_node;
	if (debug == TRUE)
	    {
	    printf("\nstarting a hierarchical triple quarantine period\n");
	    printf("  from t = %lf until t = %lf\n", tnow, trip_quar_t_end);
	    }
	return(FALSE);
	}
    else if (outer_rep_node != old_outer_rep_node ||    /* different triple?         */
	     inner_rep_node != old_inner_rep_node)
	{
        triple_quarantine = FALSE;	        /* start again in next round */
	return(FALSE);
	}
    else if (tnow > trip_quar_t_end)
	return(TRUE);
    else
	return(FALSE);
    }

/*-----------------------------------------------------------------------------
 *  how_freely_escaping  --  returns the ratio of kinetic energy of the motion
 *                           of the center of mass of node i and the potential
 *                           energy of node i with respect to all other visible
 *                           nodes.
 *                           note: the present implementation does not take
 *                                 into account the sizes of the different
 *                                 nodes. The tidal forces could be estimated
 *                                 using that information in a future version.
 *-----------------------------------------------------------------------------
 */
local real  how_freely_escaping(nodes, nnode, admin_table, dist_table, i)
bodyptr  nodes;
int  nnode;
int  admin_table[NNODE_MAX][N_ADMIN];
real  dist_table[NNODE_MAX][NNODE_MAX];
    {
    int  j;
    real  kin_energy;       /* specific kinetic energy of escaping node      */
    real  pot_energy;       /* potential with respect to other visible nodes */

    if (Visibility(admin_table, i) != VISIBLE)
	error("how_freely_escaping: node i = %d is not VISIBLE\n", i);

    DOTVP(kin_energy, Vel(nodes + i), Vel(nodes + i));
    kin_energy *= 0.5;

    pot_energy = 0.0;
    for (j = 0; j < nnode; j++)
	if (Visibility(admin_table, j) == VISIBLE && j != i)
	    pot_energy += Mass(nodes + j) / dist_table[i][j];

    return( kin_energy / pot_energy );       /* positive convention for pot. */
    }

#define  QUARANTINE_NUMBER  100.0     /* in the unit of outer orbital period */
#define  QUARANTINE_RATIO     5.0     /* outer pericenter / inner apocenter  */

/*-----------------------------------------------------------------------------
 *  find_hier_triple  --  
 *                  NOTE: rough first try; only two pointers provided, to the
 *                        first outer and the first inner node encountered;
 *                        no provisions for more than one hierarchical triple,
 *                        nor distinctions between triples and larger systems.
 *                  NOTE: quarantine is suspended for now if the ratio of
 *                        outer pericenter and inner apocenter is too large.
 *                        This implies: for FOUR bodies only right now !!!
 *-----------------------------------------------------------------------------
 */
local bool  find_hier_triple(nodes, nnode, quarantine_duration_ptr,
			     outer_rep_node, inner_rep_node,
			     admin_table, struct_table)
bodyptr  nodes;
int  nnode;
real *quarantine_duration_ptr;
int *outer_rep_node;
int *inner_rep_node;
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
    {
    int  i, j;
    int  representative_body();
    int  inner_node;
    real  a_outer;
    real  m_total;
    real  outer_orbital_period;
    real  inner_apocenter;
    real  outer_pericenter;

    for (i = 0; i < nnode; i++)
	if (Ndaughter(admin_table, i) != 0)
	    for (j = 0; j < Ndaughter(admin_table, i); j++)
		if (Ndaughter(admin_table, Daughter(admin_table, i, j)) != 0)
		    {
		    inner_node = Daughter(admin_table, i, j);
		    *outer_rep_node = representative_body(admin_table, i);
		    *inner_rep_node = representative_body(admin_table,
						    Daughter(admin_table,i,j));
		    a_outer = Apair(struct_table, i);
		    m_total = Mass(nodes + i);
		    outer_pericenter = a_outer *
			(1.0 - Epair(struct_table, i));
		    inner_apocenter = Apair(struct_table, inner_node) *
			 (1.0 + Epair(struct_table, inner_node));
		    if (outer_pericenter / inner_apocenter >
			QUARANTINE_RATIO)
			*quarantine_duration_ptr = 0.0;
		    else
			{
			outer_orbital_period = TWO_PI *
			    sqrt(a_outer * a_outer * a_outer / m_total);
		        *quarantine_duration_ptr = outer_orbital_period *
			    QUARANTINE_NUMBER;
			}
		    return(TRUE);
		    }

    return (FALSE);
    }


#define  LARGE_INTEGER   100000
/*-----------------------------------------------------------------------------
 *  representative_body  --  returns the lowest of the numbers of the atoms
 *                           which make up the node  i .
 *-----------------------------------------------------------------------------
 */
local int  representative_body(admin_table, i)
int  admin_table[NNODE_MAX][N_ADMIN];
int  i;
    {
    int  j;
    int  number_of_daughters;
    int  daughter_j_representative;
    int  lowest_number;

    number_of_daughters = Ndaughter(admin_table, i);

    if (number_of_daughters != 0)
	{
	lowest_number = LARGE_INTEGER;
	for (j = 0; j < number_of_daughters; j++)
	    {
	    daughter_j_representative = representative_body(admin_table,
						    Daughter(admin_table,i,j));
	    if (lowest_number > daughter_j_representative)
		lowest_number = daughter_j_representative;
	    }
	return(lowest_number);
	}
    else
	return(i);
    }

/*-----------------------------------------------------------------------------
 *  init_nodes  --  initializes a flat tree, consisting of a row of bodies
 *                  i.e. each node is a real body.
 *-----------------------------------------------------------------------------
 */
local void  init_nodes(nodes, nnodeptr, bodies, nbody)
bodyptr  nodes;
int  *nnodeptr;
bodyptr  bodies;
int  nbody;
    {
    bodyptr  body_i, node_i;

    if (nbody > NNODE_MAX)
	error("init_nodes: nbody = %d > NNODE_MAX = %d\n", nbody, NNODE_MAX);

    for (body_i = bodies, node_i = nodes; body_i - bodies < nbody;
                                 	  body_i++, node_i++)
	cp_dynamics(node_i, body_i);

    *nnodeptr = nbody;
    }

/*-----------------------------------------------------------------------------
 *  cp_dynamics  --  copies the mass and orbital parameters, for a single body
 *-----------------------------------------------------------------------------
 */
local void  cp_dynamics(to_body, from_body)
bodyptr  to_body, from_body;
    {
    Mass(to_body) = Mass(from_body);
    SETV(Pos(to_body), Pos(from_body));
    SETV(Vel(to_body), Vel(from_body));
    }

/*-----------------------------------------------------------------------------
 *  init_dist_table  --  computes all node-node distances between their
 *                       centers of mass, and installs these values in the
 *                       distance table.
 *                       note: a moderate amount of inefficiency is allowed,
 *                             for the sake of simplicity.
 *                                 For example, the table is  n x n ,  rather 
 *                                 than n x (n-1) / 2  [symmetry is not used];
 *                                 and the distances, rather than the cheaper
 *                                 squared distances, are computed.
 *-----------------------------------------------------------------------------
 */
local void  init_dist_table(nodes, nnode, dist_table)
bodyptr  nodes;
int  nnode;
real  dist_table[NNODE_MAX][NNODE_MAX];
    {
    int  i, j;
    real  distance;

    for (i = 0; i < nnode; i++)
	dist_table[i][i] = 0.0;

    for (i = 0; i < nnode - 1; i++)
	for (j = i + 1; j < nnode; j++)
	    {
	    DISTV(distance, Pos(nodes + i), Pos(nodes + j));
	    dist_table[i][j] = dist_table[j][i] = distance;
	    }
    }

/*-----------------------------------------------------------------------------
 *  init_admin_table  --  
 *-----------------------------------------------------------------------------
 */
local void  init_admin_table(admin_table)
int  admin_table[NNODE_MAX][N_ADMIN];
    {
    int  i;

    for (i = 0; i < NNODE_MAX; i++)      /* initially the following holds:   */
	{
	Ndaughter(admin_table, i) = 0;         /* no daughter nodes present  */
	Visibility(admin_table, i) = VISIBLE;  /* and all nodes are visible  */
	Stability(admin_table, i) = STABLE;    /* and stable point particles */
	}
    }

/*-----------------------------------------------------------------------------
 *  init_struct_table  --  
 *-----------------------------------------------------------------------------
 */
local void  init_struct_table(struct_table)
real  struct_table[NNODE_MAX][N_STRUCT];
    {
    int  i;

    for (i = 0; i < NNODE_MAX; i++)
	Radius(struct_table, i) = 0.0;
    }

/*-----------------------------------------------------------------------------
 *  count_visible_nodes  --  
 *-----------------------------------------------------------------------------
 */
local int  count_visible_nodes(admin_table, nnode)
int  nnode;
int  admin_table[NNODE_MAX][N_ADMIN];
    {
    int  i;
    int  n_visible_nodes;

    n_visible_nodes = 0;

    for (i = 0; i < nnode; i++)
	if (Visibility(admin_table, i) == VISIBLE)
	    n_visible_nodes++;

    return( n_visible_nodes );
    }

/*-----------------------------------------------------------------------------
 *  init_tuple  -  find the indices of the first  n  visible nodes, and list
 *                 these indices at the beginning of the  tuple[]  array.
 *-----------------------------------------------------------------------------
 */
local void  init_tuple(tuple, admin_table, n)
int  tuple[NBODY_MAX];
int  admin_table[NNODE_MAX][N_ADMIN];
int  n;
    {
    int  i, j;

    if (n > NBODY_MAX)
	error("init_tuple: n = %d > NBODY_MAX = %d\n", n, NBODY_MAX);

    i = j = 0;

    while(i < NNODE_MAX && j < n)
	{
	if (Visibility(admin_table, i) == VISIBLE)
	    tuple[j++] = i;
	i++;
	}
    if (j != n)
	error("init_tuple: j = %d != n = %d\n", j, n);
    }
	
/*-----------------------------------------------------------------------------
 *  next_tuple  -  search for the next k-tuple from among the ordered list of
 *                 indices of visible nodes (provided in  ordered_tuple[] );
 *                 if there is a next k-tuple,
 *                     then replace the first k entrees in tuple[] with the
 *                     new k-tuple, and return TRUE;
 *                 else
 *                     return FALSE.
 *           note: the contents of tuple[] beyond the first k entrees, i.e.
 *                     tuple[k], tuple[k+1], ..., tuple[NBODY_MAX -1]
 *                 is undetermined, and should not be accessed.
 *-----------------------------------------------------------------------------
 */
local bool  next_tuple(tuple, k, n, ordered_tuple)
int  tuple[NBODY_MAX];
int  k;
int  n;
int  ordered_tuple[NBODY_MAX];
    {
    int  next_element();

    if (k > n)
	error("next_tuple: k = %d > n = %d\n", k, n);

    if (k < 1)
	return(FALSE);

    if (tuple[k-1] < ordered_tuple[n-1])
	{
	tuple[k-1] = next_element(tuple[k-1], ordered_tuple, n);
	return(TRUE);
	}
    else if (next_tuple(tuple, k-1, n-1, ordered_tuple) == TRUE)
	{
        tuple[k-1] = next_element(tuple[k-2], ordered_tuple, n);
	return(TRUE);
	}
    else
	return(FALSE);
    }

/*-----------------------------------------------------------------------------
 *  next_element  -  return the next element from an array tuple[]
 *-----------------------------------------------------------------------------
 */
local int  next_element(element, tuple, n)
int  element;
int  tuple[NBODY_MAX];
int  n;
    {
    int  i;

    for (i = 0; i < n - 1; i++)
	if (element == tuple[i])
	    return(tuple[i+1]);

    if (element == tuple[n-1])
	error("next_element: element provided is the last element of tuple\n");
    else
	error("next_element: element provided is not an element of tuple\n");
    }

/*-----------------------------------------------------------------------------
 *  tuple_radius  --  defined as the largest value of
 *                    (distance from the center of mass of tuple to a member
 *                     + internal radius of that member);
 *                    as a side effect, the center-of-mass positon is returned
 *                    in the argument com_pos[].
 *           new note:
 *                    for a k_tuple with  k = 2 but unbound, same treatment as
 *                    for k > 2.
 *           old note:
 *                    for a k_tuple with  k > 2  the instantaneous radius
 *                    is computed; 
 *                    for a pair (k = 2), the maximum radius is computed in the
 *                    approximation of an unperturbed ellipse,
 *                    i.e. the maximum value of
 *                    (distance from the center of mass of the binary to the
 *                    apocenter of a member + internal radius of that member).
 *-----------------------------------------------------------------------------
 */
local real  tuple_radius(tuple, nodes, struct_table, dist_table, com_pos, k)
int  tuple[NBODY_MAX];
bodyptr  nodes;
real  struct_table[NNODE_MAX][N_STRUCT];
real  dist_table[NNODE_MAX][NNODE_MAX];
real  com_pos[NDIM];
int  k;
    {
    int  i;
    real  member_radius;
    real  max_member_radius;
    real  total_mass;
    real  lever_arm_factor;          /* in case of a binary, i.e. k = 2      */
    real  mass_times_pos[NDIM];
    real  a, e;                      /* binary parameters for the case k = 2 */
    bodyptr  member;

    total_mass = 0.0;
    CLRV(com_pos);
    for (i = 0; i < k; i++)
	{
        member = nodes + tuple[i];
	total_mass += Mass(member);
	MULVS(mass_times_pos, Pos(member), Mass(member));
	INCADDV(com_pos, mass_times_pos);
	}
    INCDIVVS(com_pos, total_mass);               /* the real c.o.m. position */

    if (k == 2)
	get_binary_parameters(nodes + tuple[0], nodes + tuple[1], &a, &e);

    max_member_radius = 0.0;

    for (i = 0; i < k; i++)
	{
	if (k == 2 && a > 0.0)
	    {
	    if (i == 0)
		lever_arm_factor = Mass(nodes + tuple[1]) / total_mass;
	    else
		lever_arm_factor = Mass(nodes + tuple[0]) / total_mass;
	    member_radius = a * (1.0 + e) * lever_arm_factor;
	    }
	else
	    DISTV(member_radius, Pos(nodes + tuple[i]), com_pos);
	member_radius += Radius(struct_table, tuple[i]);
	if (member_radius > max_member_radius)
	    max_member_radius = member_radius;
	}

    return(max_member_radius);
    }

/*-----------------------------------------------------------------------------
 *  install_node  --  does all the bookkeeping needed for the introduction of
 *                    a new node.
 *                    note: incrementing the node counter should be postponed
 *                          to the last line, since all functions called
 *                          presume that the new node is  nodes[*nnodeptr] .
 *                     note:
 *                          update_struct_table()  has to be invoked before
 *                          invoking  update_admin_table() , so that the
 *                          physical parameters of the new node are availalbe
 *                          to determine stability of the new node.
 *-----------------------------------------------------------------------------
 */
local void  install_node(tuple, nodes, nnodeptr, admin_table, struct_table,
			 dist_table, radius, k)
int  tuple[NBODY_MAX];
bodyptr  nodes;
int  *nnodeptr;
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
real  dist_table[NNODE_MAX][NNODE_MAX];
real  radius;
int  k;
    {
    real  a;                                   /* semimajor axis of a binary */
    real  e;                                   /* eccentricity of a binary   */

    if (*nnodeptr >= NNODE_MAX)
	error("install_node: *nnodeptr = %d >= NNODE_MAX = %d\n",
	      *nnodeptr, NNODE_MAX);

    if (k == 2)
	get_binary_parameters(nodes + tuple[0], nodes + tuple[1], &a, &e);

    install_com_dynamics(nodes, *nnodeptr, tuple, k);

    update_dist_table(nodes, *nnodeptr, dist_table);
    update_struct_table(*nnodeptr, struct_table, radius, a, e, k);
    update_admin_table(tuple, nodes, *nnodeptr, admin_table, struct_table,
		       dist_table, k);
    *nnodeptr += 1;
    }

/*-----------------------------------------------------------------------------
 *  update_dist_table  --  
 *-----------------------------------------------------------------------------
 */
local void  update_dist_table(nodes, nnode, dist_table)
bodyptr  nodes;
int  nnode;
real  dist_table[NNODE_MAX][NNODE_MAX];
    {
    int  i;
    real  distance;

    dist_table[nnode][nnode] = 0.0;

    for (i = 0; i < nnode; i++)
	{
	DISTV(distance, Pos(nodes + i), Pos(nodes + nnode));
	dist_table[i][nnode] = dist_table[nnode][i] = distance;
	}
    }

#define  STABLE_SEPARATION_FACTOR     2.0              /* somewhat arbitrary */

/*-----------------------------------------------------------------------------
 *  update_admin_table  --  determines the values for the visibility and
 *                          stability of a new node, and enters those values
 *                          in the  admin_table[] .
 *                          Check the node for being barely or strongly 
 *                          unbound; if not, then bound, so check for
 *                          stable or unstable.
 *                     note:
 *                          the use of how_bound_tuple() may well be overkill,
 *                          but as long as it works, no problem for now.
 *    note: oldversion:
 *                     note:
 *                          a k-tuple with  k > 2  is considered to be always
 *                          unstable; a binary (k = 2) is stable if both
 *                          members are internally stable and in addition
 *                          the pericenter exceeds the sum of the radii of the
 *                          members by a factor  STABLE_SEPARATION_FACTOR .
 *-----------------------------------------------------------------------------
 */
local void  update_admin_table(tuple, nodes, nnode, admin_table, struct_table,
			       dist_table, k)
int  tuple[NBODY_MAX];
bodyptr  nodes;
int  nnode;
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
real  dist_table[NNODE_MAX][NNODE_MAX];
int  k;
    {
    int  i;
    bool  isbound;
    bool  isbarelyunbound;
    real  pericenter_distance;
    real  sum_of_radii;

    Ndaughter(admin_table, nnode) = k;
    Visibility(admin_table, nnode) = VISIBLE;

    for (i = 0; i < k; i++)
	{
        Daughter(admin_table, nnode, i) = tuple[i];
	Visibility(admin_table, tuple[i]) = COVERED;
	}

    how_bound_tuple(tuple, nodes, dist_table, k, &isbound, &isbarelyunbound);

    if (isbarelyunbound == TRUE)
	Stability(admin_table, nnode) = BARELY_UNBOUND;
    else if (isbound == FALSE)
	Stability(admin_table, nnode) = STRONGLY_UNBOUND;
    else if (k > 2)
	Stability(admin_table, nnode) = UNSTABLE;
    else
	{
	Stability(admin_table, nnode) = STABLE;
	for (i = 0; i < k; i++)
	    if (Ndaughter(admin_table, Daughter(admin_table, nnode, i)) > 0)
	        if (Stability(admin_table, Daughter(admin_table, nnode, i))
		    != STABLE)
		    Stability(admin_table, nnode) = UNSTABLE;

	pericenter_distance =
	    Apair(struct_table, nnode) * (1.0 - Epair(struct_table, nnode));
	sum_of_radii = Radius(struct_table, Daughter(admin_table, nnode, 0))
	             + Radius(struct_table, Daughter(admin_table, nnode, 1));

	if (pericenter_distance < STABLE_SEPARATION_FACTOR * sum_of_radii)
	    Stability(admin_table, nnode) = UNSTABLE;
	}
    }

/*-----------------------------------------------------------------------------
 *  update_struct_table  --  
 *-----------------------------------------------------------------------------
 */
local void  update_struct_table(nnode, struct_table, radius, a, e, k)
int  nnode;
real  struct_table[NNODE_MAX][N_STRUCT];
real  radius;
real  a, e;
int  k;
    {
    int  i;

    Radius(struct_table, nnode) = radius;

    if (k == 2)
	{
	Apair(struct_table, nnode) = a;
	Epair(struct_table, nnode) = e;
	}
    }

/*-----------------------------------------------------------------------------
 *  get_binary_parameters  --  
 *-----------------------------------------------------------------------------
 */
get_binary_parameters(body1, body2, aptr, eptr)
bodyptr  body1;
bodyptr  body2;
realptr  aptr;               /* pointer to a, the semimajor axis of a binary */
realptr  eptr;               /* pointer to e, the eccentricity of a binary   */
    {
    real  delta_r, delta_v, m_sum;
    real  r_rel[NDIM], v_rel[NDIM], r_out_v[NDIM];
    real  r_out_v_squared;

    DISTV(delta_r, Pos(body1), Pos(body2));
    DISTV(delta_v, Vel(body1), Vel(body2));
    m_sum = Mass(body1) + Mass(body2);

    *aptr = 1.0 / (2.0/delta_r - delta_v*delta_v / m_sum);

    SUBV(r_rel, Pos(body1), Pos(body2));
    SUBV(v_rel, Vel(body1), Vel(body2));
    r_out_v[0] = r_rel[1] * v_rel[2] - r_rel[2] * v_rel[1];
    r_out_v[1] = r_rel[2] * v_rel[0] - r_rel[0] * v_rel[2];
    r_out_v[2] = r_rel[0] * v_rel[1] - r_rel[1] * v_rel[0];
    DOTVP(r_out_v_squared, r_out_v, r_out_v);

    *eptr = sqrt(1.0 - (r_out_v_squared / (m_sum * *aptr)));
    }

/*-----------------------------------------------------------------------------
 *  iterate_dissolve_visible_strongly_unbound_nodes  --  
 *-----------------------------------------------------------------------------
 */
local void  iterate_dissolve_visible_strongly_unbound_nodes(nodes, nnodeptr,
							    admin_table,
							    struct_table,
							    dist_table)
bodyptr  nodes;
int *nnodeptr;
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
real  dist_table[NNODE_MAX][NNODE_MAX];
    {
    bool  dissolve_visible_strongly_unbound_nodes();

    while (dissolve_visible_strongly_unbound_nodes(nodes, nnodeptr,admin_table,
						   struct_table, dist_table)
	   == TRUE)
	;
    }	  

/*-----------------------------------------------------------------------------
 *  dissolve_visible_strongly_unbound_nodes  --  
 *-----------------------------------------------------------------------------
 */
local bool  dissolve_visible_strongly_unbound_nodes(nodes, nnodeptr,
						    admin_table, struct_table,
						    dist_table)
bodyptr  nodes;
int *nnodeptr;
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
real  dist_table[NNODE_MAX][NNODE_MAX];
    {
    int  i;
    bool  new_activity;

    new_activity = FALSE;
    for (i = 0; i < *nnodeptr; i++)
	if (Visibility(admin_table, i) == VISIBLE &&
	    Stability(admin_table, i) == STRONGLY_UNBOUND)
	    {
	    dissolve_node(nodes, nnodeptr, i, admin_table, struct_table,
			  dist_table);
	    new_activity = TRUE;
	    break;                        /* only one dissolution at a time, */
	    }                             /* to retain integrity of nnodeptr */

    return(new_activity);
    }

/*-----------------------------------------------------------------------------
 *  dissolve_node  --  
 *-----------------------------------------------------------------------------
 */
local void  dissolve_node(nodes, nnodeptr, i, admin_table, struct_table,
			  dist_table)
bodyptr  nodes;
int *nnodeptr;
int  i;
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
real  dist_table[NNODE_MAX][NNODE_MAX];
    {
    if (i != *nnodeptr-1)
	swap_nodes(nodes, *nnodeptr, i, *nnodeptr-1, admin_table, struct_table,
		   dist_table);

    drop_tailnode(nnodeptr, admin_table);
    }

/*-----------------------------------------------------------------------------
 *  swap_nodes  --  
 *-----------------------------------------------------------------------------
 */
local void  swap_nodes(nodes, nnode, i1, i2, admin_table, struct_table,
		       dist_table)
bodyptr  nodes;
int  nnode;
int  i1;
int  i2;
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
real  dist_table[NNODE_MAX][NNODE_MAX];
    {
    int  j;
    int  i;
    int  tmp_int;
    real  tmp_real;
    body  tmp_body_allocation;
    bodyptr  tmp_body;

    tmp_body = &tmp_body_allocation; /* easier to think in terms of pointers */

    if (i1 >= nnode)
	error("swap_nodes: i1 = %d > nnode = %d\n", i1, nnode);
    if (i2 >= nnode)
	error("swap_nodes: i2 = %d > nnode = %d\n", i2, nnode);

    cp_dynamics(tmp_body, nodes + i1);
    cp_dynamics(nodes + i1, nodes + i2);
    cp_dynamics(nodes + i2, tmp_body);

    tmp_int = Ndaughter(admin_table, i1);
    Ndaughter(admin_table, i1) = Ndaughter(admin_table, i2);
    Ndaughter(admin_table, i2) = tmp_int;

    tmp_int = Visibility(admin_table, i1);
    Visibility(admin_table, i1) = Visibility(admin_table, i2);
    Visibility(admin_table, i2) = tmp_int;

    tmp_int = Stability(admin_table, i1);
    Stability(admin_table, i1) = Stability(admin_table, i2);
    Stability(admin_table, i2) = tmp_int;
    
    j = MAX(Ndaughter(admin_table, i1), Ndaughter(admin_table, i2));
    while (j-- > 0)           /* this swaps also unused numbers, if the two  */
	{                     /* Ndaughter values are unequal, but who cares */
	tmp_int = Daughter(admin_table, i1, j);
	Daughter(admin_table, i1, j) = Daughter(admin_table, i2, j);
	Daughter(admin_table, i2, j) = tmp_int;
	}
/*
 * Now switch also the pointers to i1 and i2 from the outside,
 * including the cases where i = i1 and i = i2 !
 */
    for (i = 0; i < nnode; i++)
	if (Ndaughter(admin_table, i) > 0)
	    for (j = 0; j < Ndaughter(admin_table, i); j++)
		{
		if (Daughter(admin_table, i, j) == i1)
		    Daughter(admin_table, i, j) == i2;
		else if (Daughter(admin_table, i, j) == i2)
		    Daughter(admin_table, i, j) == i1;
		}

    tmp_real = Radius(struct_table, i1);
    Radius(struct_table, i1) = Radius(struct_table, i2);
    Radius(struct_table, i2) = tmp_real;
    
    tmp_real = Apair(struct_table, i1);
    Apair(struct_table, i1) = Apair(struct_table, i2);
    Apair(struct_table, i2) = tmp_real;
    
    tmp_real = Epair(struct_table, i1);
    Epair(struct_table, i1) = Epair(struct_table, i2);
    Epair(struct_table, i2) = tmp_real;
    
    for (i = 0; i < nnode; i++)
	{
	tmp_real = dist_table[i][i1];
	dist_table[i][i1] = dist_table[i][i2];
	dist_table[i][i2] = tmp_real;
	}

    for (j = 0; j < nnode; j++)
	{
	tmp_real = dist_table[i1][j];
	dist_table[i1][j] = dist_table[i2][j];
	dist_table[i2][j] = tmp_real;
	}
    }

/*-----------------------------------------------------------------------------
 *  drop_tailnode  --  
 *-----------------------------------------------------------------------------
 */
local void  drop_tailnode(nnodeptr, admin_table)
int *nnodeptr;
int  admin_table[NNODE_MAX][N_ADMIN];
    {
    int  j;

    if (Visibility(admin_table, *nnodeptr-1) != VISIBLE)
	error("drop_tailnode: can't drop an invisible tail\n");

    for (j = 0; j < Ndaughter(admin_table, *nnodeptr-1); j++)
	Visibility(admin_table, Daughter(admin_table, *nnodeptr-1, j))
	    = VISIBLE;

    (*nnodeptr)--;
    }

/*-----------------------------------------------------------------------------
 *  install_com_dynamics  --  compute and install the values of the  mass,
 *                            position and velocity for a new node, with the
 *                            index  nnode , and having  k  daughters,
 *                            whose indices are contained in the first  k
 *                            entrees of  tuple[] .
 *-----------------------------------------------------------------------------
 */
local void  install_com_dynamics(nodes, nnode, tuple, k)
bodyptr  nodes;
int  nnode;
int  tuple[NBODY_MAX];
int  k;
    {
    int  i;
    real  total_mass;
    real  com_pos[NDIM];
    real  com_vel[NDIM];
    real  mass_times_pos[NDIM];
    real  mass_times_vel[NDIM];
    bodyptr  member;
    bodyptr  new_node;

    total_mass = 0.0;
    CLRV(com_pos);
    CLRV(com_vel);
    for (i = 0; i < k; i++)
	{
        member = nodes + tuple[i];
	total_mass += Mass(member);
	MULVS(mass_times_pos, Pos(member), Mass(member));
	INCADDV(com_pos, mass_times_pos);
	MULVS(mass_times_vel, Vel(member), Mass(member));
	INCADDV(com_vel, mass_times_vel);
	}
    INCDIVVS(com_pos, total_mass);               /* the real c.o.m. position */
    INCDIVVS(com_vel, total_mass);               /* the real c.o.m. velocity */

    new_node = nodes + nnode;
    Mass(new_node) = total_mass;
    SETV(Pos(new_node), com_pos);
    SETV(Vel(new_node), com_vel);
    }

/*-----------------------------------------------------------------------------
 *  print_hierarchical_structure  --  
 *     for now I have removed the following extra if statement which would 
 *     limit the output verbiage:
 * 	    if ((Visibility(admin_table, i) == VISIBLE) || (debug == TRUE))
 *-----------------------------------------------------------------------------
 */
local void  print_hierarchical_structure(nnode, admin_table, struct_table,
					 dist_table, tnow)
int  nnode;
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
real  dist_table[NNODE_MAX][NNODE_MAX];
real  tnow;
    {
    int  i, j;
    int  blankfield_counter;
    bool  printflag;
    char  node_report[BUFF_LENGTH];
    real  distance;

    printf("\n");
    if (debug == TRUE)
	printf("  t = %8lf\n", tnow);

    for (i = 0; i < nnode; i++)
	if (Ndaughter(admin_table, i) > 0)
	    {
	    if (debug == TRUE)
		sprintf(node_report, "#%2d = ", i);
	    else
	        sprintf(node_report, "  ");

	    print_node(node_report, admin_table, struct_table, i);

	    if ((Ndaughter(admin_table, i) == 2))
		print_binary_parameters(node_report, admin_table, struct_table,
					i);
	    else
		print_radius(node_report, admin_table, struct_table, i);

	    printf(node_report);
	    }

    if (debug == TRUE)
	for (i = 0; i < nnode - 1; i++)
	    if (Visibility(admin_table, i) == VISIBLE)
		{
		printflag = FALSE;
		blankfield_counter = 0;
		if (i > 0)
		    for (j = 0; j < i; j++)
			if (Visibility(admin_table, j) == VISIBLE)
			    blankfield_counter++;
		for (j = i+1; j < nnode; j++)
		    if (Visibility(admin_table, j) == VISIBLE)
			{
			while (blankfield_counter-- >0)
			    printf("                   ");
			if (printflag == FALSE)
			    printf("  ");
			printflag = TRUE;

			if (i >= 10 && j >= 10)
			    printf("  r[%d,%d]=", i, j);
			else if (i >= 10 || j >= 10)
			    printf("  r[%d,%d]= ", i, j);
			else
			    printf("  r[%d,%d] = ", i, j);

			if ((distance = dist_table[i][j]) < 10.0)
			    printf("%.6lf", distance);
			else if (distance < 100.0)
			    printf("%.5lf", distance);
			else if (distance < 1000.0)
			    printf("%.4lf", distance);
			else if (distance < 10000.0)
			    printf("%.3lf", distance);
			else if (distance < 100000.0)
			    printf("%.2lf", distance);
			else if (distance < 1000000.0)
			    printf("%.1lf", distance);
			else
			    printf("%.0lf", distance);
			}
		if (printflag == TRUE)
		    printf("\n");
		}
    }

/*-----------------------------------------------------------------------------
 *  print_node  --  
 *-----------------------------------------------------------------------------
 */
local void  print_node(node_report, admin_table, struct_table, member)
char  node_report[BUFF_LENGTH];
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
int  member;
    {
    int  i;
    char *end_of_string();

    if (Ndaughter(admin_table, member) == 0)
	sprintf(end_of_string(node_report), "%d", member);
    else
	{
	if (Stability(admin_table, member) == STABLE)
	    sprintf(end_of_string(node_report), "[");
	else if (Stability(admin_table, member) == UNSTABLE)
	    sprintf(end_of_string(node_report), "(");
	else if (Stability(admin_table, member) == BARELY_UNBOUND)
	    sprintf(end_of_string(node_report), "<");
	else if (Stability(admin_table, member) == STRONGLY_UNBOUND)
	    sprintf(end_of_string(node_report), "{");
	else if (Stability(admin_table, member) == UNKNOWN)
	    error("print_node: stability of member %d is UNKNOWN\n", member);
	else
	    error("print_node: invalid value %d for stability of member %d\n",
		  Stability(admin_table, member), member);
	
	for (i = 0; i < Ndaughter(admin_table, member); i++)
	    {
	    if (i > 0)
		sprintf(end_of_string(node_report), ", ");
	    print_node(node_report, admin_table, struct_table,
	    	       Daughter(admin_table, member, i));
	    }

	if (Stability(admin_table, member) == STABLE)
	    sprintf(end_of_string(node_report), "]");
	else if (Stability(admin_table, member) == UNSTABLE)
	    sprintf(end_of_string(node_report), ")");
	else if (Stability(admin_table, member) == BARELY_UNBOUND)
	    sprintf(end_of_string(node_report), ">");
	else if (Stability(admin_table, member) == STRONGLY_UNBOUND)
	    sprintf(end_of_string(node_report), "}");
	else if (Stability(admin_table, member) == UNKNOWN)
	    error("print_node: stability of member %d is UNKNOWN\n", member);
	else
	    error("print_node: invalid value %d for stability of member %d\n",
		  Stability(admin_table, member), member);
	}
    }

/*-----------------------------------------------------------------------------
 *  print_binary_parameters  --  
 *-----------------------------------------------------------------------------
 */
local void  print_binary_parameters(node_report, admin_table, struct_table,
				    member)
char  node_report[BUFF_LENGTH];
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
int  member;
    {
    int  offset_in_node_report;
    char *end_of_string();
    real  semimajor_axis;

    offset_in_node_report = end_of_string(node_report) - node_report;
/*
 * if the radius information starts in the 49th column, the output will be
 * nicely lined up with the distances printed under the debug option.
 */
    while (offset_in_node_report < 49)
	node_report[offset_in_node_report++] = ' ';

    semimajor_axis = Apair(struct_table, member);
    if (semimajor_axis < 10.0 && semimajor_axis > 0.0)
	sprintf(node_report + offset_in_node_report,"a = %.6lf  ;  e = %8lf\n",
		semimajor_axis, Epair(struct_table, member));
    else if (semimajor_axis < 100.0 && semimajor_axis > -10.0)
	sprintf(node_report + offset_in_node_report,"a = %.5lf  ;  e = %8lf\n",
		semimajor_axis, Epair(struct_table, member));
    else if (semimajor_axis < 1000.0 && semimajor_axis > -100.0)
	sprintf(node_report + offset_in_node_report,"a = %.4lf  ;  e = %8lf\n",
		semimajor_axis, Epair(struct_table, member));
    else if (semimajor_axis < 10000.0 && semimajor_axis > -1000.0)
	sprintf(node_report + offset_in_node_report,"a = %.3lf  ;  e = %8lf\n",
		semimajor_axis, Epair(struct_table, member));
    else if (semimajor_axis < 100000.0 && semimajor_axis > -10000.0)
	sprintf(node_report + offset_in_node_report,"a = %.2lf  ;  e = %8lf\n",
		semimajor_axis, Epair(struct_table, member));
    else if (semimajor_axis < 1000000.0 && semimajor_axis > -100000.0)
	sprintf(node_report + offset_in_node_report,"a = %.1lf  ;  e = %8lf\n",
		semimajor_axis, Epair(struct_table, member));
    else
	sprintf(node_report + offset_in_node_report,"a = %.0lf  ;  e = %8lf\n",
		semimajor_axis, Epair(struct_table, member));
    }

/*-----------------------------------------------------------------------------
 *  print_radius  --  
 * NOTE: CLEANUP THE SILLY OUTPUT SWITCHES WITH A NEW FIXED_WIDTH %lf PROCEDURE
 *-----------------------------------------------------------------------------
 */
local void  print_radius(node_report, admin_table, struct_table, member)
char  node_report[BUFF_LENGTH];
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
int  member;
    {
    int  offset_in_node_report;
    char *end_of_string();
    real  r_node;

    offset_in_node_report = end_of_string(node_report) - node_report;
/*
 * if the radius information starts in the 49th column, the output will be
 * nicely lined up with the distances printed under the debug option.
 */
    while (offset_in_node_report < 49)
	node_report[offset_in_node_report++] = ' ';

    r_node = Radius(struct_table, member);
    if (r_node < 10.0)
	sprintf(node_report + offset_in_node_report, "R = %.6lf\n", r_node);
    else if (r_node < 100.0)
	sprintf(node_report + offset_in_node_report, "R = %.5lf\n", r_node);
    else if (r_node < 1000.0)
	sprintf(node_report + offset_in_node_report, "R = %.4lf\n", r_node);
    else if (r_node < 10000.0)
	sprintf(node_report + offset_in_node_report, "R = %.3lf\n", r_node);
    else if (r_node < 100000.0)
	sprintf(node_report + offset_in_node_report, "R = %.2lf\n", r_node);
    else if (r_node < 1000000.0)
	sprintf(node_report + offset_in_node_report, "R = %.1lf\n", r_node);
    else
	sprintf(node_report + offset_in_node_report, "R = %.0lf\n", r_node);
    }

/*-----------------------------------------------------------------------------
 *  how_much_offspring  --  computes the number of real particles below the
 *                          the node  member .
 *-----------------------------------------------------------------------------
 */
local int  how_much_offspring(admin_table, member)
int  admin_table[NNODE_MAX][N_ADMIN];
int  member;
    {
    int  i, k, n_total;

    k = Ndaughter(admin_table, member);
    n_total = 0;
    
    if (k > 0)
	for (i = 0; i < k; i++)
	    n_total += how_much_offspring(admin_table,
					  Daughter(admin_table, member, i));
    else
	n_total = 1;

    return(n_total);
    }

/******  #define  BUFF_LENGTH    128  *******  this is now in file  diag.h   */

/*-----------------------------------------------------------------------------
 *  unique_decomposition  --  
 *-----------------------------------------------------------------------------
 */
local void  unique_decomposition(diags, nnode, admin_table, struct_table)
diagptr  diags;
int  nnode;
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
    {
    int  i;
    char  hier_string[BUFF_LENGTH];

    hier_string[0] = NULL;               /* initialize with a null_string,  */
                                         /* for  end_of_string()  in        */
                                         /* hier_string()  to work properly */

    for (i = 0; i < nnode; i++)
	if (Visibility(admin_table, i) == VISIBLE)
	    sprint_node(hier_string, admin_table, struct_table,	i);

    sprint_node_end(hier_string);
/*
 *  printf("ARBITRARY DECOMPOSITION:  %s\n", hier_string);
 */
    bring_to_normal_form(hier_string);

    printf("\nnormal form:  %s\n", hier_string);

    sprintf(Hier_string(diags), "%s", hier_string);
    }

/*-----------------------------------------------------------------------------
 *  sprint_node  --  
 *-----------------------------------------------------------------------------
 */
local void  sprint_node(hier_string, admin_table, struct_table,	member)
string  hier_string;
int  admin_table[NNODE_MAX][N_ADMIN];
real  struct_table[NNODE_MAX][N_STRUCT];
int  member;
    {
    int  i;
    char *end_of_string();

    insert_delimiter_if_needed(hier_string);

    if (Ndaughter(admin_table, member) == 0)
	sprintf(end_of_string(hier_string), "%d", member);
    else
	{
	if (Stability(admin_table, member) == STABLE)
	    sprintf(end_of_string(hier_string), "[");
	else if (Stability(admin_table, member) == UNSTABLE)
	    sprintf(end_of_string(hier_string), "(");
	else if (Stability(admin_table, member) == BARELY_UNBOUND)
	    sprintf(end_of_string(hier_string), "<");
	else if (Stability(admin_table, member) == STRONGLY_UNBOUND)
	    sprintf(end_of_string(hier_string), "{");
	else if (Stability(admin_table, member) == UNKNOWN)
	    error("sprint_node: stability of member %d is UNKNOWN\n", member);
	else
	    error("sprint_node: invalid value %d for stability of member %d\n",
		  Stability(admin_table, member), member);
	
	for (i = 0; i < Ndaughter(admin_table, member); i++)
	    sprint_node(hier_string, admin_table, struct_table,
		        Daughter(admin_table, member, i));

	if (Stability(admin_table, member) == STABLE)
	    sprintf(end_of_string(hier_string), "]");
	else if (Stability(admin_table, member) == UNSTABLE)
	    sprintf(end_of_string(hier_string), ")");
	else if (Stability(admin_table, member) == BARELY_UNBOUND)
	    sprintf(end_of_string(hier_string), ">");
	else if (Stability(admin_table, member) == STRONGLY_UNBOUND)
	    sprintf(end_of_string(hier_string), "}");
	else if (Stability(admin_table, member) == UNKNOWN)
	    error("sprint_node: stability of member %d is UNKNOWN\n", member);
	else
	    error("sprint_node: invalid value %d for stability of member %d\n",
		  Stability(admin_table, member), member);
	}
    }

/*-----------------------------------------------------------------------------
 *  sprint_node_end  --  appends the terminator symbol  ;  to the hierarchical
 *                       decomposition string  hier_string .
 *-----------------------------------------------------------------------------
 */
local void  sprint_node_end(hier_string)
string  hier_string;
    {
    char *end_of_string();

    sprintf(end_of_string(hier_string), ";");
    }

/*-----------------------------------------------------------------------------
 *  insert_delimiter_if_needed  --  appends the delimiter symbol  ,  if needed,
 *                                  i.e. when the last character of the 
 *                                  hierarchical decomposition string
 *                                   hier_string  is not an indicator of the
 *                                  beginning of the string (a NULL) or the
 *                                  beginning of the substring (a '[' or '('
 *                                  or '{' or '<').
 *-----------------------------------------------------------------------------
 */
local void  insert_delimiter_if_needed(hier_string)
string  hier_string;
    {
    int  char_position;
    char  last_c;

    if (hier_string[0] == NULL)    /* beginning of an, as yet empty, string: */
	return;                    /*   no delimiter needed.                 */

    char_position = 1;

    while(hier_string[char_position] != NULL)
	char_position++;

    if(char_position > BUFF_LENGTH - 2)
	error("insert_delimiter_if_needed: too little room left in buffer\n");

    last_c = hier_string[char_position - 1];
    if (last_c == '(' || last_c == '[' || last_c == '{' || last_c == '<')
	return;            /* beginning of substring:  no delimiter needed.  */
    else
	sprintf(hier_string + char_position, ",");
    }


#define  BUFF_SAFETY_MARGIN  28

/*-----------------------------------------------------------------------------
 *  end_of_string  --  
 *-----------------------------------------------------------------------------
 */
local char *end_of_string(the_string)
char *the_string;
    {
    int  char_position;

    char_position = 0;
    while(the_string[char_position] != NULL)
	char_position++;

    if(char_position > BUFF_LENGTH - 1 - BUFF_SAFETY_MARGIN)
	error("end_of_string: too little room left in buffer\n");

    return(the_string + char_position);
    }

/*-----------------------------------------------------------------------------
 *  bring_to_normal_form  --  
 *                   example: 3,(7,[5,[4,8]],0),6,(2,1);  ==>
 *                            (0,[[4,8],5],7),(1,2),3,6;
 *                      NOTE: only one-digit star numbering implemented, i.e.
 *                            only valid for an N-body sub-system with N < 11
 *                  OH, WELL: a much cleaner way, of course, would be to first
 *                            map the string to an integer array of symbols,
 *                            where, e.g.,  ";" --> -1 , "(" --> -2 , 
 *                            ")" --> -3 , "[" --> -4 , "]" --> -5 , etc,
 *                            and the delimiting , is simply skipped since the
 *                            array elements are automatically separated.
 *-----------------------------------------------------------------------------
 */
local void  bring_to_normal_form(old_hier_string)
char  old_hier_string[BUFF_LENGTH];
    {
    int  i, j;
    int  number_of_objects;
    int  ordered[BUFF_LENGTH];
    int  unordered[BUFF_LENGTH];
    char  new_hier_string[BUFF_LENGTH];
    char  substring[BUFF_LENGTH];
    char  head_char, tail_char;

    if (old_hier_string[1] == ';')
	return;                     /* a singleton is already in normal form */

    new_hier_string[0] = NULL;         /* initialize with a null_string, for */
    substring[0] = NULL;              /* end_of_string()  to work properly  */

    map_to_integers(old_hier_string, unordered, &number_of_objects);

    sort_int_array(unordered, ordered, number_of_objects);

    for (i = 0; i < number_of_objects; i++)
	{
	j = old_ranking(unordered, ordered, i, number_of_objects);
        select_object(substring, old_hier_string, j);
	unwrap_string(substring, &head_char, &tail_char);
	bring_to_normal_form(substring);
	wrap_string(substring, head_char, tail_char);
	append_object(new_hier_string, substring);
	}
/*
 * copy the new string back unto the old string:
 */
    i = 0;    
    while (new_hier_string[i] != NULL && i < BUFF_LENGTH)
	{
	old_hier_string[i] = new_hier_string[i];
	i++;
	}
    if (i >= BUFF_LENGTH)
	error("bring_to_normal_form: buffer length exceeded\n");
    if (old_hier_string[i] != NULL)
	error("bring_to_normal_form: new and old hier. string diff. length\n");
    }

/*-----------------------------------------------------------------------------
 *  map_to_integers  --  
 *                 NOTE: only one-digit star numbering implemented, i.e.
 *                       only valid for an N-body sub-system with N < 11
 *-----------------------------------------------------------------------------
 */
local void  map_to_integers(hier_string, unordered, nptr)
char  hier_string[BUFF_LENGTH];
int  unordered[BUFF_LENGTH];
int *nptr;                             /* will count total number of objects */
    {
    int  i;
    int  lowest_number;
    int  level;
    char  c;

    *nptr = 0;
    level = 0;
    lowest_number = 10;        /* higher than highest allowed number */
    i = 0;
    while (i < BUFF_LENGTH)
	{
	c = hier_string[i++];

	if (c >= '0' && c <= '9')
	    {
	    if (c - '0' < lowest_number)
		lowest_number = c - '0';
	    }
	else if (c == '(' || c == '[' || c == '{' || c == '<')
	    level++;
	else if (c == ')' || c == ']' || c == '}' || c == '>')
	    level--;
	else if (level == 0 && (c == ',' || c == ';'))
	    {
	    unordered[*nptr] = lowest_number;
            (*nptr)++;
	    lowest_number = 10;
	    }

	if (c == ';')
	    break;
	}
    if (i >= BUFF_LENGTH)
	error("map_to_integers: buffer length exceeded\n");
    }

/*-----------------------------------------------------------------------------
 *  sort_int_array  --  
 *                NOTE: only one-digit star numbering implemented, i.e.
 *                      only valid for an N-body sub-system with N < 11
 *-----------------------------------------------------------------------------
 */
local void  sort_int_array(old_array, new_array, n)
int  old_array[BUFF_LENGTH];
int  new_array[BUFF_LENGTH];
int  n;
    {
    int  i, j;
    int  dummy;

    for (i = 0; i < n; i++)
	new_array[i] = old_array[i];

    for (i = 0; i < n - 1; i++)         /* lazy programmer's sorting method  */
	for (j = i+1; j < n; j++)
	    if (new_array[j] < new_array[i])
		{
		dummy = new_array[i];
		new_array[i] = new_array[j];
		new_array[j] = dummy;
		}
    }

/*-----------------------------------------------------------------------------
 *  old_ranking  --  
 *               NOTE: only one-digit star numbering implemented, i.e.
 *                     only valid for an N-body sub-system with N < 11
 *-----------------------------------------------------------------------------
 */
local void  old_ranking(old_array, new_array, i, n)
int  old_array[BUFF_LENGTH];
int  new_array[BUFF_LENGTH];
int  i;
int  n;
    {
    int  j;

    j = 0;
    while (j < n)
	if (old_array[j++] == new_array[i])
	    break;
    if (j > n)
	error("old_ranking: used buffer length exceeded\n");

    return(--j);
    }

/*-----------------------------------------------------------------------------
 *  select_object  --  
 *               NOTE: only one-digit star numbering implemented, i.e.
 *                     only valid for an N-body sub-system with N < 11
 *-----------------------------------------------------------------------------
 */
local void  select_object(tmp_string, hier_string, member)
char  tmp_string[BUFF_LENGTH];
char  hier_string[BUFF_LENGTH];
int  member;
    {
    int  i, j;
    int  level;
    char  c;

    level = 0;
    i = 0;
    while (member > 0)
	{
	c = hier_string[i];
	if (c == '(' || c == '[' || c == '{' || c == '<')
	    level++;
	else if (c == ')' || c == ']' || c == '}' || c == '>')
	    level--;
	else if (c == ',' && level == 0)
	    member--;
	else if (c == ';')
	    error("select_object: reached end of string\n");
        if (++i >= BUFF_LENGTH)
	    error("select_object: at position #1: buffer length exceeded\n");
	}

    if (level != 0)
	error("select_object: number of (,[,{,< and of ),],},> unbalanced\n");
    j = 0;
    while (c = hier_string[i])
	{
	if (c == '(' || c == '[' || c == '{' || c == '<')
	    level++;
	else if (c == ')' || c == ']' || c == '}' || c == '>')
	    level--;
	else if ((c == ',' && level == 0) || c == ';')
	    {
	    sprintf(tmp_string + j, ";");      /* end with ';' and  NULL  !! */
	    if (level != 0)
		error("select_object: number of (,[, etc. unbalanced\n");
	    break;
	    }
	tmp_string[j++] = c;

        if (++i >= BUFF_LENGTH)
	    error("select_object: at position #2: buffer length exceeded\n");
	}
    }

/*-----------------------------------------------------------------------------
 *  append_object  --  
 *               NOTE: only one-digit star numbering implemented, i.e.
 *                     only valid for an N-body sub-system with N < 11
 *-----------------------------------------------------------------------------
 */
local void  append_object(hier_string, tmp_string)
char  hier_string[BUFF_LENGTH];
char  tmp_string[BUFF_LENGTH];
    {
    char *last_endptr;
    char *new_startingptr;
    char *end_of_string();

    if (hier_string[0] == NULL)
	new_startingptr = hier_string;
    else
	{
	last_endptr = end_of_string(hier_string) - 1;
	if (*last_endptr == ';')
	    *last_endptr = ',';                      /* insert new delimiter */
	else	    
	    error("append_object: hier_string not properly terminated\n");
	new_startingptr = last_endptr + 1;
	}

    sprintf(new_startingptr, tmp_string);
    }

/*-----------------------------------------------------------------------------
 *  unwrap_string  --  takes off the leading and trailing bracket,
 *                     if present, and only if the string 'substring' contains
 *                     only one object (the first two 'return;' statements
 *                     take care of these two exceptions).
 *                     If no enclosing brackets are present, the head and tail
 *                     characters '*head_ptr' and '*tail_ptr' are assigned the
 *                     value '*' to indicate the absence of brackets.
 *               NOTE: only one-digit star numbering implemented, i.e.
 *                     only valid for an N-body sub-system with N < 11
 *-----------------------------------------------------------------------------
 */

local void  unwrap_string(substring, head_ptr, tail_ptr)
char  substring[BUFF_LENGTH];
char *head_ptr;
char *tail_ptr;
    {
    int  i;
    int  level;
    int  n;
    char  c;

    c = substring[0];

    if (c == '(' || c == '[' || c == '{' || c == '<')
	*head_ptr = c;
    else
	{
	*head_ptr = *tail_ptr = '*';
	return;                         /* no enclosing brackets of any kind */
	}

    level = 1;
    i = 1;
    for (;;)
	{
	c = substring[i];
	if (c == '(' || c == '[' || c == '{' || c == '<')
	    level++;
	else if (c == ')' || c == ']' || c == '}' || c == '>')
	    level--;

	if (level == 0)
	    {
	    if (substring[i+1] != ';')    /* more than one object,           */
		{                         /* therefore no enclosing brackets */
		*head_ptr = *tail_ptr = '*';
		return;
		}
	    else
		break;
	    }

        if (++i >= BUFF_LENGTH)
	    error("unwrap_string: buffer length exceeded\n");
	}

    i = 1;
    while ((c = substring[i]) != ';')
	{
	substring[i-1] = c;
	i++;
	}
    *tail_ptr = substring[i-2];
    substring[i-2] = ';';
    substring[i-1] = NULL;                    /* proper string ending */
    }

/*-----------------------------------------------------------------------------
 *  wrap_string  --  puts the leading and trailing bracket in place again,
 *                   that is, if there are any brackets.
 *       convention: the character '*' for the value of 'head_char' and 
 *                   'tail_char' indicates the absence of brackets.
 *             NOTE: only one-digit star numbering implemented, i.e.
 *                   only valid for an N-body sub-system with N < 11
 *-----------------------------------------------------------------------------
 */
local void  wrap_string(substring, head_char, tail_char)
char  substring[BUFF_LENGTH];
char head_char;
char tail_char;
    {
    int  i;
    char  last_c;
    char  c;

    if (head_char == '*')
	return;

    last_c = head_char;
    i = 0;
    while((c = substring[i]) != ';')
	{
	substring[i] = last_c;
	last_c = c;
	i++;
	}
    substring[i] = last_c;
    substring[++i] = tail_char;
    substring[++i] = ';';
    substring[++i] = NULL;
    }

/* endof: scatdiagnose.c */
