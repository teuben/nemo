/* scatter.c - */

/*
 *  scatter.c:  orchestrates a simulation of one scattering event
 *
 *     July 1988  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
#include  "newton0.h"
#include  "scat.h"
/* #include  <getparam.h> */
#include  "getparam.h"

bool  debug;              /* flag for output control, for debugging purposes */
bool  turbo;              /* flag for integration spee-up, through shortcuts */

/*-----------------------------------------------------------------------------
 *  scatter: performs a scattering event with one target, located originally
 *           in the origin, and one projectile. The target as well as the
 *           projectile can consist of a single star or a binary; in the latter
 *           case each of the members of the binary themselves can be a single
 *           star or a binary, recursively. After setting up the whole system
 *           of particles, a coordinate transformation is performed to the
 *           center-of-mass frame of the system as a whole.
 *  notation:
 *           rho:   impact parameter
 *           psi:   impact angle
 *           r_in:  initial separation of single star from binary c.o.m.
 *           v_inf: relative velocity at infinity of single star vs. binary
 *           theta: spherical coordinate of incoming direction of single star
 *           phi:   see theta.
 *           
 *  UNITS: the gravitational constant is unity:                   G = 1
 *         semimajor axis of the target binary is unity:          a = 1
 *         total mass of the target binary is unity:          M + M = 1
 *                                                             1   2
 *  Space homogeneity: The center of mass of the target binary is initially at
 *                     rest in the origin.
 *  Space isotropy:    The target binary lies in the {x,y} plane, with the
 *                     periastron position of star #2 on the positive x-axis.
 *  Time homogeneity:  The zero of time is chosen to coincide with the 
 *                     pericenter passage of the internal binary of the target.
 *
 *  Implementation:  scatter.c           replaces      newton0.c
 *                   scatdiagnose.c      replaces      diagnose.c
 *
 *                  in the compilation of "newton0", resulting in "scatter".
 *-----------------------------------------------------------------------------
 */

string  defv[] =            /* DEFAULT INPUT PARAMETERS                      */
    {           
    "lastout=",             /* file name for one single last output          */
    "out=",		    /* output file name                              */
    "t_end=1.e9",           /* time at which the integration should halt     */
    "dt_minor=10.0",        /* time interval between minor outputs           */
    "dt_major=100.0",       /* time interval between major outputs           */
    "dt_max=1.e9",          /* maximum size of integration time step         */
    "nstep=1000000000",     /* maximum allowed number of integration steps   */
    "de_rel_max=0.1",       /* maximum allowed rel drift in total energy     */
    "nstep_de=1000000000",  /* number of steps between checking  de_rel_max  */
    "dn_major=1000000000",  /* number of integration steps between           */
                            /* major outputs                                 */
    "min_pairdist=0.0",     /* minimum allowed pair distance (softened)      */
    "cpu_max=2.0e9",        /* maximum allowed amount of CPU time (minutes)  */
    "cpu_last_call=1.0e9",  /* amount of CPU time after which integration    */
                            /* will be halted at next major-output time.     */
    "eta_acc=0.01",         /* dimensionless integration accuracy parameter, */
			    /* used in determining the size of the time step.*/
    "r0_soft=0.0",          /* potential softening length                    */
    "soft_focus=plummer",   /* type of softening used                        */
    "diagnostics=scatter",  /* set of diagnostics provided at output times   */
    "scheme=runge_kutta_4", /* integration scheme                            */
#ifndef REGULARIZATION
    "timestep=collision_time",  /* timestep criterion                        */
#else
    "timestep=constant",    /* timestep criterion                            */
    "niter=3",              /* number of iterations in arriving at the       */
                            /* correct output times                          */
#endif
    "prompt=TRUE",          /* for interactive parameter input from terminal */
    "debug=TRUE",           /* output information for debugging purposes     */
    "turbo=TRUE",           /* integration speed-up through shortcuts        */
    "headline=",            /* verbiage for output                           */
    NULL,
    };

/*-----------------------------------------------------------------------------
 *  main  --  read the command line arguments, and pass control to  scatter()
 *-----------------------------------------------------------------------------
 */
main(argc, argv)
int  argc;
string  argv[];
    {
    initparam(argv, defv);     		           /* setup parameter access */

    debug = getbparam("debug");
    turbo = getbparam("turbo");

    scatter();
    }

/*-----------------------------------------------------------------------------
 *  scatter  --  currently orchestrates only one gravitational evolution
 *-----------------------------------------------------------------------------
 */
local void  scatter()
    {
    stateptr  new_state;
    stateptr  setup_experiment();

    new_state = setup_experiment();

    perform_experiment(new_state);

    report_experiment(new_state);
    }

/*-----------------------------------------------------------------------------
 *  setup_experiment  --  set variables in diverse state substructures,
 *                        for a many-body scattering experiment.
 *                accept: new_state: a state pointer.
 *                  note: the fourth state substructure, containing the
 *                        diagnostics, cannot yet be set at this point, since
 *                        that would require energy diagnostics, which requires
 *                        revving the engines, which could not have been done
 *                        before the other three substructures are set.
 *                        Therefore, the diagnostics are currently set in the
 *                        file  diagnose.c .
 *-----------------------------------------------------------------------------
 */
local stateptr  setup_experiment()
    {
    stateptr  new_state;
    stateptr  mk_state();
    systptr  setup_scene();
    string local_announcement = "\tSetting up a scattering experiment";

    new_state = mk_state();       /* in  statealgebra.c . The new state does */
                                  /* not yet contain a System or Regsystem.  */
    System(new_state) = setup_scene();
    set_specs(Specs(new_state));
    set_controls(Ctrls(new_state), Tnow(System(new_state)));

    Announcement(Ctrls(new_state)) = local_announcement;

    return(new_state);
    }

/*-----------------------------------------------------------------------------
 *  perform_experiment  --  performs a scattering experiment
 *-----------------------------------------------------------------------------
 */
local void  perform_experiment(some_state)
stateptr  some_state;
    {
    rev_engine(some_state);                 /* determination of total energy;*/
                                            /* resides in file  orbit.c .    */

    write_initial_output(some_state, TRUE); /* resides in file  out.c .      */
                                           
    evolve(some_state);                     /* resides in file  orbit.c .    */
    }

/*-----------------------------------------------------------------------------
 *  report_experiment  --  reports the outcome of a scattering experiment
 *-----------------------------------------------------------------------------
 */
local void  report_experiment(old_state)
stateptr  old_state;
    {
    write_final_output(old_state, TRUE, TRUE);            /* in file  out.c */
    }

/*-----------------------------------------------------------------------------
 *  set_specs  --  set spec variables
 *                 accept: specs: a spec pointer.
 *-----------------------------------------------------------------------------
 */
#define  BIGNUMBER  1.0e20

local void  set_specs(specs)
specptr  specs;
    {
    Softfocus(specs) = getparam("soft_focus");
    Timestepmethod(specs) = getparam("timestep");
    Integrationscheme(specs) = getparam("scheme");
    Diagnostics(specs) = getparam("diagnostics");

    Stepparam(specs) = getdparam("eta_acc");
    Softparam(specs) = getdparam("r0_soft");
#ifdef EXTRAPOLATION
    if (Softparam(specs) > 1.0/BIGNUMBER)
	error("set_specs: scatter_ext: not yet non-zero softening length\n");
#endif
    Nstep_de(specs) = getiparam("nstep_de");
    SPECmin_pair(specs) = BIGNUMBER;

    if (turbo == TRUE)
	{
	if (!streq("plummer", Softfocus(specs)))
	    error("set_specs: turbo option requires  ->plummer<- n");
	else if (!streq("collision_time", Timestepmethod(specs)))
	    error("set_specs: turbo option requires  ->collision_time<- \n");
	else if (!streq("runge_kutta_4", Integrationscheme(specs)))
	    error("set_specs: turbo option requires ->runge_kutta_4<- \n");
	if (Softparam(specs) > 1.0/BIGNUMBER)
	    error("set_specs: turbo option not yet for non-zero softening\n");
	}
    }

/*-----------------------------------------------------------------------------
 *  set_controls  --  set control variables
 *                    accept:         ctr: a control pointer;
 *                           initial_time: needed to determine the boolean
 *                                         value of Forwards(ctr).
 *-----------------------------------------------------------------------------
 */
local void  set_controls(ctr, initial_time)
ctrlptr  ctr;
real  initial_time;
    {
    string  headline_helper;
#ifndef REGULARIZATION
#  ifndef EXTRAPOLATION
    string  headline_default ="Scatter code: shared time steps, turbo version";

    if (turbo == FALSE)                            /* override defined value */
	sprintf(headline_default, "Scatter code: shared time steps");
#  endif
#endif
#ifdef REGULARIZATION
    string  headline_default =
            "scatter code: shared time steps & regularization";
#endif
#ifdef EXTRAPOLATION
    string  headline_default = 
         "scatter code: individual time steps: polynomial orbit extrapolation";
#endif

/*
 * determine the temporal relation between begin and end point of integration:
 *   Note: the control variable "Forwards()" is an essential ingredient
 *         in the temporal testing macros "earlier" and "later" (in newton0.h)
 */
    Tend(ctr) = getdparam("t_end");
    Forwards(ctr) = (Tend(ctr) > initial_time);
    if (! Forwards(ctr))
	error("scatter: only scattering forwards in time\n");
/*
 * initialize the other control variables which can directly receive outside
 * information through the command line arguments:
 */
    Lastoutfile(ctr) = getparam("lastout");
    Outfile(ctr) = getparam("out");
    DTminor(ctr) = getdparam("dt_minor");
    DTmajor(ctr) = getdparam("dt_major");
    DTsave(ctr) = 0.0;                         /* no savings */
    DTmax(ctr) = getdparam("dt_max");
    Nmaxstep(ctr) = getiparam("nstep");
    DNmajor(ctr) = getdparam("dn_major");
#ifdef REGULARIZATION
    Ntimingiter(ctr) = getiparam("niter");
#endif
    DErelmax(ctr) = getdparam("de_rel_max");
    CPUmax(ctr) = getdparam("cpu_max");
    CPUlastcall(ctr) = getdparam("cpu_last_call");
    Min_pairdist(ctr) = getdparam("min_pairdist");

    headline_helper = getparam("headline");
    if (*headline_helper != NULL)                /* 3 possible headlines:    */
	Headline(ctr) = headline_helper;         /* 1st choice: command-line */
    else if (Headline(ctr) == NULL)              /* 2nd choice: input file   */
	Headline(ctr) = headline_default;        /* 3rd choice: default      */
/*
 * initialize the remaining control variables:
 */
    Tmajor(ctr) = Tminor(ctr) = Tsave(ctr) = initial_time;
    Nmajor(ctr) = 0;
    DTstep(ctr) = 0.0;                       /* not necessary, but good form */
    }

#define  NPBODY_MAX  128       /* maximum number of particles which can be   */
                               /* read in for a single scattering experiment */
#define  UNKNOWN   1                     /* these are the possible      */
#define  SINGLE    2                     /* contents of the elements    */
#define  BINARY    3                     /* pseudo_body_table[...][0]   */

#define  ROOT_ADDRESS          0    /* these are the offsets of the root,    */
#define  TARGET_ADDRESS        1    /* target and projectile in the array:   */
#define  PROJECTILE_ADDRESS    2    /* pseudo_body_table[some_address][...]  */

/*-----------------------------------------------------------------------------
 *  setup_scene  --  initializes the masses, time and orbital parameters
 *                   for a many-body scattering experiment.
 *                   flow:
 *                        first the orbit of the center of mass of a projectile
 *                        is read in, and the clock is set at a desired time
 *                        before closest encounter of the projectile and the
 *                        target.
 *                        note: i) first the time of projected encounter is
 *                                 read in (i.e. the time at which the 
 *                                 encounter would take place if the projectile
 *                                 would keep approaching the target in a 
 *                                 perfect two-body hyperbola). This time is
 *                                 taken with respect to the time t=0 defined
 *                                 by the time of pericenter passage of the 
 *                                 internal binary in the target.
 *                                 The motivation for this choice is to
 *                                 minimize the freedom in setting up the
 *                                 target binary.
 *                             ii) then the time is transformed to a different
 *                                 time frame in which the projected encounter
 *                                 occurs at time t=0.
 *                                 The motivation is to provide a more 
 *                                 convenient time frame for performing the
 *                                 scattering experiment.
 *                        then internal structure is added, first to the target
 *                        and then to the projectile; 
 *                   note: I: a SECOND choice has to be made to fix
 *                            the time scale, namely the time at which
 *                            the numerical orbit integration starts,
 *                            i.e. the time at which the initial
 *                            conditions should be given for all the
 *                            particles involved in the scattering
 *                            event. This has nothing to do with 
 *                            spacetime symmetries, but would occur
 *                            for every simulation with a less-then-
 *                            infinite duration. It has a purely
 *                            practical, and APPROXIMATE nature,
 *                            since it determines the point in time
 *                            before which the interactions between
 *                            target and projectiles are neglected.
 *                            As a consequence, chosing a different
 *                            starting time will lead to (slightly)
 *                            different results -- it is therefore up
 *                            to the experimenter to chose the initial
 *                            time judiciously.

 *                        II: The zero point of our time scale, defined
 *                            as the time of projected closest
 *                            encounter, differs slightly from the
 *                            actual time of closest encounter because
 *                            of the interactions between the internal
 *                            degrees of freedom of the target and 
 *                            projectiles -- besides, the actual time 
 *                            of closest encounter is likely to be
 *                            ill-defined when the internal structures
 *                            are taken into account. However, in the
 *                            limit of a computation starting
 *                            infinitely early, our choice of temporal
 *                            zero point is well defined. Of course,
 *                            as discussed under point I, truncation
 *                            of the computed part of the orbits to a
 *                            finite time interval introduces (slight)
 *                            perturbations in the exact definition of
 *                            time zero.
 *         implementation: 
 *                        A pseudo body is defined as either a single star or
 *                        the center of mass of a binary star (members of which
 *                        may or may not be composite themselves). First an
 *                        array of pseudo bodies is constructed, starting with
 *                        the root (the unbound "binary" consisting of target
 *                        and  projectile), the target and the projectile, in 
 *                        that order. The members of each binary present are
 *                        then added subsequently to the pseudo body array.
 *                        Finally the single stars are copied to the array
 *                        of real bodies.
 *             body table:
 *                        The body table has three entries for each body.
 *                        If a body has an array number  address , i.e. if it
 *                        is stored in  pseudo_bodies[address] , then:
 *                        pseudo_body_table[address][0] contains the character
 *                        of the body, either UNKNOWN, or SINGLE, or BINARY.
 *                        In the latter case, pseudo_body_table[address][1]
 *                        contains the address of the first member of the
 *                        binary, while pseudo_body_table[address][2] contains
 *                        the address of the second member of the binary.
 *                 prompt:
 *                        When the flag "prompt" is true, the procedure will
 *                        put the requests for the successive input values on
 *                        the screen. 
 *                        When the flag is false, nothing will be written; this
 *                        may be more convenient for input from a file (by
 *                        redirecting the standard input).
 *-----------------------------------------------------------------------------
 */

local systptr  setup_scene()
    {
    int  n_pseudobody;
    int  nbody;
    int  pseudo_body_table[NPBODY_MAX][3];
    real  time_of_projected_encounter;     /* between projectile and target  */
    real  time_before_projected_encounter; /* when scattering expt. starts   */
    body  pseudo_bodies[NPBODY_MAX];       /* pointer to pseudobody array    */
    systptr  new_system;                   /* to be initialized              */
    systptr  mk_system();
    bool  prompt;                          /* flag for screen output control */

    if (NDIM != 3)
	error("setup_scene: NDIM = %d != 3\n", NDIM);

    prompt = getbparam("prompt");

    nbody = n_pseudobody = 0;                     /* start with empty space */

    setup_empty_root(&n_pseudobody, pseudo_body_table);

    setup_target_external_orbit(pseudo_bodies, &n_pseudobody, &nbody,
				pseudo_body_table);

    setup_projectile_external_orbit(pseudo_bodies, &n_pseudobody, &nbody,
                                    &time_of_projected_encounter,
                                    &time_before_projected_encounter, 
                                    pseudo_body_table, prompt);

    setup_root_orbit(pseudo_bodies, n_pseudobody);

    transform_to_com_coordinates(pseudo_bodies, n_pseudobody);

    setup_target_internal_structure(pseudo_bodies, &n_pseudobody, &nbody,
                                    time_of_projected_encounter,
                                    time_before_projected_encounter, 
                                    pseudo_body_table, prompt);

    setup_projectile_internal_structure(pseudo_bodies, &n_pseudobody, &nbody,
                                        time_before_projected_encounter, 
                                        pseudo_body_table, prompt);
/*
 * transform from {z,x,y} coordinates, as provided by transkepler(),
 * to {x,y,z} coordinates:
 */
    kep2norm(pseudo_bodies, n_pseudobody);
/*
 * finally, package the results:
 */
    new_system = mk_system(nbody);             /* in file  systemalgebra.c  */
    Tnow(new_system) = - time_before_projected_encounter;
    Nbody(new_system) = nbody;
    make_body_array(Bodies(new_system), pseudo_bodies, nbody, n_pseudobody,
		    pseudo_body_table);
    return(new_system);
    }

/*-----------------------------------------------------------------------------
 *  setup_empty_root  --  since the root is the first pseudobody, the slot
 *                         pseudo_bodies[0]  has to be reserved before the
 *                        target and projectiles are installed in external
 *                        orbit form.  The relevant root information is added
 *                        later by the procedure  insert_root_orbit() .
 *-----------------------------------------------------------------------------
 */
local void  setup_empty_root(n_pseudobody, pseudo_body_table)
int  *n_pseudobody;
int  pseudo_body_table[NPBODY_MAX][3];
    {
    if (*n_pseudobody != 0)
	error("setup_empty_root: *n_pseudobody = %d != 0\n",
	      *n_pseudobody);

    pseudo_body_table[ROOT_ADDRESS][0] = BINARY;
    pseudo_body_table[ROOT_ADDRESS][1] = TARGET_ADDRESS;
    pseudo_body_table[ROOT_ADDRESS][2] = PROJECTILE_ADDRESS;
    
    *n_pseudobody +=1;

    if (debug)
	internal_structure_message(ROOT_ADDRESS, pseudo_body_table);
    }

/*-----------------------------------------------------------------------------
 *  setup_target_external_orbit  --  initializes the center of mass system
 *                                   of the target object, which by convention
 *                                   is centered on the origin of the
 *                                   coordinate system.
 *-----------------------------------------------------------------------------
 */
local void  setup_target_external_orbit(pseudo_bodies, n_pseudobody, nbody,
					pseudo_body_table)
bodyptr  pseudo_bodies;
int  *n_pseudobody;
int  *nbody;
int  pseudo_body_table[NPBODY_MAX][3];
    {
    real  orbital_elements[2*NDIM];

    if (*n_pseudobody != 1)
	error("setup_target_external_orbit: *n_pseudobody = %d != 1\n",
	      *n_pseudobody);	
    if (*nbody != 0)
	error("setup_target_external_orbit: *nbody = %d != 0\n", *nbody);

    orbital_elements[0] = TARGET_COM_RX;
    orbital_elements[1] = TARGET_COM_RY;
    orbital_elements[2] = TARGET_COM_RZ;
    orbital_elements[3] = TARGET_COM_VX;
    orbital_elements[4] = TARGET_COM_VY;
    orbital_elements[5] = TARGET_COM_VZ;

    setup_body_orbit(pseudo_bodies + *n_pseudobody,
		     TARGET_MASS, orbital_elements);

    *n_pseudobody +=1;
    *nbody +=1;

    pseudo_body_table[TARGET_ADDRESS][0] = UNKNOWN;
    }

/*-----------------------------------------------------------------------------
 *  setup_projectile_external_orbit  --  initializes the mass and orbital
 *                                       parameters for the motion of the
 *                                       center of mass of a the projectile.
 *                      note: When the flag  prompt  is true, the procedure 
 *                            will put the requests for the successive input
 *                            values on the screen. 
 *                            When the flag is false, nothing will be written;
 *                            this may be more convenient for input from a file
 *                            (by redirecting the standard input).
 *-----------------------------------------------------------------------------
 */
local void  setup_projectile_external_orbit(pseudo_bodies, n_pseudobody, nbody,
                                            time_of_projected_encounter,
                                            time_before_projected_encounter, 
                    	                    pseudo_body_table, prompt)
bodyptr  pseudo_bodies;
int  *n_pseudobody;
int  *nbody;
realptr  time_of_projected_encounter;
realptr  time_before_projected_encounter;
int  pseudo_body_table[NPBODY_MAX][3];
bool  prompt;
    {
    real  projectile_mass, total_mass;
    real  ratio;
    real  orbital_elements[2*NDIM];
    real  initial_separation;

    if (*n_pseudobody != 2)
	error("setup_projectile_external_orbit: *n_pseudobody = %d != 2\n",
	      *n_pseudobody);	
    if (*nbody != 1)
	error("setup_target_external_orbit: *nbody = %d != 1\n", *nbody);

    if (prompt)
	printf("==> Projectile: External Orbit Specification:\n");

    if (prompt)
	printf("  [ projectile mass ] / [ target mass ] = ");
    scanf("%lf", &ratio);
    if (ratio < 0.0)
	error("setup_projectile_external_orbit: negative mass\n");	
    projectile_mass = ratio * TARGET_MASS;

    if (prompt)
	printf("  time of projected encounter: T_pericenter = ");
    scanf("%lf", time_of_projected_encounter);

    if (prompt)
	printf("  impact parameter:  rho = ");
    scanf("%lf", orbital_elements);
    if (orbital_elements[0] < 0.0)
	error("setup_projectile_external_orbit: negative impact parameter\n");

    if (prompt)
	printf("  impact angle:  psi = ");
    scanf("%lf", orbital_elements + 1);

    if (prompt)
	printf("  initial separation between projectile and target: r_in = ");
    scanf("%lf", orbital_elements + 2);
    if (orbital_elements[2] <= 0.0)
	error("setup_projectile_external_orbit: initial separation <= 0.0\n");

    if (prompt)
	printf("  asymptotic incoming velocity:  v_inf = ");
    scanf("%lf", orbital_elements + 3);
    if (orbital_elements[3] < 0.0)
	error("setup_projectile_external_orbit: velocity at infinity < 0.0\n");

    if (prompt)
	printf("  incoming direction spherical polar angle:  theta = ");
    scanf("%lf", orbital_elements + 4);

    if (prompt)
	printf("  incoming direction spherical azimuth angle:  phi = ");
    scanf("%lf", orbital_elements + 5);

    total_mass = projectile_mass + TARGET_MASS;

    initial_separation = orbital_elements[2];
    transkepler(total_mass, orbital_elements, SCATTERING, TSCATTERING);
    *time_before_projected_encounter = orbital_elements[2];
    orbital_elements[2] = initial_separation;

    transkepler(total_mass, orbital_elements, SCATTERING, CARTESIAN);

    if (prompt)
	print_cartesian(orbital_elements);

    setup_body_orbit(pseudo_bodies + *n_pseudobody,
                     projectile_mass, orbital_elements);

    *n_pseudobody +=1;
    *nbody +=1;

    pseudo_body_table[PROJECTILE_ADDRESS][0] = UNKNOWN;
    }

/*-----------------------------------------------------------------------------
 *  setup_target_internal_structure  --  reads in optional binary structure
 *                                       of the target object
 *-----------------------------------------------------------------------------
 */
local void  setup_target_internal_structure(pseudo_bodies, n_pseudobody, nbody,
                                    time_of_projected_encounter,
                                    time_before_projected_encounter, 
				    pseudo_body_table, prompt)
bodyptr  pseudo_bodies;
int  *n_pseudobody;
int  *nbody;
real  time_of_projected_encounter;
real  time_before_projected_encounter;
int  pseudo_body_table[NPBODY_MAX][3];
bool  prompt;
    {
    char  c;
    int  primary_address;
    int  secondary_address;
    real  secondary_mass_fraction;
    real  primary_mass;
    real  secondary_mass;
    real  orbital_elements[2*NDIM];
    bodyptr  target;
    body  primary;
    body  secondary;

    if (*n_pseudobody != 3)
	error("setup_target_internal_structure: *n_pseudobody = %d != 3\n",
	      *n_pseudobody);
    if (*nbody != 2)
	error("setup_target_internal_structure: *nbody = %d != 2\n", *nbody);

    if (pseudo_body_table[TARGET_ADDRESS][0] != UNKNOWN)
	error("setup_target_internal_structure: TARGET not of type UNKNOWN\n");

    if (prompt)
	printf("is the target a binary? (y/n)  ");
    while ((c = getchar()) != 'y' && c != 'n') 
	    ;
    if (c == 'n')           /* single star? */
	{
        pseudo_body_table[TARGET_ADDRESS][0] = SINGLE;
	return;             /* then return, leaving the body count unchanged */
                            /* otherwise, initialize binary components: */
	}
    if (prompt)
	printf("==> Target: Internal Structure Specification:\n");

    target = pseudo_bodies + TARGET_ADDRESS;

    if (prompt)
	printf("  mass fraction in secondary:  m[2]/(m[1]+m[2]) = ");
    scanf("%lf", &secondary_mass_fraction);
    if (secondary_mass_fraction < 0.0)
	error("setup_target_internal_structure: secondary_mass_fraction <0\n");
    if (secondary_mass_fraction > 1.0)
	error("setup_target_internal_structure: secondary_mass_fraction >1\n");
    primary_mass = (1.0 - secondary_mass_fraction) * TARGET_MASS;
    secondary_mass = secondary_mass_fraction * TARGET_MASS;

    orbital_elements[0] = TARGET_SEMIMAJOR_AXIS;

    if (prompt)
	printf("  eccentricity:  e = ");
    scanf("%lf", orbital_elements + 1);
    if (orbital_elements[1] < 0.0)
	error("setup_target_internal_structure: eccentricity < 0.0\n");
    if (orbital_elements[1] > 1.0)
	error("setup_target_internal_structure: eccentricity > 1.0\n");

    orbital_elements[2] = TARGET_INCLINATION;

    orbital_elements[3] = TARGET_LONGITUDE_OF_ASCENDING_NODE;

    orbital_elements[4] = TARGET_ARGUMENT_OF_PERICENTER;

    orbital_elements[5] = TARGET_TIME_OF_PERICENTER_PASSAGE;
                                                   /* in original time frame */
    orbital_elements[5] -= time_of_projected_encounter;
            /* in time frame in which the projected encounter is at time 0.0 */
    orbital_elements[5] += time_before_projected_encounter;
                /* in time frame in which start of simulation is at time 0.0 */

    transkepler(TARGET_MASS, orbital_elements, PERIPASSAGE, CARTESIAN);

    if (prompt)
	print_cartesian(orbital_elements);

    setup_offset(&primary, &secondary, orbital_elements,
		 primary_mass, secondary_mass);

    add_to_com(&primary, target);
    add_to_com(&secondary, target);

    insert_body_orbit(pseudo_bodies + *n_pseudobody, &primary);
    *n_pseudobody +=1;

    insert_body_orbit(pseudo_bodies + *n_pseudobody, &secondary);
    *n_pseudobody +=1;

    *nbody +=1;

    primary_address = *n_pseudobody - 2;
    secondary_address = *n_pseudobody - 1;
    pseudo_body_table[TARGET_ADDRESS][0] = BINARY;
    pseudo_body_table[TARGET_ADDRESS][1] = primary_address;
    pseudo_body_table[TARGET_ADDRESS][2] = secondary_address;
    pseudo_body_table[primary_address][0] = UNKNOWN;
    pseudo_body_table[secondary_address][0] = UNKNOWN;

    if (debug)
	internal_structure_message(TARGET_ADDRESS, pseudo_body_table);

    read_in_internal_structure(pseudo_bodies, primary_address, n_pseudobody,
			       nbody, time_before_projected_encounter, 
			       pseudo_body_table, prompt);
    read_in_internal_structure(pseudo_bodies, secondary_address, n_pseudobody,
			       nbody, time_before_projected_encounter, 
			       pseudo_body_table, prompt);
    }

/*-----------------------------------------------------------------------------
 *  setup_projectile_internal_structure  --  reads in optional binary structure
 *                                           of the projectile
 *-----------------------------------------------------------------------------
 */
local void  setup_projectile_internal_structure(pseudo_bodies,
						n_pseudobody, nbody,
                                               time_before_projected_encounter,
				               pseudo_body_table, prompt)
bodyptr  pseudo_bodies;
int  *n_pseudobody;
int  *nbody;
real  time_before_projected_encounter;
int  pseudo_body_table[NPBODY_MAX][3];
bool  prompt;
    {
    read_in_internal_structure(pseudo_bodies, PROJECTILE_ADDRESS, n_pseudobody,
			       nbody, time_before_projected_encounter, 
			       pseudo_body_table, prompt);
    }

/*-----------------------------------------------------------------------------
 *  read_in_internal_structure  --  
 *-----------------------------------------------------------------------------
 */
local void  read_in_internal_structure(pseudo_bodies, object_address,
				       n_pseudobody, nbody,
                                       time_before_projected_encounter,
				       pseudo_body_table, prompt)
bodyptr  pseudo_bodies;
int  object_address;
int  *n_pseudobody;
int  *nbody;
real  time_before_projected_encounter;
int  pseudo_body_table[NPBODY_MAX][3];
bool  prompt;
    {
    char  c;
    int  primary_address;
    int  secondary_address;
    real  secondary_mass_fraction;
    real  primary_mass;
    real  secondary_mass;
    real  orbital_elements[2*NDIM];
    bodyptr  object;
    body  primary;
    body  secondary;

    if (*n_pseudobody > NPBODY_MAX - 2)
	error("read_in_internal_structure: *n_pseudobody = %d >NPBODY_MAX-2\n",
	      *n_pseudobody);
    if (*nbody >= *n_pseudobody)
	error("setup_target_external_orbit: *nbody= %d >= *n_pseudobody= %d\n",
	      *nbody, *n_pseudobody);

    if (pseudo_body_table[object_address][0] != UNKNOWN)
	error("read_in_internal_structure: object not of type UNKNOWN\n");

    if (prompt)
	printf("is pseudo body #%d a binary? (y/n)  ", object_address);
    while ((c = getchar()) != 'y' && c != 'n') 
	    ;
    if (c == 'n')           /* single star? */
	{
        pseudo_body_table[object_address][0] = SINGLE;
	return;             /* then return, leaving the body count unchanged */
                            /* otherwise, initialize binary components: */
	}

    if (prompt)
	printf("==> Pseudo body #%d: Internal Structure Specification:\n",
	       object_address);

    object = pseudo_bodies + object_address;

    if (prompt)
	printf("  mass fraction in secondary:  m[2]/(m[1]+m[2]) = ");
    scanf("%lf", &secondary_mass_fraction);
    if (secondary_mass_fraction < 0.0)
	error("read_in_internal_structure: secondary_mass_fraction < 0.0\n");
    if (secondary_mass_fraction > 1.0)
	error("read_in_internal_structure: secondary_mass_fraction < 1.0\n");
    primary_mass = (1.0 - secondary_mass_fraction) * Mass(object);
    secondary_mass = secondary_mass_fraction * Mass(object);

    if (prompt)
	printf("  semimajor axis:  a = ");
    scanf("%lf", orbital_elements);
    if (orbital_elements[0] < 0.0)
	error("read_in_internal_structure: semimajor axis a < 0.0\n");

    if (prompt)
	printf("  eccentricity:  e = ");
    scanf("%lf", orbital_elements + 1);
    if (orbital_elements[1] < 0.0)
	error("read_in_internal_structure: eccentricity < 0.0\n");
    if (orbital_elements[1] > 1.0)
	error("read_in_internal_structure: eccentricity > 1.0\n");

    if (prompt)
	printf("  inclination:  i = ");
    scanf("%lf", orbital_elements + 2);

    if (prompt)
	printf("  longitude of the ascending node:  Omega = ");
    scanf("%lf", orbital_elements + 3);

    if (prompt)
	printf("  argument of pericenter:  omega = ");
    scanf("%lf", orbital_elements + 4);

    if (prompt)
	printf("  mean anomaly at projected time of encounter:  m = ");
    scanf("%lf", orbital_elements + 5);

    transkepler(Mass(object), orbital_elements, MEANANOMALY, PERIPASSAGE);

    orbital_elements[5] += time_before_projected_encounter;
                /* in time frame in which start of simulation is at time 0.0 */

    transkepler(Mass(object), orbital_elements, PERIPASSAGE, CARTESIAN);

    if (prompt)
	print_cartesian(orbital_elements);

    setup_offset(&primary, &secondary, orbital_elements,
		 primary_mass, secondary_mass);

    add_to_com(&primary, object);
    add_to_com(&secondary, object);

    insert_body_orbit(pseudo_bodies + *n_pseudobody, &primary);
    *n_pseudobody +=1;

    insert_body_orbit(pseudo_bodies + *n_pseudobody, &secondary);
    *n_pseudobody +=1;

    *nbody +=1;

    primary_address = *n_pseudobody - 2;
    secondary_address = *n_pseudobody - 1;
    pseudo_body_table[object_address][0] = BINARY;
    pseudo_body_table[object_address][1] = primary_address;
    pseudo_body_table[object_address][2] = secondary_address;
    pseudo_body_table[primary_address][0] = UNKNOWN;
    pseudo_body_table[secondary_address][0] = UNKNOWN;

    if (debug)
	internal_structure_message(object_address, pseudo_body_table);

    read_in_internal_structure(pseudo_bodies, primary_address, n_pseudobody,
			       nbody, time_before_projected_encounter, 
			       pseudo_body_table, prompt);
    read_in_internal_structure(pseudo_bodies, secondary_address, n_pseudobody,
			       nbody, time_before_projected_encounter, 
			       pseudo_body_table, prompt);
    }

/*-----------------------------------------------------------------------------
 *  make_body_array  --  from among the array of pseudo bodies, copy the real
 *                       bodies to the body array. The numbering is in the
 *                       order of the body array, except that first all the
 *                       target bodies are listed, and then all the projectile
 *                       bodies.
 *-----------------------------------------------------------------------------
 */
local void  make_body_array(some_bodies, pseudo_bodies, nbody, n_pseudobody,
			    pseudo_body_table)
bodyptr  some_bodies;
bodyptr  pseudo_bodies;
int  nbody;
int  n_pseudobody;
int  pseudo_body_table[NPBODY_MAX][3];
    {
    int  i;                                           /* body address        */
    int  j;                                           /* pseudo body address */

    if (n_pseudobody < 3)
	error("make_body_array: n_pseudobody = %d < 3\n", n_pseudobody);
    if (nbody >= n_pseudobody)
	error("make_body_array: nbody= %d >= n_pseudobody= %d\n",
	      nbody, n_pseudobody);

    i = j = 0;
    for (j = 0; j < n_pseudobody; j++)
	{
	if (pseudo_body_table[j][0] == UNKNOWN)
	    error("make_body_array: pseudo_body_table[%d][0] == UNKNOWN\n", j);
	if (pseudo_body_table[j][0] == SINGLE)
	    if (j != PROJECTILE_ADDRESS)
		{
		insert_body_orbit(some_bodies + i, pseudo_bodies + j);
		i++;
		}
	}

    if (pseudo_body_table[PROJECTILE_ADDRESS][0] == SINGLE)
	{
	insert_body_orbit(some_bodies + i, pseudo_bodies + PROJECTILE_ADDRESS);
	i++;
	}

    if (i != nbody)
	error("make_body_array: i = %d != nbody = %d\n", i, nbody);
    }

/*-----------------------------------------------------------------------------
 *  setup_root_orbit  --  initializes the mass and orbital parameters for the
 *                        root pseudo body, in a Cartesian coordinate system
 *-----------------------------------------------------------------------------
 */
local void  setup_root_orbit(some_body, n)
bodyptr  some_body;
int  n;
    {
    real  target_weighted_pos[NDIM];
    real  target_weighted_vel[NDIM];
    real  projectile_weighted_pos[NDIM];
    real  projectile_weighted_vel[NDIM];
    
    if (n < 3)
	error("setup_root_orbit: n = %d < 3\n", n);

    Mass(some_body) = Mass(some_body + TARGET_ADDRESS)
	              + Mass(some_body + PROJECTILE_ADDRESS);

    MULVS(target_weighted_pos,
	  Pos(some_body + TARGET_ADDRESS),
	  Mass(some_body + TARGET_ADDRESS));
    MULVS(target_weighted_vel,
	  Vel(some_body + TARGET_ADDRESS),
	  Mass(some_body + TARGET_ADDRESS));
    MULVS(projectile_weighted_pos,
	  Pos(some_body + PROJECTILE_ADDRESS),
	  Mass(some_body + PROJECTILE_ADDRESS));
    MULVS(projectile_weighted_vel,
	  Vel(some_body + PROJECTILE_ADDRESS),
	  Mass(some_body + PROJECTILE_ADDRESS));

    ADDV(Pos(some_body), target_weighted_pos, projectile_weighted_pos);
    ADDV(Vel(some_body), target_weighted_vel, projectile_weighted_vel);

    INCDIVVS(Pos(some_body), Mass(some_body));
    INCDIVVS(Vel(some_body), Mass(some_body));
/*
 * no pseudo_body_table update necessary here, since  setup_empty_root()  has
 * already taken care of that before.
 */
    }

/*-----------------------------------------------------------------------------
 *  setup_body_orbit  --  initializes the mass and orbital parameters for a
 *                        single body, in a Cartesian coordinate system
 *-----------------------------------------------------------------------------
 */
local void  setup_body_orbit(some_body, some_mass, pos_and_vel)
bodyptr  some_body;
real  some_mass;
realptr  pos_and_vel;
    {
    Mass(some_body) = some_mass;
    SETV(Pos(some_body), pos_and_vel);
    SETV(Vel(some_body), pos_and_vel + NDIM);
    }

/*-----------------------------------------------------------------------------
 *  insert_body_orbit  --  copies the mass and orbital parameters, for a
 *                         single body
 *-----------------------------------------------------------------------------
 */
local void  insert_body_orbit(to_body, from_body)
bodyptr  to_body, from_body;
    {
    Mass(to_body) = Mass(from_body);
    SETV(Pos(to_body), Pos(from_body));
    SETV(Vel(to_body), Vel(from_body));
    }

/*-----------------------------------------------------------------------------
 *  add_to_com  --  adds the center-of-mass orbit to a member of a binary
 *-----------------------------------------------------------------------------
 */
local void  add_to_com(member, com)
bodyptr  member;
bodyptr  com;
    {
    INCADDV(Pos(member), Pos(com));
    INCADDV(Vel(member), Vel(com));
    }

/*-----------------------------------------------------------------------------
 *  transform_to_com_coordinates  --  transforms to a coordinate system in
 *                                    which the center of mass remains at
 *                                    rest in the origin.
 *-----------------------------------------------------------------------------
 */
local void  transform_to_com_coordinates(pseudo_bodies, n_pseudobody)
bodyptr  pseudo_bodies;
int  n_pseudobody;
    {
    real  pos_com[NDIM];
    real  vel_com[NDIM];
    bodyptr  body_i;

    SETV(pos_com, Pos(pseudo_bodies));
    SETV(vel_com, Vel(pseudo_bodies));

    for (body_i = pseudo_bodies; body_i-pseudo_bodies < n_pseudobody; body_i++)
        {
	INCSUBV(Pos(body_i), pos_com);
	INCSUBV(Vel(body_i), vel_com);
        }
    }

/*-----------------------------------------------------------------------------
 *  setup_offset  --  initializes the mass, and orbital parameters with respect
 *                    to the center of mass, for each of two bodies, starting
 *                    with the relative orbital parameters for the separation
 *                    vector pointing from the primary to the secondary star
 *                    in a binary, in a Cartesian coordinate system.
 *-----------------------------------------------------------------------------
 */
local void  setup_offset(primary, secondary, pos_and_vel,
	         	 primary_mass, secondary_mass)
bodyptr  primary;
bodyptr  secondary;
realptr  pos_and_vel;
real  primary_mass;
real  secondary_mass;
    {
    real  scale_factor;
    real  total_mass;
    real  helper_pos_and_vel[2*NDIM];

    Mass(primary) = primary_mass;
    Mass(secondary) = secondary_mass;

    total_mass = primary_mass + secondary_mass;
    scale_factor = primary_mass / total_mass;

    SETV(helper_pos_and_vel, pos_and_vel);
    SETV(helper_pos_and_vel + NDIM, pos_and_vel + NDIM);
    INCMULVS(helper_pos_and_vel, scale_factor);
    INCMULVS(helper_pos_and_vel + NDIM, scale_factor);

    SETV(Pos(secondary), helper_pos_and_vel);
    SETV(Vel(secondary), helper_pos_and_vel + NDIM);

    scale_factor = scale_factor - 1.0;                        /*  < 0  */

    SETV(helper_pos_and_vel, pos_and_vel);
    SETV(helper_pos_and_vel + NDIM, pos_and_vel + NDIM);
    INCMULVS(helper_pos_and_vel, scale_factor);
    INCMULVS(helper_pos_and_vel + NDIM, scale_factor);

    SETV(Pos(primary), helper_pos_and_vel);
    SETV(Vel(primary), helper_pos_and_vel + NDIM);
    }

/*-----------------------------------------------------------------------------
 *  internal_structure_message  --  prints the composition of a pseudo particle
 *-----------------------------------------------------------------------------
 */
local void internal_structure_message(binary_address, pseudo_body_table)
int  binary_address;
int  pseudo_body_table[NPBODY_MAX][3];
    {
    int  first_member_address;
    int  second_member_address;

    if (pseudo_body_table[binary_address][0] != BINARY)
	error("internal_structure_message: pseudo body is not a binary\n");

    first_member_address = pseudo_body_table[binary_address][1];
    second_member_address = pseudo_body_table[binary_address][2];

    if (binary_address == ROOT_ADDRESS)
        {
	printf("  the ROOT (pseudo particle #%d)", ROOT_ADDRESS);
        printf(" has received two internal components:\n");
	}
    else if (binary_address == TARGET_ADDRESS)
        {
	printf("  the TARGET (pseudo particle #%d)", TARGET_ADDRESS);
        printf(" has received two internal components:\n");
	}
    else if (binary_address == PROJECTILE_ADDRESS)
        {
	printf("  the PROJECTILE (pseudo particle #%d)", PROJECTILE_ADDRESS);
        printf(" has received two internal components:\n");
	}
    else
	printf("  pseudo particle #%d has received two internal components:\n",
	       binary_address);

    if (first_member_address == ROOT_ADDRESS)
	error("internal_structure_message: ROOT cannot be a binary member\n");
    else if (first_member_address == TARGET_ADDRESS)
	printf("    the TARGET (pseudo particle #%d)  &  ", TARGET_ADDRESS);
    else if (first_member_address == PROJECTILE_ADDRESS)
	printf("    the PROJECTILE (pseudo particle #%d)  &  ",
	       PROJECTILE_ADDRESS);
    else
	printf("    pseudo particle #%d  &  ", first_member_address);

    if (second_member_address == ROOT_ADDRESS)
	error("internal_structure_message: ROOT cannot be a binary member\n");
    else if (second_member_address == TARGET_ADDRESS)
	printf("the TARGET (pseudo particle #%d)\n", TARGET_ADDRESS);
    else if (second_member_address == PROJECTILE_ADDRESS)
	printf("the PROJECTILE (pseudo particle #%d)\n",
	       PROJECTILE_ADDRESS);
    else
	printf("pseudo particle #%d\n", second_member_address);
    }

/*-----------------------------------------------------------------------------
 *  kep2norm  --  transform from {z,x,y} coordinates, as provided by
  *               transkepler(), to {x,y,z} coordinates.
 *-----------------------------------------------------------------------------
 */
local void  kep2norm(some_bodies, n)
bodyptr  some_bodies;
int  n;
    {
    bodyptr  body_i;
    real  z_component;

    if (n <0)
	error("kep2norm: particle number is negative\n");

    for (body_i = some_bodies; body_i - some_bodies < n; body_i++)
	{
	z_component = Pos(some_bodies)[0];
	Pos(some_bodies)[0] = Pos(some_bodies)[1];
	Pos(some_bodies)[1] = Pos(some_bodies)[2];
	Pos(some_bodies)[2] = z_component;

	z_component = Vel(some_bodies)[0];
	Vel(some_bodies)[0] = Vel(some_bodies)[1];
	Vel(some_bodies)[1] = Vel(some_bodies)[2];
	Vel(some_bodies)[2] = z_component;
	}
    }

/*-----------------------------------------------------------------------------
 *  print_cartesian_elements  --  print a kepler orbit in Cartesian coordinates
 *-----------------------------------------------------------------------------
 */
local void  print_cartesian(elementptr)
realptr  elementptr;
    {
    printf("\n");
    printf("    x = %f    y = %f    z = %f\n",
           elementptr[1], elementptr[2], elementptr[0]);
    printf("   vx = %f   vy = %f   vz = %f\n",
           elementptr[4], elementptr[5], elementptr[3]);
    printf("\n");
    }

/* endof: scatter.c */
