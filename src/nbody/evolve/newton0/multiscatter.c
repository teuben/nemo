/* multiscatter.c - */

/*
 *  multiscatter.c:  orchestrates a simple series of scattering events
 *
 *     July 1988  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
#include  "newton0.h"                /* which includes  state.h  */
#include  "multiscat.h"              /* which includes  scat.h   */
#include  <getparam.h>

#define  NBODY_MAX    10                     /* maximum number of particles  */

bool  debug;              /* flag for output control, for debugging purposes */
bool  turbo;              /* flag for integration spee-up, through shortcuts */

/*-----------------------------------------------------------------------------
 *  multiscatter: 
 *  TO BE DOCUMENTED   TO BE DOCUMENTED    TO BE DOCUMENTED    TO BE DOCUMENTED
 * 
 *  simple version which performes a series of scattering experiments with
 *  parameters chosen randomly within given boundaries, currently haphazardly
 *  implemented just for three- and four-body scattering.
 *
 *  TO BE DOCUMENTED   TO BE DOCUMENTED    TO BE DOCUMENTED    TO BE DOCUMENTED
 *           performs a scattering event with one target, located originally
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
    "de_rel_max=0.001",     /* maximum allowed rel drift in total energy     */
    "nstep_de=100",         /* number of steps between checking  de_rel_max  */
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
    "prompt=FALSE",         /* for interactive parameter input from terminal */
    "debug=TRUE",           /* output information for debugging purposes     */
    "turbo=TRUE",           /* integration speed-up through shortcuts        */
    "headline=",            /* verbiage for output                           */
    NULL,
    };

/*-----------------------------------------------------------------------------
 *  main  --  read the command line arguments, pass control to  multiscatter()
 *-----------------------------------------------------------------------------
 */
main(argc, argv)
int  argc;
string  argv[];
    {
    initparam(argv, defv);     		           /* setup parameter access */

    debug = getbparam("debug");
    turbo = getbparam("turbo");

    multiscatter();
    }

/*-----------------------------------------------------------------------------
 *  multiscatter  --  orchestrates a simple series of scattering experiments
 *                    preliminary version: just to make it work, not yet
 *                    any attempt to make it modular, etc.
 *-----------------------------------------------------------------------------
 */
local void  multiscatter()
    {
    int  i;
    int  n_experiments;
    int  seed;
    int  nbody;
    bool  answers[2*NBODY_MAX-2];    /* information about binary hierarchies */
    real  mproj;                     /* projectile mass                      */
    real  v_inf;                 /* target-projectile rel. vel. at infinity  */
    real  mfrac[NBODY_MAX-2];    /* mass fraction of secondaries in binaries */
    real  a[NBODY_MAX-2];        /* semimajor axes of binaries               */
    real  e[NBODY_MAX-2];        /* eccentricities of binaries               */
    real  rho_min;
    real  rho_max;
    real  r_in;
    multiclass  report;

    read_physical_parameters(&nbody, answers, &mproj, &v_inf, mfrac, a, e);

    read_run_control_parameters(&n_experiments, &rho_min, &rho_max, &r_in,
				&seed);

    setup_a_series_of_experiments(seed, &report);

    run_a_series_of_experiments(nbody, n_experiments, answers, mproj, v_inf,
				mfrac, a, e, rho_min, rho_max, r_in, &report);

    report_a_series_of_experiments(&report);
    }

/*-----------------------------------------------------------------------------
 *  read_physical_parameters  --  
 *                           note: if negative eccentricity values are given,
 *                                 that will be interpreted as requesting 
 *                                 random eccentricity values, drawn from an
 *                                 isothermal  f(e) = 2e  distribution.
 *-----------------------------------------------------------------------------
 */
local void  read_physical_parameters(nbody, answers, mproj, v_inf, mfrac, a, e)
int *nbody;                      /* total number of particles                */
bool  answers[2*NBODY_MAX-2];    /* information about binary hierarchies     */
realptr  mproj;                  /* projectile mass                          */
realptr  v_inf;                  /* target-projectile rel. vel. at infinity  */
real  mfrac[NBODY_MAX-2];        /* mass fraction of secondaries in binaries */
real  a[NBODY_MAX-2];            /* semimajor axes of binaries               */
real  e[NBODY_MAX-2];            /* eccentricities of binaries               */
    {
    char  c;
    bool  prompt;                          /* flag for screen output control */
    int  answer_counter;
    int  mfrac_counter;
    int  a_counter;
    int  e_counter;
    int  how_many_questions_left;

    *nbody = 0;
    how_many_questions_left = 2;
    answer_counter = 0;
    mfrac_counter = a_counter = e_counter = 0;

    if (NDIM != 3)
	error("read_physical_parameters: NDIM = %d != 3\n", NDIM);

    prompt = getbparam("prompt");

    if (prompt)
	printf("==> Projectile: External Orbit Specification:\n");

    if (prompt)
	printf("  [ projectile mass ] / [ target mass ] = ");
    scanf("%lf", mproj);
    if (*mproj < 0.0)
	error("read_physical_parameters: negative projectile mass\n");	

    if (prompt)
	printf("  asymptotic incoming velocity:  v_inf = ");
    scanf("%lf", v_inf);
    if (*v_inf < 0.0)
	error("read_physical_parameters: velocity at infinity < 0.0\n");

    *nbody += 2;

    if (prompt)
	printf("is the target a binary? (y/n)  ");
    while ((c = getchar()) != 'y' && c != 'n') 
	;
    if (c == 'n')           /* single star? */
	{
	answers[answer_counter++] = FALSE;
	how_many_questions_left--;
	}
    else
	{
	answers[answer_counter++] = TRUE;
	how_many_questions_left++;
	(*nbody)++;

	if (prompt)
	    printf("==> Target: Internal Structure Specification:\n");

	if (prompt)
	    printf("  mass fraction in secondary:  m[2]/(m[1]+m[2]) = ");
	scanf("%lf", mfrac + mfrac_counter++);
	if (*(mfrac + mfrac_counter - 1) < 0.0)
	    error("read_physical_parameters: target mfrac < 0\n");
	else if (*(mfrac + mfrac_counter - 1) > 1.0)
	    error("read_physical_parameters: target mfrac > 1\n");

	if (prompt)
	    printf("  eccentricity:  e = ");
	scanf("%lf", e + e_counter++);
	if (*(e + e_counter - 1) > 1.0)       /* e < 0 for isothermal values */
	    error("read_physical_parameters: target eccentricity > 1.0\n");
	}

    while (how_many_questions_left > 0)
	{
	if (prompt)
	    printf("is the next node a binary? (y/n)  ");
	while ((c = getchar()) != 'y' && c != 'n') 
	    ;
	if (c == 'n')           /* single star? */
	    {
	    answers[answer_counter++] = FALSE;
	    how_many_questions_left--;
	    }
	else
	    {
	    answers[answer_counter++] = TRUE;
	    how_many_questions_left++;
	    if (++*nbody > NBODY_MAX)
		error("read_physical_parameters: nbody > %d\n", NBODY_MAX);

	    if (prompt)
	        printf("==> Next Node: Internal Structure Specification:\n");

	    if (prompt)
	        printf("  mass fraction in secondary: m_sec/(m_sec+m_prim)= ");
	    scanf("%lf", mfrac + mfrac_counter++);
	    if (*(mfrac + mfrac_counter - 1) < 0.0)
	        error("read_physical_parameters: mfrac < 0\n");
	    else if (*(mfrac + mfrac_counter - 1) > 1.0)
	        error("read_physical_parameters: mfrac > 1\n");

	    if (prompt)
	        printf("  semimajor axis:  a = ");
	    scanf("%lf", a + a_counter++);
	    if (*(a + a_counter - 1) < 0.0)
	        error("read_physical_parameters: semimajor axis < 0.0\n");

	    if (prompt)
	        printf("  eccentricity:  e = ");
	    scanf("%lf", e + e_counter++);
	    if (*(e + e_counter - 1) > 1.0)   /* e < 0 for isothermal values */
	        error("read_physical_parameters: eccentricity > 1.0\n");
	    }
	}

    if (answer_counter != 2*(*nbody)-2)
	error("read_physical_parameters: answer_counter has illegal value\n");
    }

/*-----------------------------------------------------------------------------
 *  read_run_control_parameters  --  and also initializes the random number
 *                                   generator (with a random seed if *seed==0)
 *-----------------------------------------------------------------------------
 */
local void  read_run_control_parameters(n_experiments, rho_min, rho_max, r_in,
					seed)
int *n_experiments;
realptr  rho_min;
realptr  rho_max;
realptr  r_in;
int  *seed;
    {
    bool  prompt;                          /* flag for screen output control */

    if (NDIM != 3)
	error("read_run_control_parameters: NDIM = %d != 3\n", NDIM);

    prompt = getbparam("prompt");

    if (prompt)
	printf("==> Run Control Parameters:\n");

    if (prompt)
	printf("Number of experiments N = ");
    scanf("%d", n_experiments);
    if (*n_experiments < 0)
	error("read_run_control_parameters: N < 0\n");

    if (prompt)
	printf("  minimum impact parameter:  rho_min = ");
    scanf("%lf", rho_min);
    if (*rho_min < 0.0)
	error("read_run_control_parameters: minimum impact parameter < 0.0\n");

    if (prompt)
	printf("  maximum impact parameter:  rho_max = ");
    scanf("%lf", rho_max);
    if (*rho_max < 0.0)
	error("read_run_control_parameters: maximum impact parameter < 0.0\n");

    if (prompt)
	printf("  initial separation between projectile and target: r_in = ");
    scanf("%lf", r_in);
    if (*r_in <= 0.0)
	error("read_run_control_parameters: initial separation <= 0.0\n");

    if (prompt)
	printf("  for the random number generator: seed = ");
    scanf("%d", seed);
    if (*seed < 0)
	error("read_run_control_parameters: seed < 0\n");
    if (*seed == 0)         /* no particular positive seed provided?      */
	*seed = time(0);    /* give a random value, different each second */

    srandom(*seed);         /* initialize the random number generator     */
    }

/*-----------------------------------------------------------------------------
 *  setup_a_series_of_experiments  --  
 *-----------------------------------------------------------------------------
 */
local void  setup_a_series_of_experiments(seed, report)
int  seed;
multiclassptr  report;
    {
    int  i;
/*
 * for now, I will violate the rule that all ascii output occurs only in
 * the file  out.c :
 */
    printf("\n\tStarting a series of scattering experiments, from seed = %d\n",
	   seed);
/*
 * initialize the array of qualitative outcome strings and their frequencies:
 */
    for (i = 0; i < MAX_DIVERSITY; i++)
	{
	Scattype(report)[i][0] = NULL;
	Scatsubtotal(report)[i] = 0;
	}
    Nscattype(report) = 0;
    }

/*-----------------------------------------------------------------------------
 *  run_a_series_of_experiments  --  
 *-----------------------------------------------------------------------------
 */
local void  run_a_series_of_experiments(nbody, n_experiments, answers, mproj,
					v_inf, mfrac, a, e, rho_min, rho_max,
					r_in, report)
int  nbody;
int  n_experiments;
bool  answers[2*NBODY_MAX-2];        /* information about binary hierarchies */
real  mproj;                         /* projectile mass                      */
real  v_inf;                     /* target-projectile rel. vel. at infinity  */
real  mfrac[NBODY_MAX-2];        /* mass fraction of secondaries in binaries */
real  a[NBODY_MAX-2];            /* semimajor axes of binaries               */
real  e[NBODY_MAX-2];            /* eccentricities of binaries               */
real  rho_min;
real  rho_max;
real  r_in;
multiclassptr  report;
    {
    int  i;
/*
 * check whether all arguments have appropriate values:
 */
    if (NDIM != 3)
	error("run_a_series_of_experiments: NDIM = %d != 3\n", NDIM);

    if (nbody < 2)
	error("run_a_series_of_experiments: %d-body scattering impossible\n",
	      nbody);
    else if (nbody > NBODY_MAX)
	error("run_a_series_of_experiments: %d-body scattering not provided\n",
	      nbody);

    if (n_experiments < 1)
	error("run_a_series_of_experiments: wanna run %d experiments -- uh?\n",
	      n_experiments);

    if (mproj < 0.0)
	error("run_a_series_of_experiments: negative projectile mass\n");

    if (v_inf < 0.0)
	error("run_a_series_of_experiments: relative velocity at inf. < 0\n");

    if (mfrac == NULL && nbody > 2)
	error("run_a_series_of_experiments: mfrac = NULL (null pointer)\n");
    else if (a == NULL && (nbody > 3 || (nbody == 3 && answers[0] == FALSE)))
	error("run_a_series_of_experiments: a = NULL (null pointer)\n");
    else if (e == NULL && nbody > 2)
	error("run_a_series_of_experiments: e = NULL (null pointer)\n");

    if (rho_min < 0.0)
	error("run_a_series_of_experiments: rho_min < 0.0\n");
    else if (rho_max < rho_min)
	error("run_a_series_of_experiments: rho_max < rho_min\n");

    if (r_in < 0.0)
	error("run_a_series_of_experiments: r_in < 0.0\n");

    if (report == NULL)
	error("run_a_series_of_experiments: report = NULL (null pointer)\n");

/*
 * first run:
 */
    run_a_scattering_experiment(nbody, answers, mproj, v_inf, mfrac, a, e,
				rho_min, rho_max, r_in, TRUE, FALSE, report);
/*
 * second to one-to-last runs:
 */
    for (i = 1; i < n_experiments - 1; i++)
	run_a_scattering_experiment(nbody, answers, mproj, v_inf, mfrac, a, e,
				    rho_min, rho_max, r_in, FALSE, FALSE,
				    report);
/*
 * last run:
 */
    run_a_scattering_experiment(nbody, answers, mproj, v_inf, mfrac, a, e,
				rho_min, rho_max, r_in, FALSE, TRUE, report);
    }

/*-----------------------------------------------------------------------------
 *  report_a_series_of_experiments  --  
 *-----------------------------------------------------------------------------
 */
local void  report_a_series_of_experiments(report)
multiclassptr  report;
    {
    int  i;

    printf("\n  report of qualitative outcomes,");
    printf(" and their number of occurrences:\n\n");

    i = 0;
    while (i < MAX_DIVERSITY)
	{
	if (Scattype(report)[i][0] != NULL)
	    printf("     %s    %d\n", Scattype(report)[i],
		   Scatsubtotal(report)[i]);
	else
	    break;
	i++;
	}
    printf("\n");
    }

/*-----------------------------------------------------------------------------
 *  run_a_scattering_experiment  --
 *                             note: the case that the target is single and the
 *                                   projectile is composite is an unnatural
 *                                   choice, since the orbit parameters for the
 *                                   toplevel target binary are chosen so as to
 *                                   take out all ambiguity arising from
 *                                   scaling, homogeneity and isotropy; this
 *                                   advantage is lost with such an unnatural
 *                                   choice. However, it is a legal choice,
 *                                   which could be used for, e.g., debugging
 *                                   purposes, and it will work.
 *-----------------------------------------------------------------------------
 */
local void  run_a_scattering_experiment(nbody, answers, mproj, v_inf, mfrac, a,
					e, rho_min, rho_max, r_in, first_run,
					last_run, report)
int  nbody;
bool  answers[2*NBODY_MAX-2];        /* information about binary hierarchies */
real  mproj;                         /* projectile mass                      */
real  v_inf;                     /* target-projectile rel. vel. at infinity  */
real  mfrac[NBODY_MAX-2];        /* mass fraction of secondaries in binaries */
real  a[NBODY_MAX-2];            /* semimajor axes of binaries               */
real  e[NBODY_MAX-2];            /* eccentricities of binaries               */
real  rho_min;
real  rho_max;
real  r_in;
bool  first_run;
bool  last_run;
multiclassptr  report;
    {
    int  answer_counter;
    int  e_counter;
    int  internal_angles_counter;     /* = e_counter - 1,  unless target is  */
                                      /* single and projectile is composite, */
                                      /* in which (unnatural!) case          */
                                      /* internal_angles_counter = e_counter */
    int  how_many_questions_left;
    real  theta, phi, psi;
    real  t_enc;
    real  rho;
    real  xrandom();                  /* random number generator             */
    real  inclination[NBODY_MAX-3];
    real  longitude[NBODY_MAX-3];     /* of the ascending node               */
    real  argument[NBODY_MAX-3];      /* of pericenter                       */
    real  meananomaly[NBODY_MAX-3];

    answer_counter = 0;
    how_many_questions_left = 2;
    e_counter = internal_angles_counter = 0;
/*
 *  Projectile: remaining specification for the external orbit:
 */
    theta = acos(xrandom(-1.0, 1.0));
    phi = xrandom(0.0, TWO_PI);
    psi = xrandom(0.0, TWO_PI);

    t_enc = xrandom(0.0, TWO_PI);     /* period of target binary is T = 2 pi */
                                      /* and target daughter binaries are    */
                                      /* given a random phase below          */
    rho = sqrt(xrandom(rho_min*rho_min, rho_max*rho_max));

/*
 *  Target: remaining specification for the internal top-level binary orbit:
 */
    if (answers[answer_counter++] == FALSE)
	how_many_questions_left--;
    else
	{
	how_many_questions_left++;
	if (e[e_counter] < 0.0)
	    e[e_counter] = sqrt(xrandom(0.0, 1.0));
	e_counter++;
	}

/*
 *  Remaining specification for all other internal binary orbits:
 */
    while (how_many_questions_left > 0)
	{
	if (answer_counter >= 2*nbody-2)
	    error("run_a_scattering_experiment: answer_counter too large\n");

        if (answers[answer_counter++] == FALSE)
	    how_many_questions_left--;
        else
	    {
	    how_many_questions_left++;

	    if (e[e_counter] < 0.0)
	        e[e_counter] = sqrt(xrandom(0.0, 1.0));
	    e_counter++;

	    inclination[internal_angles_counter] = acos(xrandom(-1.0, 1.0));
	    longitude[internal_angles_counter] = xrandom(0.0, TWO_PI);
	    argument[internal_angles_counter] = xrandom(0.0, TWO_PI);
	    meananomaly[internal_angles_counter] = xrandom(0.0, TWO_PI);
	    internal_angles_counter++;
	    }
	}
	    
    one_scatter_experiment(nbody, answers, mproj, v_inf, theta, phi, psi,
			   t_enc, rho, r_in, mfrac, a, e, inclination,
			   longitude, argument, meananomaly, first_run,
			   last_run, report);
    }

/*-----------------------------------------------------------------------------
 *  one_scatter_experiment  --  orchestrates one gravitational scattering
 *                              experiment within a series.
 *-----------------------------------------------------------------------------
 */
local void  one_scatter_experiment(nbody, answers, mproj, v_inf, theta, phi,
				   psi, t_enc, rho, r_in, mfrac, a, e,
				   inclination, longitude, argument,
				   meananomaly, first_run, last_run, report)
int  nbody;
bool  answers[2*NBODY_MAX-2];        /* information about binary hierarchies */
real  mproj;                         /* projectile mass                      */
real  v_inf;                     /* target-projectile rel. vel. at infinity  */
real  theta;
real  phi;
real  psi;
real  t_enc;
real  rho;
real  r_in;
real  mfrac[NBODY_MAX-2];        /* mass fraction of secondaries in binaries */
real  a[NBODY_MAX-2];            /* semimajor axes of binaries               */
real  e[NBODY_MAX-2];            /* eccentricities of binaries               */
real  inclination[NBODY_MAX-3];
real  longitude[NBODY_MAX-3];    /* of the ascending node                    */
real  argument[NBODY_MAX-3];     /* of pericenter                            */
real  meananomaly[NBODY_MAX-3];
bool  first_run;
bool  last_run;
multiclassptr  report;
    {
    stateptr  new_state;
    stateptr  setup_one_experiment();

    new_state = setup_one_experiment(nbody, answers, mproj, v_inf, theta, phi,
				   psi, t_enc, rho, r_in, mfrac, a, e,
				   inclination, longitude, argument,
				   meananomaly);

    perform_one_experiment(new_state, first_run);

    report_one_experiment(new_state, first_run, last_run, report);

    clean_up_after_one_experiment(new_state);
    }

/*-----------------------------------------------------------------------------
 *  setup_one_experiment  --  set variables in diverse state substructures,
 *                            for one many-body scattering experiment, within
 *                            a series.
 *                      note: the fourth state substructure, containing the
 *                            diagnostics, cannot yet be set at this point,
 *                            since that would require energy diagnostics,
 *                            which requires revving the engines, which could
 *                            not have been done before the other three
 *                            substructures are set.
 *                            Therefore, the diagnostics are currently set in
 *                            the file  diagnose.c .
 *-----------------------------------------------------------------------------
 */
local stateptr  setup_one_experiment(nbody, answers, mproj, v_inf, theta, phi,
				   psi, t_enc, rho, r_in, mfrac, a, e,
				   inclination, longitude, argument,
				   meananomaly)
int  nbody;
bool  answers[2*NBODY_MAX-2];        /* information about binary hierarchies */
real  mproj;                         /* projectile mass                      */
real  v_inf;                     /* target-projectile rel. vel. at infinity  */
real  theta;
real  phi;
real  psi;
real  t_enc;
real  rho;
real  r_in;
real  mfrac[NBODY_MAX-2];        /* mass fraction of secondaries in binaries */
real  a[NBODY_MAX-2];            /* semimajor axes of binaries               */
real  e[NBODY_MAX-2];            /* eccentricities of binaries               */
real  inclination[NBODY_MAX-3];
real  longitude[NBODY_MAX-3];    /* of the ascending node                    */
real  argument[NBODY_MAX-3];     /* of pericenter                            */
real  meananomaly[NBODY_MAX-3];
    {
    stateptr  new_state;
    stateptr  mk_state();
    systptr  setup_one_scene();
    string local_announcement = "\tSetting up one more scattering experiment";

    new_state = mk_state();       /* in  statealgebra.c . The new state does */
                                  /* not yet contain a System or Regsystem.  */
    System(new_state) = setup_one_scene(nbody, answers, mproj, v_inf, theta,
					phi, psi, t_enc, rho, r_in, mfrac, a,
					e, inclination, longitude, argument,
					meananomaly);
    set_specs(Specs(new_state));
    set_controls(Ctrls(new_state), Tnow(System(new_state)));

    Announcement(Ctrls(new_state)) = local_announcement;

    return(new_state);
    }

/*-----------------------------------------------------------------------------
 *  perform_one_experiment  --  performs a scattering experiment
 *-----------------------------------------------------------------------------
 */
local void  perform_one_experiment(some_state, first_run)
stateptr  some_state;
bool  first_run;
    {
    rev_engine(some_state);                 /* determination of total energy;*/
                                            /* resides in file  orbit.c .    */
    write_initial_output(some_state, first_run);
                                            /* resides in file  out.c .      */
    evolve(some_state);                     /* resides in file  orbit.c .    */
    }

/*-----------------------------------------------------------------------------
 *  report_one_experiment  --  reports the outcome of a scattering experiment
 *-----------------------------------------------------------------------------
 */
local void  report_one_experiment(old_state, first_run, last_run, report)
stateptr  old_state;
bool  first_run;
bool  last_run;
multiclassptr  report;
    {
    int  i;
    char *hier_string;
    char  error_string[BUFF_LENGTH];
    diagptr  diags;
    bool  too_much_energy_drift();

    write_final_output(old_state, first_run, last_run);    /* in file  out.c */

    diags = Diags(old_state);

    if (too_much_energy_drift(diags, Ctrls(old_state),
			      DIAGetot(diags)[CURRENT]) == TRUE)
	{
	sprintf(error_string, "TOO MUCH ENERGY DRIFT");
	hier_string = error_string;
	}
    else
	hier_string = Hier_string(diags);

    i = 0;
    while (i < MAX_DIVERSITY)
	{
	if (Scattype(report)[i][0] == NULL)
	    {
	    sprintf(Scattype(report)[i], "%s", hier_string);
	    Scatsubtotal(report)[i] = 1;
	    Nscattype(report) += 1;
	    break;
	    }	    
	else if (streq(Scattype(report)[i], hier_string))
	    {
	    Scatsubtotal(report)[i] += 1;
	    break;
	    }
	i++;
	}
    if (i >= MAX_DIVERSITY)
	error("report_one_experiment: too many types of hier_string\n");
    }

/*-----------------------------------------------------------------------------
 *  clean_up_after_one_experiment  --  free up memory, after one many-body
 *                                     scattering experiment, within a series.
 *                             accept: old_state: a state pointer.
 *-----------------------------------------------------------------------------
 */
local void  clean_up_after_one_experiment(old_state)
stateptr  old_state;
    {
    rm_state(old_state);
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
    string  headline_default =
	"Multiscatter code: shared time steps, turbo version";

    if (turbo == FALSE)                            /* override defined value */
	sprintf(headline_default, "Multiscatter code: shared time steps");
#  endif
#endif
#ifdef REGULARIZATION
    string  headline_default =
            "Multiscatter code: equal time steps & regularization";
#endif
#ifdef EXTRAPOLATION
    string  headline_default = 
      "Multiscatter code: individual time steps: polynom. orbit extrapolation";
#endif
/*
 * determine the temporal relation between begin and end point of integration:
 *   Note: the control variable "Forwards()" is an essential ingredient
 *         in the temporal testing macros "earlier" and "later" (in newton0.h)
 */
    Tend(ctr) = getdparam("t_end");
    Forwards(ctr) = (Tend(ctr) > initial_time);
    if (! Forwards(ctr))
	error("multiscatter: only scattering forwards in time\n");

/*
 * initialize the other control variables which can directly receive outside
 * information through the command line arguments:
 */
    Lastoutfile(ctr) = getparam("lastout");
    Outfile(ctr) = getparam("out");
    DTminor(ctr) = getdparam("dt_minor");
    DTmajor(ctr) = getdparam("dt_major");
    DTsave(ctr) = 0.0;                           /* no savings */
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

#define  NPBODY_MAX   2*NBODY_MAX-1     /* maximum number of pseudoparticles */
                                        /* which can be read in for a single */
                                        /* scattering experiment             */
#define  UNKNOWN   1                    /* these are the possible      */
#define  SINGLE    2                    /* contents of the elements    */
#define  BINARY    3                    /* pseudo_body_table[...][0]   */

#define  ROOT_ADDRESS          0    /* these are the offsets of the root,    */
#define  TARGET_ADDRESS        1    /* target and projectile in the array:   */
#define  PROJECTILE_ADDRESS    2    /* pseudo_body_table[some_address][...]  */

/*-----------------------------------------------------------------------------
 *  setup_one_scene  --  initializes the masses, time and orbital parameters
 *                       for one many-body scattering experiment, within a
 *                       series
 *                  flow:
 *                       first the orbit of the center of mass of a projectile
 *                       is read in, and the clock is set at a desired time
 *                       before closest encounter of the projectile and the
 *                       target.
 *                       note: i) first the time of projected encounter is
 *                                read in (i.e. the time at which the 
 *                                encounter would take place if the projectile
 *                                would keep approaching the target in a 
 *                                perfect two-body hyperbola). This time is
 *                                taken with respect to the time t=0 defined
 *                                by the time of pericenter passage of the 
 *                                internal binary in the target.
 *                                The motivation for this choice is to
 *                                minimize the freedom in setting up the
 *                                target binary.
 *                            ii) then the time is transformed to a different
 *                                time frame in which the projected encounter
 *                                occurs at time t=0.
 *                                The motivation is to provide a more 
 *                                convenient time frame for performing the
 *                                scattering experiment.
 *                       then internal structure is added, first to the target
 *                       and then to the projectile; 
 *                  note: I: a SECOND choice has to be made to fix
 *                           the time scale, namely the time at which
 *                           the numerical orbit integration starts,
 *                           i.e. the time at which the initial
 *                           conditions should be given for all the
 *                           particles involved in the scattering
 *                           event. This has nothing to do with 
 *                           spacetime symmetries, but would occur
 *                           for every simulation with a less-then-
 *                           infinite duration. It has a purely
 *                           practical, and APPROXIMATE nature,
 *                           since it determines the point in time
 *                           before which the interactions between
 *                           target and projectiles are neglected.
 *                           As a consequence, chosing a different
 *                           starting time will lead to (slightly)
 *                           different results -- it is therefore up
 *                           to the experimenter to chose the initial
 *                           time judiciously.

 *                       II: The zero point of our time scale, defined
 *                           as the time of projected closest
 *                           encounter, differs slightly from the
 *                           actual time of closest encounter because
 *                           of the interactions between the internal
 *                           degrees of freedom of the target and 
 *                           projectiles -- besides, the actual time 
 *                           of closest encounter is likely to be
 *                           ill-defined when the internal structures
 *                           are taken into account. However, in the
 *                           limit of a computation starting
 *                           infinitely early, our choice of temporal
 *                           zero point is well defined. Of course,
 *                           as discussed under point I, truncation
 *                           of the computed part of the orbits to a
 *                           finite time interval introduces (slight)
 *                           perturbations in the exact definition of
 *                           time zero.
 *        implementation: 
 *                       A pseudo body is defined as either a single star or
 *                       the center of mass of a binary star (members of which
 *                       may or may not be composite themselves). First an
 *                       array of pseudo bodies is constructed, starting with
 *                       the root (the unbound "binary" consisting of target
 *                       and  projectile), the target and the projectile, in 
 *                       that order. The members of each binary present are
 *                       then added subsequently to the pseudo body array.
 *                       Finally the single stars are copied to the array
 *                       of real bodies.
 *            body table:
 *                       The body table has three entries for each body.
 *                       If a body has an array number  address , i.e. if it
 *                       is stored in  pseudo_bodies[address] , then:
 *                       pseudo_body_table[address][0] contains the character
 *                       of the body, either UNKNOWN, or SINGLE, or BINARY.
 *                       In the latter case, pseudo_body_table[address][1]
 *                       contains the address of the first member of the
 *                       binary, while pseudo_body_table[address][2] contains
 *                       the address of the second member of the binary.
 *-----------------------------------------------------------------------------
 */

local systptr  setup_one_scene(nbody, answers, mproj, v_inf, theta, phi, psi,
			       t_enc, rho, r_in, mfrac, a, e, inclination,
			       longitude, argument, meananomaly)
int  nbody;
bool  answers[2*NBODY_MAX-2];        /* information about binary hierarchies */
real  mproj;                         /* projectile mass                      */
real  v_inf;                     /* target-projectile rel. vel. at infinity  */
real  theta;
real  phi;
real  psi;
real  t_enc;
real  rho;
real  r_in;
real  mfrac[NBODY_MAX-2];        /* mass fraction of secondaries in binaries */
real  a[NBODY_MAX-2];            /* semimajor axes of binaries               */
real  e[NBODY_MAX-2];            /* eccentricities of binaries               */
real  inclination[NBODY_MAX-3];
real  longitude[NBODY_MAX-3];    /* of the ascending node                    */
real  argument[NBODY_MAX-3];     /* of pericenter                            */
real  meananomaly[NBODY_MAX-3];
    {
    int  n_pseudobody;
    int  n_realbody;
    int  answer_counter;
    int  a_counter;                   /* counts also the four angles         */
    int  e_counter;                   /* counts also the fractional masses;  */
                                      /* = a_counter + 1,  unless target is  */
                                      /* single and projectile is composite, */
                                      /* in which (unnatural!) case          */
                                      /* e_counter = a_counter               */
    int  pseudo_body_table[NPBODY_MAX][3];
    real  time_before_projected_encounter; /* when scattering expt. starts   */
    body  pseudo_bodies[NPBODY_MAX];       /* pointer to pseudobody array    */
    systptr  new_system;                   /* to be initialized              */
    systptr  mk_system();

    if (NDIM != 3)
	error("setup_one_scene: NDIM = %d != 3\n", NDIM);

    n_realbody = n_pseudobody = 0;                 /* start with empty space */
    answer_counter = a_counter = e_counter = 0;    /* and counters cleared   */

    setup_empty_root(&n_pseudobody, pseudo_body_table);

    setup_target_external_orbit(pseudo_bodies, &n_pseudobody, &n_realbody,
				pseudo_body_table);

    setup_projectile_external_orbit(pseudo_bodies, &n_pseudobody, &n_realbody,
                                    &time_before_projected_encounter, 
                                    pseudo_body_table, mproj, v_inf, theta,
				    phi, psi, t_enc, rho, r_in);

    setup_root_orbit(pseudo_bodies, n_pseudobody);

    transform_to_com_coordinates(pseudo_bodies, n_pseudobody);

    setup_target_internal_structure(pseudo_bodies, &n_pseudobody, &n_realbody,
				    answers, &answer_counter, &a_counter,
				    &e_counter, t_enc,
				    time_before_projected_encounter, 
                                    pseudo_body_table, mfrac, a, e,
				    inclination, longitude, argument,
				    meananomaly);

    setup_projectile_internal_structure(pseudo_bodies, &n_pseudobody,
					&n_realbody, answers, &answer_counter,
					&a_counter, &e_counter,
                                        time_before_projected_encounter,
					pseudo_body_table, mfrac, a, e,
					inclination, longitude, argument,
					meananomaly);
/*
 * transform from {z,x,y} coordinates, as provided by transkepler(),
 * to {x,y,z} coordinates:
 */
    kep2norm(pseudo_bodies, n_pseudobody);
/*
 * finally, package the results:
 */
    new_system = mk_system(n_realbody);         /* in file  systemalgebra.c  */
    Tnow(new_system) = - time_before_projected_encounter;
    Nbody(new_system) = n_realbody;
    make_body_array(Bodies(new_system), pseudo_bodies, n_realbody,
		    n_pseudobody, pseudo_body_table);
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
local void  setup_target_external_orbit(pseudo_bodies, n_pseudobody,
					n_realbody, pseudo_body_table)
bodyptr  pseudo_bodies;
int  *n_pseudobody;
int  *n_realbody;
int  pseudo_body_table[NPBODY_MAX][3];
    {
    real  orbital_elements[2*NDIM];

    if (*n_pseudobody != 1)
	error("setup_target_external_orbit: *n_pseudobody = %d != 1\n",
	      *n_pseudobody);	
    if (*n_realbody != 0)
	error("setup_target_external_orbit: *n_realbody = %d != 0\n",
	      *n_realbody);

    orbital_elements[0] = TARGET_COM_RX;
    orbital_elements[1] = TARGET_COM_RY;
    orbital_elements[2] = TARGET_COM_RZ;
    orbital_elements[3] = TARGET_COM_VX;
    orbital_elements[4] = TARGET_COM_VY;
    orbital_elements[5] = TARGET_COM_VZ;

    setup_body_orbit(pseudo_bodies + *n_pseudobody,
		     TARGET_MASS, orbital_elements);

    *n_pseudobody +=1;
    *n_realbody +=1;

    pseudo_body_table[TARGET_ADDRESS][0] = UNKNOWN;
    }

/*-----------------------------------------------------------------------------
 *  setup_projectile_external_orbit  --  initializes the mass and orbital
 *                                       parameters for the motion of the
 *                                       center of mass of a the projectile.
 *-----------------------------------------------------------------------------
 */
local void  setup_projectile_external_orbit(pseudo_bodies, n_pseudobody,
					    n_realbody, 
                                            time_before_projected_encounter,
					    pseudo_body_table, mproj, v_inf,
					    theta, phi, psi, t_enc, rho, r_in)
bodyptr  pseudo_bodies;
int  *n_pseudobody;
int  *n_realbody;
realptr  time_before_projected_encounter;
int  pseudo_body_table[NPBODY_MAX][3];
real  mproj;                     /* projectile mass                          */
real  v_inf;                     /* target-projectile rel. vel. at infinity  */
real  theta;
real  phi;
real  psi;
real  t_enc;
real  rho;
real  r_in;
    {
    real  total_mass;
    real  orbital_elements[2*NDIM];

    if (*n_pseudobody != 2)
	error("setup_projectile_external_orbit: *n_pseudobody = %d != 2\n",
	      *n_pseudobody);	
    if (*n_realbody != 1)
	error("setup_target_external_orbit: *n_realbody = %d != 1\n",
	      *n_realbody);

    if (mproj < 0.0)
	error("setup_projectile_external_orbit: negative mass\n");	

    orbital_elements[0] = rho;
    if (orbital_elements[0] < 0.0)
	error("setup_projectile_external_orbit: negative impact parameter\n");

    orbital_elements[1] = psi;

    orbital_elements[2] = r_in;
    if (orbital_elements[2] <= 0.0)
	error("setup_projectile_external_orbit: initial separation <= 0.0\n");

    orbital_elements[3] = v_inf;
    if (orbital_elements[3] < 0.0)
	error("setup_projectile_external_orbit: velocity at infinity < 0.0\n");

    orbital_elements[4] = theta;

    orbital_elements[5] = phi;

    total_mass = mproj + TARGET_MASS;

    r_in = orbital_elements[2];
    transkepler(total_mass, orbital_elements, SCATTERING, TSCATTERING);
    *time_before_projected_encounter = orbital_elements[2];
    orbital_elements[2] = r_in;

    transkepler(total_mass, orbital_elements, SCATTERING, CARTESIAN);

    setup_body_orbit(pseudo_bodies + *n_pseudobody, mproj, orbital_elements);

    *n_pseudobody +=1;
    *n_realbody +=1;

    pseudo_body_table[PROJECTILE_ADDRESS][0] = UNKNOWN;
    }

/*-----------------------------------------------------------------------------
 *  setup_target_internal_structure  --  reads in optional binary structure
 *                                       of the target object
 *-----------------------------------------------------------------------------
 */
local void  setup_target_internal_structure(pseudo_bodies, n_pseudobody,
					    n_realbody, answers,
					    answer_counter, a_counter,
					    e_counter,
					    time_of_projected_encounter,
					    time_before_projected_encounter,
					    pseudo_body_table, mfrac, a, e,
					    inclination, longitude, argument,
					    meananomaly)
bodyptr  pseudo_bodies;
int  *n_pseudobody;
int  *n_realbody;
bool  answers[2*NBODY_MAX-2];        /* information about binary hierarchies */
int *answer_counter;
int *a_counter;                      /* counts also the four angles          */
int *e_counter;                      /* counts also the fractional masses;   */
                                     /* = *a_counter + 1,  unless target is  */
                                     /* single and projectile is composite,  */
                                     /* in which (unnatural!) case           */
                                     /*  *e_counter = *a_counter             */
real  time_of_projected_encounter;
real  time_before_projected_encounter;
int  pseudo_body_table[NPBODY_MAX][3];
real  mfrac[NBODY_MAX-2];        /* mass fraction of secondaries in binaries */
real  a[NBODY_MAX-2];                /* semimajor axes of binaries           */
real  e[NBODY_MAX-2];                /* eccentricities of binaries           */
real  inclination[NBODY_MAX-3];
real  longitude[NBODY_MAX-3];        /* of the ascending node                */
real  argument[NBODY_MAX-3];         /* of pericenter                        */
real  meananomaly[NBODY_MAX-3];
    {
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
    if (*n_realbody != 2)
	error("setup_target_internal_structure: *n_realbody = %d != 2\n",
	      *n_realbody);

    if (pseudo_body_table[TARGET_ADDRESS][0] != UNKNOWN)
	error("setup_target_internal_structure: TARGET not of type UNKNOWN\n");

    if (answers[(*answer_counter)++] == FALSE)               /* single star? */
	{
	pseudo_body_table[TARGET_ADDRESS][0] = SINGLE;
	return;            /* then return, leaving the body count unchanged; */
        }                  /* otherwise, initialize binary components:       */

    target = pseudo_bodies + TARGET_ADDRESS;

    secondary_mass_fraction = mfrac[*e_counter];    /* counter updated below */
    if (secondary_mass_fraction < 0.0)
	error("setup_target_internal_structure: secondary_mass_fraction <0\n");
    if (secondary_mass_fraction > 1.0)
	error("setup_target_internal_structure: secondary_mass_fraction >1\n");
    primary_mass = (1.0 - secondary_mass_fraction) * TARGET_MASS;
    secondary_mass = secondary_mass_fraction * TARGET_MASS;

    orbital_elements[0] = TARGET_SEMIMAJOR_AXIS;

    orbital_elements[1] = e[*e_counter];            /* counter updated below */
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

    setup_offset(&primary, &secondary, orbital_elements,
		 primary_mass, secondary_mass);

    add_to_com(&primary, target);
    add_to_com(&secondary, target);

    insert_body_orbit(pseudo_bodies + *n_pseudobody, &primary);
    *n_pseudobody +=1;

    insert_body_orbit(pseudo_bodies + *n_pseudobody, &secondary);
    *n_pseudobody +=1;

    *n_realbody +=1;

    primary_address = *n_pseudobody - 2;
    secondary_address = *n_pseudobody - 1;
    pseudo_body_table[TARGET_ADDRESS][0] = BINARY;
    pseudo_body_table[TARGET_ADDRESS][1] = primary_address;
    pseudo_body_table[TARGET_ADDRESS][2] = secondary_address;
    pseudo_body_table[primary_address][0] = UNKNOWN;
    pseudo_body_table[secondary_address][0] = UNKNOWN;

    if (debug)
	internal_structure_message(TARGET_ADDRESS, pseudo_body_table);

    (*e_counter)++;        /* only eccentricity and mass fraction used here; */
                           /* no semimajor axis and angles                   */

    setup_node_internal_structure(pseudo_bodies, primary_address,
				  n_pseudobody, n_realbody, answers,
				  answer_counter, a_counter, e_counter,
				  time_before_projected_encounter,
				  pseudo_body_table, mfrac, a, e, inclination,
				  longitude, argument, meananomaly);

    setup_node_internal_structure(pseudo_bodies, secondary_address,
				  n_pseudobody, n_realbody, answers,
				  answer_counter, a_counter, e_counter,
				  time_before_projected_encounter,
				  pseudo_body_table, mfrac, a, e, inclination,
				  longitude, argument, meananomaly);
    }

/*-----------------------------------------------------------------------------
 *  setup_projectile_internal_structure  --  reads in optional binary structure
 *                                           of the projectile
 *-----------------------------------------------------------------------------
 */
local void setup_projectile_internal_structure(pseudo_bodies, n_pseudobody,
					       n_realbody, answers,
					       answer_counter, a_counter,
					       e_counter,
					       time_before_projected_encounter,
					       pseudo_body_table, mfrac, a, e,
					       inclination, longitude,
					       argument, meananomaly)
bodyptr  pseudo_bodies;
int  *n_pseudobody;
int  *n_realbody;
bool  answers[2*NBODY_MAX-2];        /* information about binary hierarchies */
int *answer_counter;
int *a_counter;                      /* counts also the four angles          */
int *e_counter;                      /* counts also the fractional masses;   */
                                     /* = *a_counter + 1,  unless target is  */
                                     /* single and projectile is composite,  */
                                     /* in which (unnatural!) case           */
                                     /*  *e_counter = *a_counter             */
real  time_before_projected_encounter;
int  pseudo_body_table[NPBODY_MAX][3];
real  mfrac[NBODY_MAX-2];        /* mass fraction of secondaries in binaries */
real  a[NBODY_MAX-2];                /* semimajor axes of binaries           */
real  e[NBODY_MAX-2];                /* eccentricities of binaries           */
real  inclination[NBODY_MAX-3];
real  longitude[NBODY_MAX-3];        /* of the ascending node                */
real  argument[NBODY_MAX-3];         /* of pericenter                        */
real  meananomaly[NBODY_MAX-3];
    {
    setup_node_internal_structure(pseudo_bodies, PROJECTILE_ADDRESS,
				  n_pseudobody, n_realbody, answers,
				  answer_counter, a_counter, e_counter,
				  time_before_projected_encounter,
				  pseudo_body_table, mfrac, a, e, inclination,
				  longitude, argument, meananomaly);
    }

/*-----------------------------------------------------------------------------
 *  setup_node_internal_structure  --  
 *-----------------------------------------------------------------------------
 */
local void  setup_node_internal_structure(pseudo_bodies, object_address,
					  n_pseudobody, n_realbody, answers,
					  answer_counter, a_counter, e_counter,
					  time_before_projected_encounter,
					  pseudo_body_table, mfrac, a, e,
					  inclination, longitude, argument,
					  meananomaly)
bodyptr  pseudo_bodies;
int  object_address;
int  *n_pseudobody;
int  *n_realbody;
bool  answers[2*NBODY_MAX-2];        /* information about binary hierarchies */
int *answer_counter;
int *a_counter;                      /* counts also the four angles          */
int *e_counter;                      /* counts also the fractional masses;   */
                                     /* = *a_counter + 1,  unless target is  */
                                     /* single and projectile is composite,  */
                                     /* in which (unnatural!) case           */
                                     /*  *e_counter = *a_counter             */
real  time_before_projected_encounter;
int  pseudo_body_table[NPBODY_MAX][3];
real  mfrac[NBODY_MAX-2];        /* mass fraction of secondaries in binaries */
real  a[NBODY_MAX-2];                /* semimajor axes of binaries           */
real  e[NBODY_MAX-2];                /* eccentricities of binaries           */
real  inclination[NBODY_MAX-3];
real  longitude[NBODY_MAX-3];        /* of the ascending node                */
real  argument[NBODY_MAX-3];         /* of pericenter                        */
real  meananomaly[NBODY_MAX-3];
    {
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
	error("setup_node_internal_structure: *n_pseudobody=%d>NPBODY_MAX-2\n",
	      *n_pseudobody);
    if (*n_realbody >= *n_pseudobody)
	error("setup_target_external_orbit:*n_realbody=%d>=*n_pseudobody=%d\n",
	      *n_realbody, *n_pseudobody);

    if (pseudo_body_table[object_address][0] != UNKNOWN)
	error("setup_node_internal_structure: object not of type UNKNOWN\n");

    if (answers[(*answer_counter)++] == FALSE)               /* single star? */
	{
	pseudo_body_table[object_address][0] = SINGLE;
	return;            /* then return, leaving the body count unchanged; */
        }                  /* otherwise, initialize binary components:       */

    object = pseudo_bodies + object_address;

    secondary_mass_fraction = mfrac[*e_counter];    /* counter updated below */
    if (secondary_mass_fraction < 0.0)
	error("setup_node_internal_structure: secondary_mass_fraction < 0\n");
    if (secondary_mass_fraction > 1.0)
	error("setup_node_internal_structure: secondary_mass_fraction > 1\n");
    primary_mass = (1.0 - secondary_mass_fraction) * Mass(object);
    secondary_mass = secondary_mass_fraction * Mass(object);

    orbital_elements[0] = a[*a_counter];            /* counter updated below */
    if (orbital_elements[0] < 0.0)
	error("setup_node_internal_structure: semimajor axis a < 0.0\n");

    orbital_elements[1] = e[*e_counter];            /* counter updated below */
    if (orbital_elements[1] < 0.0)
	error("setup_node_internal_structure: eccentricity < 0.0\n");
    if (orbital_elements[1] > 1.0)
	error("setup_node_internal_structure: eccentricity > 1.0\n");

    orbital_elements[2] = inclination[*a_counter];  /* counter updated below */
    orbital_elements[3] = longitude[*a_counter];
    orbital_elements[4] = argument[*a_counter];
    orbital_elements[5] = meananomaly[*a_counter];

    transkepler(Mass(object), orbital_elements, MEANANOMALY, PERIPASSAGE);

    orbital_elements[5] += time_before_projected_encounter;
                /* in time frame in which start of simulation is at time 0.0 */

    transkepler(Mass(object), orbital_elements, PERIPASSAGE, CARTESIAN);

    setup_offset(&primary, &secondary, orbital_elements,
		 primary_mass, secondary_mass);

    add_to_com(&primary, object);
    add_to_com(&secondary, object);

    insert_body_orbit(pseudo_bodies + *n_pseudobody, &primary);
    *n_pseudobody +=1;

    insert_body_orbit(pseudo_bodies + *n_pseudobody, &secondary);
    *n_pseudobody +=1;

    *n_realbody +=1;

    primary_address = *n_pseudobody - 2;
    secondary_address = *n_pseudobody - 1;
    pseudo_body_table[object_address][0] = BINARY;
    pseudo_body_table[object_address][1] = primary_address;
    pseudo_body_table[object_address][2] = secondary_address;
    pseudo_body_table[primary_address][0] = UNKNOWN;
    pseudo_body_table[secondary_address][0] = UNKNOWN;

    if (debug)
	internal_structure_message(object_address, pseudo_body_table);

    (*e_counter)++;
    (*a_counter)++;

    setup_node_internal_structure(pseudo_bodies, primary_address,
				  n_pseudobody, n_realbody, answers,
				  answer_counter, a_counter, e_counter,
				  time_before_projected_encounter,
				  pseudo_body_table, mfrac, a, e, inclination,
				  longitude, argument, meananomaly);

    setup_node_internal_structure(pseudo_bodies, secondary_address,
				  n_pseudobody, n_realbody, answers,
				  answer_counter, a_counter, e_counter,
				  time_before_projected_encounter,
				  pseudo_body_table, mfrac, a, e, inclination,
				  longitude, argument, meananomaly);
    }

/*-----------------------------------------------------------------------------
 *  make_body_array  --  from among the array of pseudo bodies, copy the real
 *                       bodies to the body array. The numbering is in the
 *                       order of the body array, except that first all the
 *                       target bodies are listed, and then all the projectile
 *                       bodies.
 *-----------------------------------------------------------------------------
 */
local void  make_body_array(some_bodies, pseudo_bodies, n_realbody,
			    n_pseudobody, pseudo_body_table)
bodyptr  some_bodies;
bodyptr  pseudo_bodies;
int  n_realbody;
int  n_pseudobody;
int  pseudo_body_table[NPBODY_MAX][3];
    {
    int  i;                                           /* body address        */
    int  j;                                           /* pseudo body address */

    if (n_pseudobody < 3)
	error("make_body_array: n_pseudobody = %d < 3\n", n_pseudobody);
    if (n_realbody >= n_pseudobody)
	error("make_body_array: n_realbody= %d >= n_pseudobody= %d\n",
	      n_realbody, n_pseudobody);

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

    if (i != n_realbody)
	error("make_body_array: i = %d != n_realbody = %d\n", i, n_realbody);
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
	printf("\n=::=*=::=*=::=*=::=*=::=*=::=*=::=*=::=");
	printf("*=::=*=::=*=::=*=::=*=::=*=::=*=::=*=::=\n\n");
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

/* endof: multiscatter.c */
