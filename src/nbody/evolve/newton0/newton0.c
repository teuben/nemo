/* newton0.c - BIGNUMBER, CORNEROFFSET, create_and_set, create_system, finish, 
               introduce_and_set, main, newton0, read_and_set, resume_and_set, 
               set_controls, set_specs, set_system, setup, start */

/*
 *  newton0.c:  simple gravitational many-body code: equal time steps
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 *      Feb  2001  - resurrected for the current NEMO release
 *      Feb  2004  - added to NEMO's official release
 */
   
#include  "newton0.h"
#include  <getparam.h>

static void newton0(void);
static stateptr start(void);
static void finish(stateptr old_state);
static stateptr read_and_set(void);
static stateptr resume_and_set(string resumefile);
static stateptr introduce_and_set(string infile);
static stateptr create_and_set(void);
static systptr create_system(string *announceptr);
static void setup(stateptr new_state);
static void set_system(systptr sys);
static void set_specs(specptr specs);
static void set_controls(ctrlptr ctr, real initial_time);


/*-----------------------------------------------------------------------------
 *  newton0: every collective timestep every force from each particle on each
 *           other particle is computed, and integrated according to the 
 *           integration scheme specified.
 *  newton0tree: as newton0, but every particle feels only the forces from a
 *               number of selected nodes in an Eulerian tree; a node is a 
 *               either an individual particle or a cell which represents a
 *               cluster of particles.
 *  newton0reg: as newton0, but using four-dimensional regularization for each
 *              particle pair simultaneously.
 *
 *  newton0 ,  newton0tree  and  newton0reg  each can be invoked without any
 *  command line arguments: in that case no binary output is provided but
 *  only output on the standard outstream (by default on the screen), and
 *  a Plummer model is created from scratch starting with a random number
 *  taken from the UNIX clock, with the effect that two different commands
 *  newton0 will give different results unless given within one second
 *  interval.
 *_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  _ _ _ _ _ _ _ 
 *
 *  The array of strings  defv[]  given on the next page contains the default
 *  values for all command line arguments (default values may be empty).  Using
 *   getparam()  and related procedures, these values are assigned to the
 *  appropriate variables in the following procedures below:
 *
 *  read_and_set(), create_system(), set_controls(), set_specs(), set_system().
 *
 *-----------------------------------------------------------------------------
 */

string  defv[] =            /* DEFAULT INPUT PARAMETERS                      */
    {           
    "in=\n	    	    input file name",
    "lastout=\n             file name for one single last output",
    "out=\n		    output file name",
    "save=\n                save file, in case restore is needed",
    "resume=\n              resume file from previous save, with new control pars",
    "restore=\n             restore file , no new pars",
    "restart=\n             file in which to write the restart command",
    "nbody=3\n  	    number of particles",
    "seed=0\n               seed for the random number generator",
    "t_begin=0.0\n           time at which the integration starts",
    "t_end=1.0\n             time at which the integration should halt",
    "dt_minor=0.5\n          time interval between minor outputs",
    "dt_major=1.0\n          time interval between major outputs",
    "dt_save=0.5\n           time interval between system state saving",
    "dt_max=1000.0\n         maximum size of integration time step",
    "nstep=10000000\n        maximum allowed number of integration steps",
    "de_max=1.0e7\n          maximum allowed drift in total energy",
    "nstep_de=100000000\n    number of steps after which  de_max  checked",
    "dn_major=100000000\n    number of integration steps between major outputs",
    "min_pairdist=0.0\n      minimum allowed pair distance (softened)",
    "cpu_max=1.0e7\n         maximum allowed amount of CPU time (minutes)",
    "cpu_last_call=1.0e7\n   amount of CPU time after which integration halted at next major-output",
    "eta_acc=0.01\n          dimensionless integration accuracy parameter (for time step)",
    "r0_soft=0.0\n           potential softening length",
    "soft_focus=plummer\n    type of softening used",
    "timestep=constant\n     timestep criterion",
    "diagnostics=standard\n  set of diagnostics provided at output times",
#ifndef TREE
    "scheme=runge_kutta_4\n  integration scheme",
#endif
#ifdef TREE
    "celldivision=constant_theta\n       cell subdivision criterion",
    "scheme=leapfrog\n       integration scheme",
    "tol=1.0\n               cell subdivision tolerance",
#endif
#ifdef REGULARIZATION
    "niter=3\n               number of iterations in arriving at the correct output times",
#endif
    "headline=\n             verbiage for output",
    "VERSION=2.0a\n          20-feb-04 PJT",
    NULL,
    };


/*-----------------------------------------------------------------------------
 *  main  --  read the command line arguments, and pass control to newton0
 *  nemo_main in NEMO V3
 *-----------------------------------------------------------------------------
 */
void nemo_main()
    {
        newton0();
    }

/*-----------------------------------------------------------------------------
 *  newton0  --  currently orchestrates only one gravitational evolution
 *-----------------------------------------------------------------------------
 */
local void  newton0()
    {
    stateptr  the_state;

    the_state = start();

    evolve(the_state);                         /* in file  orbit.c   */

    finish(the_state);
    }

/*-----------------------------------------------------------------------------
 *  start  --  initializes state variables, and prepares for orbit integration
 *-----------------------------------------------------------------------------
 */
local stateptr  start()
    {
    stateptr  new_state;

    new_state = read_and_set();        /* manages all types of input.        */
                                            /* rev the engines, to enable    */
    rev_engine(new_state);                  /* determination of total energy;*/
                                            /* resides in file  orbit.c .    */
    write_initial_output(new_state);   /* resides in file  out.c .           */

    return(new_state);
    }

/*-----------------------------------------------------------------------------
 *  finish  --  orchestrates final diagnostics & output, and closes all files
 *-----------------------------------------------------------------------------
 */
local void  finish(old_state)
stateptr  old_state;
    {
    write_final_output(old_state);         /* in file  out.c */
    }

/*-----------------------------------------------------------------------------
 *  read_and_set  --  supervises all input:
 *                    read refers to reading binary input; 
 *                      binary input files can contain either a snapshot or a
 *                      previously saved system state;
 *                    set refers to setting state variables;
 *                      these will be set to values provided as command line
 *                      arguments, or to default values, or to fixed 
 *                      values, depending on the individual variables.
 *-----------------------------------------------------------------------------
 */
local stateptr  read_and_set()
    {
    bool  in_flag;
    bool  restore_flag;
    bool  resume_flag;
    string  infile;       /* input file                                      */
    string  restorefile;  /* file containing data saved from a previous run, */
                          /* which will be continued without any change.     */
    string  resumefile;   /* file containing data saved from a previous run, */
                          /* which will be resumed with new control          */
                          /* parameters.                                     */
    stateptr  new_state;

/*
 * The four startup possibilities are: restore, resume, introduce, create:
 */
    restorefile = getparam("restore");
    restore_flag = hasvalue("restore");

    resumefile = getparam("resume");
    resume_flag = hasvalue("resume");

    infile = getparam("in");
    in_flag = hasvalue("in");

/*
 * at most one of the choices restore, resume and introduce can be followed:
 */
    if (restore_flag && resume_flag && in_flag)
	error("read_and_set: choose one among restore, resume or input\n");
    else if (restore_flag && resume_flag)
	error("read_and_set: choose either restore or resume\n");
    else if (in_flag && restore_flag)
	error("read_and_set: choose either input or restore\n");
    else if (in_flag && resume_flag)
	error("read_and_set: choose either input or resume\n");

/*
 * carry out the startup of choice: 
 * the default choice is to create a new nbody system.
 */
    if (restore_flag)
	new_state = restore_state(restorefile);       /* in file  save.c  */
    else if (resume_flag)
	new_state = resume_and_set(resumefile);
    else if (in_flag)
	new_state = introduce_and_set(infile);
    else
	new_state = create_and_set();

    return(new_state);
    }

/*-----------------------------------------------------------------------------
 *  resume_and_set  --  set state variables, after restoring the state
 *                      of a previous run using the data in the file with the
 *                      name xxx in "resume=xxx" on the command line.
 *-----------------------------------------------------------------------------
 */
local stateptr  resume_and_set(resumefile)
string  resumefile;
    {
    stateptr  new_state;

    new_state = restore_state( resumefile );            /* in file  save.c  */

    setup(new_state);

    return(new_state);
    }

/*-----------------------------------------------------------------------------
 *  introduce_and_set  --  set state variables, after reading in the
 *                         positions, velocities and masses given in binary
 *                         form (in standard snapshot format) in the file
 *                         with the name yyy in "in=yyy" on the command line.
 *-----------------------------------------------------------------------------
 */
local stateptr  introduce_and_set(infile)
string  infile;
    {
    stateptr  new_state;

    new_state = mk_state();                           /* in  statealgebra.c  */

    System(new_state) = introduce_system( infile,          /* in  binaryin.c */
                                          &Headline(Ctrls(new_state)),
                                          &Announcement(Ctrls(new_state)) );
    setup(new_state);

    return(new_state);
    }

/*-----------------------------------------------------------------------------
 *  create_and_set  --  set state variables, and construct a new Plummer
 *                      model to provide initial conditions for all particles.
 *-----------------------------------------------------------------------------
 */
local stateptr  create_and_set()
    {
    stateptr  new_state;

    new_state = mk_state();                           /* in  statealgebra.c  */

    System(new_state) = create_system( &Announcement(Ctrls(new_state)) );

    setup(new_state);

    return(new_state);
    }

/*-----------------------------------------------------------------------------
 *  create_system  --  constructs a new system in the form of a Plummer model.
 *                     the necessary parameters are read using getparam() & co.
 *-----------------------------------------------------------------------------
 */
local systptr  create_system(announceptr)
string  *announceptr;
    {
    int  seed;
    int  npart;
    systptr  new_system;
    permanent char  local_announcement[128];

    seed = getiparam("seed");
    if (seed == 0)       /* no particular positive seed provided?            */
	seed = time(0);  /* then give a random value, different every second */
    srandom(seed);

    npart = getiparam("nbody");
    if (npart < 1)
	error("creation_state: nbody = %d, but nbody > 0 is required\n",npart);

    new_system = mk_empty_system();          /* in file  systemalgebra.c  */

    Tnow(new_system) = getdparam("t_begin");
    Nbody(new_system) = npart;
    Bodies(new_system) = create_plum(npart);       /* in file  create.c  */

    sprintf(local_announcement,
            "\tCreating a new Plummer model from seed = %d", seed);
    *announceptr = local_announcement;

    return(new_system);
    }

/*-----------------------------------------------------------------------------
 *  setup  --  set variables in diverse state substructures
 *             accept: new_state: a state pointer.
 *             note: the fourth state substructure, containing the diagnostics,
 *                   cannot yet be set at this point, since that would require
 *                   energy diagnostics, which requires revving the engines,
 *                   which could not have been done before the other three
*                    substructures are set.
 *-----------------------------------------------------------------------------
 */
local void  setup(new_state)
stateptr  new_state;
    {
    set_system(System(new_state));
    set_specs(Specs(new_state));
    set_controls(Ctrls(new_state), Tnow(System(new_state)));
    }

/*-----------------------------------------------------------------------------
 *  set_system  --  set system variables
 *                  accept: sys: a system pointer.
 *-----------------------------------------------------------------------------
 */
#ifdef TREE
#  define  CORNEROFFSET  2.0
#endif

local void  set_system(sys)
systptr  sys;
    {
#ifdef TREE
    SETVS(Potmincorner(sys), -CORNEROFFSET);          /* initial box scaling */
    Potsizesq(sys) = (2.0 * CORNEROFFSET) * (2.0 * CORNEROFFSET);
#endif                                            /* initial square box size */
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
    real  tol;

    Softfocus(specs) = getparam("soft_focus");
    Timestepmethod(specs) = getparam("timestep");
    Integrationscheme(specs) = getparam("scheme");
    Diagnostics(specs) = getparam("diagnostics");

    Stepparam(specs) = getdparam("eta_acc");
    Softparam(specs) = getdparam("r0_soft");
#ifdef EXTRAPOLATION
    if (Softparam(specs) > 1.0/BIGNUMBER)
	error("set_specs: newton0_ext: not yet non-zero softening length\n");
#endif
#ifdef TREE
    Celldivisionmethod(specs) = getparam("celldivision");
    tol = getdparam("tol");
    Tolsqparam(specs) = tol * tol;
    Nfcalc(specs) = N2bcalc(specs) = Nbccalc(specs) = 0;
#endif
    Nstep_de(specs) = getiparam("nstep_de");
    SPECmin_pair(specs) = BIGNUMBER;
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
#ifndef TREE
#  ifndef REGULARIZATION
#    ifndef EXTRAPOLATION
    string  headline_default = "Newton0 code: equal time steps";
#    endif
#  endif
#endif
#ifdef TREE
#  ifndef QUADPOLE
    string  headline_default = 
            "Newton0 code: equal time steps & tree forces, monopoles only";
#   else
    string  headline_default = 
            "Newton0 code: equal time steps & tree forces, up to quadrupoles";
#  endif
#endif
#ifdef REGULARIZATION
    string  headline_default =
            "Newton0 code: equal time steps & regularization";
#endif
#ifdef EXTRAPOLATION
    string  headline_default = 
         "Newton0 code: individual time steps: polynomial orbit extrapolation";
#endif

/*
 * determine the temporal relation between begin and end point of integration:
 *   Note: the control variable "Forwards()" is an essential ingredient
 *         in the temporal testing macros "earlier" and "later" (in newton0.h)
 */
    Tend(ctr) = getdparam("t_end");
    Forwards(ctr) = (Tend(ctr) > initial_time);
                                        /* progressing or regressing? */

/*
 * initialize the other control variables which can directly receive outside
 * information through the command line arguments:
 */
    Infile(ctr)      = getparam("in");
    Lastoutfile(ctr) = getparam("lastout");
    Outfile(ctr)     = getparam("out");
    Savefile(ctr)    = getparam("save");
    Resumefile(ctr)  = getparam("resume");
    Restorefile(ctr) = getparam("restore");
    Restartfile(ctr) = getparam("restart");

    DTmajor(ctr) = getdparam("dt_major");
    DTminor(ctr) = getdparam("dt_minor");
    DTsave(ctr) = getdparam("dt_save");
    DTmax(ctr) = getdparam("dt_max");
    if (! Forwards(ctr))
        {                                          /* leaves the user the    */
	DTmajor(ctr) = -ABS(DTmajor(ctr));         /* freedom of sign choice */
	DTminor(ctr) = -ABS(DTminor(ctr));         /* while providing        */
	DTsave(ctr) = -ABS(DTsave(ctr));           /* the correct sign       */
        DTmax(ctr) = -ABS(DTmax(ctr));             /* in all four variables  */
        }
    Nmaxstep(ctr) = getiparam("nstep");
    DNmajor(ctr) = getdparam("dn_major");
#ifdef REGULARIZATION
    Ntimingiter(ctr) = getiparam("niter");
#endif
    DEmax(ctr) = getdparam("de_max");
    CPUmax(ctr) = getdparam("cpu_max");
    CPUlastcall(ctr) = getdparam("cpu_last_call");
    Min_pairdist(ctr) = getdparam("min_pairdist");

    headline_helper = getparam("headline");
    if (*headline_helper != 0)                   /* 3 possible headlines:    */
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

/* endof: newton0.c */
