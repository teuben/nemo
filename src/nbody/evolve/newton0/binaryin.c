/* binaryin.c - NDIM_IO, SETV_IN, get_mass, get_parameters, get_phase,
                get_particles, get_snapshot, introduce_system */

/*
 *  binaryin.c: input module for newton0.c : equal time steps
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 *	July 1990  - new location of <snapshot.h>macro	PJT
 *      20 feb 2004 - skip history
 */
   
#include  "newton0.h"
#include  <filestruct.h>
#include  <archaic/snapshot.h>

#ifndef REGULARIZATION
#  define  NDIM_IO  NDIM                /* dimensionality for input/output   */
#  define  SETV_IN  SETV                /* vector arithmetic for input       */
#endif

#ifdef REGULARIZATION

#  define  NDIM_IO  3              /* dimensionality for input/output:       */
                                   /* external data are written in 3D-format */
                                   /* while internal data are all in 4D      */

#  define SETV_IN(v,u)		/* SET 4D Vector using a 3D representation */ \
       { 				/* for vectors read from input     */ \
       register int _i; 					              \
       for (_i = 0; _i < NDIM_IO; _i++) 				      \
           (v)[_i] = (u)[_i]; 						      \
       (v)[NDIM_IO] = 0;           /* fourth component vanishes identically */\
       }
#endif

local int  coordsys = CSCode(Cartesian, NDIM_IO, 2);
local stream  instr;                        	     /* input stream pointer */

static systptr get_snapshot(void);
static void get_parameters(int *nbodyptr, real *timeptr);
static bodyptr get_particles(int nbody);
static void get_mass(bodyptr bodies, int nbody);
static void get_phase(bodyptr bodies, int nbody);

extern   bodyptr  mk_bodies();
extern   systptr  mk_empty_system();


/*    char *malloc();             /* to allocate memory for the masses */


/*-----------------------------------------------------------------------------
 *  introduce_system  --  get the initial data from the input file
 *-----------------------------------------------------------------------------
 */

systptr  introduce_system(infile, headlineptr, announceptr)
string  infile;
string  *headlineptr;
string  *announceptr;
    {
    systptr  new_system;
    string local_announcement = "\tReading an old Snapshot";

    instr = stropen(infile, "r");		    /* open input stream */
    if (get_tag_ok(instr, HeadlineTag))		    /* read headline, if any */
	*headlineptr = get_string(instr, HeadlineTag);
    else
	*headlineptr = NULL;
    *announceptr = local_announcement;
    new_system = get_snapshot();		    /* read snapshot data */
    strclose(instr);				    /* close input stream */

    return(new_system);    
    }

/*-----------------------------------------------------------------------------
 *  get_snapshot  --  read in the first snapshot with the right HeadlineTag,
 *                    as selected above in  introduce_system() .
 *-----------------------------------------------------------------------------
 */
local systptr  get_snapshot()
    {
    systptr  new_system;

    new_system = mk_empty_system();             /* in file  systemalgebra.c  */

    get_history(instr);
    get_set(instr, SnapShotTag);		/* access snapshot data      */
    get_parameters(&Nbody(new_system), &Tnow(new_system));
    Bodies(new_system) = get_particles(Nbody(new_system));
    get_tes(instr, SnapShotTag);		/* finish reading snapshot   */

    return(new_system);    
    }

/*-----------------------------------------------------------------------------
 *  get_parameters  --  read in the snapshot parameters
 *-----------------------------------------------------------------------------
 */
local void  get_parameters(nbodyptr, timeptr)
int  *nbodyptr;
real  *timeptr;
    {
    get_set(instr, ParametersTag);		/* access parameters */
    if (! get_tag_ok(instr, TimeTag))		/* no time data given? */
        *timeptr = 0.0;		                /*   start at 0.0 */
    else
	get_data(instr, TimeTag, RealType, timeptr, 0);
    get_data(instr, NobjTag, IntType, nbodyptr, 0);
    if (*nbodyptr < 1)              		/* silly value for nbody? */
	error("get_parameters: %s = %d out of range\n", NobjTag, *nbodyptr);
    get_tes(instr, ParametersTag);		/* finished parameters */
    }

/*-----------------------------------------------------------------------------
 *  get_particles  --  read in the snapshot particle data
 *-----------------------------------------------------------------------------
 */
local bodyptr  get_particles(nbody)
int  nbody;
    {
    int cs;
    bodyptr  bodies;

    bodies = mk_bodies(nbody);

    get_set(instr, ParticlesTag);		/* access particles */
    get_data(instr, CoordSystemTag, IntType, &cs, 0);
    if (cs != coordsys)				/* not standard coords? */
	error("get_particles: cant handle %s = 0%o\n", CoordSystemTag, cs);
    get_mass(bodies, nbody);
    get_phase(bodies, nbody);
    get_tes(instr, ParticlesTag);		/* finished particles */

    return(bodies);    
    }

/*-----------------------------------------------------------------------------
 *  get_mass  --  for the particles in the snapshot, read in their mass, ...
 *-----------------------------------------------------------------------------
 */
local void  get_mass(bodies, nbody)
bodyptr  bodies;
int  nbody;
    {
    real *tmpmass, *massptr;
    bodyptr  body_i;

    tmpmass = (real *) malloc((unsigned)nbody*sizeof(real));
    if (tmpmass == NULL)
	error("get_mass: not enough memory left for tmpmass[%d]\n", nbody);

    get_data(instr, MassTag, RealType, tmpmass, nbody, 0);
    massptr = tmpmass;    			/* point to mass array */
    for (body_i = bodies; body_i - bodies < nbody; body_i++)
	Mass(body_i) = *massptr++;		/*   copy individual masses */
    free(tmpmass);
    }

/*-----------------------------------------------------------------------------
 *  get_phase  --  ... , and their positions and velocities
 *-----------------------------------------------------------------------------
 */
local void  get_phase(bodies, nbody)
bodyptr  bodies;
int  nbody;
    {
    real *tmpphase, *pp;
    bodyptr  body_i;

    tmpphase = (real *) malloc((unsigned)2*NDIM_IO*nbody*sizeof(real));
    if (tmpphase == NULL)
	error("get_phase: not enough memory left for tmpphase[%d]\n",
                                                              2*NDIM_IO*nbody);

    get_data(instr, PhaseSpaceTag, RealType, tmpphase, nbody, 2, NDIM_IO, 0);
    pp = tmpphase;
    for (body_i = bodies; body_i - bodies < nbody; body_i++)
        {
	SETV_IN(Pos(body_i), pp);		/* copy position data     */
	pp += NDIM_IO;
	SETV_IN(Vel(body_i), pp);		/* copy velocity data     */
	pp += NDIM_IO;
        }
    free (tmpphase);
    }

/* endof: binaryin.c */
