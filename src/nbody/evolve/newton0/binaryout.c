/* binaryout.c - NDIM_IO, SETV_OUT, close_binary_output, coordsys,
                 open_binary_output, outstr, put_acc, put_diagnostics, 
                 put_mass, put_parameters, put_particles, put_phase, put_pot, 
                 put_full_snapshot, put_diagnostics_snapshot */

/*
 *  binaryout.c: for binary output for newton0 : equal time steps
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 *	July 1990  - new location of <snapshot.h>macro	PJT
 */
   
#include  "newton0.h"
#include  <filestruct.h>
#include  <archaic/snapshot.h>

#ifndef REGULARIZATION
#  define  NDIM_IO  NDIM                /* dimensionality for input/output   */
#  define  SETV_OUT  SETV               /* vector arithmetic for input       */
#endif

#ifdef REGULARIZATION

#  define  NDIM_IO  3              /* dimensionality for input/output:       */
                                   /* external data are written in 3D-format */
                                   /* while internal data are all in 4D      */

#  define SETV_OUT(v,u)		/* SET 3D Vector using a 4D representation */ \
       { 				/* for vectors written to output   */ \
       register int _i; 					              \
       for (_i = 0; _i < NDIM_IO; _i++) 				      \
           (v)[_i] = (u)[_i]; 						      \
       }
#endif

local int coordsys = CSCode(Cartesian, NDIM_IO, 2);
local stream outstr;                                /* output stream pointer */


static void put_parameters(int *nbodyptr, real *timeptr);
static void put_particles(systptr sys);
static void put_mass(bodyptr bodies, int nbody);
static void put_phase(bodyptr bodies, int nbody);
static void put_pot(bodyptr bodies, int nbody);
static void put_acc(bodyptr bodies, int nbody);

extern real  cputime();


/*-----------------------------------------------------------------------------
 *  open_binary_output  --  opens binary output stream;
 *                          directed by  write_initial_output()  in out.c
 *-----------------------------------------------------------------------------
 */
void   open_binary_output(outfile, headline)
string  outfile;
string  headline;
    {
    outstr = stropen(outfile, "w");                   /* setup output stream */
    put_string(outstr, HeadlineTag, headline);        /* output headline     */
    }

/*-----------------------------------------------------------------------------
 *  close_binary_output  --  closes binary output stream;
 *                           directed by  write_final_output()  in out.c
 *-----------------------------------------------------------------------------
 */
void   close_binary_output()
    {
    strclose(outstr);
    }

/*-----------------------------------------------------------------------------
 *  put_full_snapshot  --  write binary output data of a complete snapshot
 *-----------------------------------------------------------------------------
 */
void  put_full_snapshot(the_state)
stateptr  the_state;
    {
    int  nbody;
    real  tnow;
    systptr  sys;
    diagptr  diags;

    sys = System(the_state);
    tnow = Tnow(sys);
    nbody = Nbody(sys);
    diags = Diags(the_state);

    put_set(outstr, SnapShotTag);		/* output snapshot set       */
    put_parameters(&nbody, &tnow);		/* output parameters         */
    put_particles(sys);	                        /* output body data          */
    put_diagnostics(diags);                     /* output diagnostics        */
    put_tes(outstr, SnapShotTag);		/* finish snapshot set       */
    fflush(outstr);                             /* write out to disk         */
/*
 * the following announcement is channeled through a procedure in the file
 * out.c  for modularity's sake, leaving final output to  out.c .
 */
    announce("\n\tFull Snapshot output done at time t_now = %lf\n", tnow);
    }

/*-----------------------------------------------------------------------------
 *  put_diagnostics_snapshot  --  write binary output data, in snapshot format,
 *                                containing only paramteters and diagnostics.
 *-----------------------------------------------------------------------------
 */
void  put_diagnostics_snapshot(the_state)
stateptr  the_state;
    {
    int  nbody;
    real  tnow;
    systptr  sys;
    diagptr  diags;

    sys = System(the_state);
    tnow = Tnow(sys);
    nbody = Nbody(sys);
    diags = Diags(the_state);

    put_set(outstr, SnapShotTag);		/* output snapshot set       */
    put_parameters(&nbody, &tnow);		/* output parameters         */
    put_diagnostics(diags);                     /* output diagnostics        */
    put_tes(outstr, SnapShotTag);		/* finish snapshot set       */
    fflush(outstr);                             /* write out to disk         */
/*
 * the following announcement is channeled through a procedure in the file
 * out.c  for modularity's sake, leaving final output to  out.c .
 */
    announce("\n\tDiagnostics Snapshot output done at time t_now = %lf\n",
	     tnow);
    }

/*-----------------------------------------------------------------------------
 *  put_parameters --  output the number of particles and the current time
 *-----------------------------------------------------------------------------
 */
local void  put_parameters(nbodyptr, timeptr)
int  *nbodyptr;
real  *timeptr;
    {
    put_set(outstr, ParametersTag);
    put_data(outstr, NobjTag, IntType, nbodyptr, 0);
    put_data(outstr, TimeTag, RealType, timeptr, 0);
    put_tes(outstr, ParametersTag);
    }

/*-----------------------------------------------------------------------------
 *  put_particles --  output the information about each particle
 *-----------------------------------------------------------------------------
 */
local void  put_particles(sys)
systptr  sys;
    {
    int  nbody;
    bodyptr  bodies;

    nbody = Nbody(sys);
    bodies = Bodies(sys);

    put_set(outstr, ParticlesTag);
    put_data(outstr, CoordSystemTag, IntType, &coordsys, 0);
    put_mass(bodies, nbody);
    put_phase(bodies, nbody);
#ifndef REGULARIZATION
    put_acc(bodies, nbody);
    put_pot(bodies, nbody);
#endif
    put_tes(outstr, ParticlesTag);
    }

/*-----------------------------------------------------------------------------
 *  put_mass  --  output the mass of each particle
 *-----------------------------------------------------------------------------
 */
local void   put_mass(bodies, nbody)
bodyptr  bodies;
int  nbody;
    {
    real *tmpmass, *massptr;
    bodyptr body_i;

    tmpmass = (real *) malloc((unsigned)nbody*sizeof(real));
    if (tmpmass == NULL)
	error("put_mass: not enought memory left for tmpmass[%d]\n", nbody);

    massptr = tmpmass;				/* point to mass array */
    for (body_i = bodies; body_i - bodies < nbody; body_i++)
	*massptr++ = Mass(body_i);		/*   copy individual masses */
    put_data(outstr, MassTag, RealType, tmpmass, nbody, 0);
    free(tmpmass);
    }

/*-----------------------------------------------------------------------------
 *  put_phase  --  output the basic phase space information of each particle:
 *                 positions and velocities
 *-----------------------------------------------------------------------------
 */
local void  put_phase(bodies, nbody)
bodyptr  bodies;
int  nbody;
    {
    real *tmpphase, *pp;
    bodyptr body_i;

    tmpphase = (real *) malloc((unsigned)2*NDIM_IO*nbody*sizeof(real));
    if (tmpphase == NULL)
	error("put_phase: not enought memory left for tmpphase[%d]\n",
                                                              2*NDIM_IO*nbody);

    pp = tmpphase;			/* point to phase array */
    for (body_i = bodies; body_i - bodies < nbody; body_i++)
        {
	SETV_OUT(pp, Pos(body_i));		       /* copy position data */
	pp += NDIM_IO;
	SETV_OUT(pp, Vel(body_i));		       /* copy velocity data */
	pp += NDIM_IO;
        }
    put_data(outstr, PhaseSpaceTag, RealType, tmpphase, nbody, 2, NDIM_IO, 0);
    free(tmpphase);
    }

/*-----------------------------------------------------------------------------
 *  put_pot  --  output the potential energy of each particle
 *-----------------------------------------------------------------------------
 */
local void   put_pot(bodies, nbody)
bodyptr  bodies;
int  nbody;
    {
    real *tmpphi, *pp;
    bodyptr body_i;

    tmpphi = (real *) malloc((unsigned)nbody*sizeof(real));
    if (tmpphi == NULL)
	error("put_phi: not enought memory left for tmpphi[%d]\n", nbody);

    pp = tmpphi;				/* point to phi array */
    for (body_i = bodies; body_i - bodies < nbody; body_i++)
	*pp++ = Pot(body_i);				/*   copy potentials */
    put_data(outstr, PotentialTag, RealType, tmpphi, nbody, 0);
    free(tmpphi);
    }

/*-----------------------------------------------------------------------------
 *  put_acc  --  output the accelerations of each particle
 *-----------------------------------------------------------------------------
 */
local void  put_acc(bodies, nbody)
bodyptr  bodies;
int  nbody;
    {
    real *tmpacc, *pp;
    bodyptr body_i;

    tmpacc = (real *) malloc((unsigned)NDIM_IO*nbody*sizeof(real));
    if (tmpacc == NULL)
	error("put_acc: not enought memory left for tmpacc[%d]\n",
                                                                NDIM_IO*nbody);

    pp = tmpacc;			/* point to acc array */
    for (body_i = bodies; body_i - bodies < nbody; body_i++)
        {
	SETV_OUT(pp, Acc(body_i));		/* copy acceleration data */
	pp += NDIM_IO;
        }
    put_data(outstr, AccelerationTag, RealType, tmpacc, nbody, NDIM_IO, 0);
    free(tmpacc);
    }

/*-----------------------------------------------------------------------------
 *  put_diagnostics  --  output some diagnostics (more to be implemented later)
 *-----------------------------------------------------------------------------
 */
void   put_diagnostics(diags)
diagptr  diags;
    {
    real  cputime_helper;
    real  etot[3];           /* [0]:total, [1]:kinetic, [2]:potential energy */
    
    cputime_helper = cputime();
    
    etot[0] = DIAGetot(diags)[CURRENT];
    etot[1] = DIAGekin(diags)[CURRENT];
    etot[2] = DIAGepot(diags)[CURRENT];

    put_set(outstr, DiagnosticsTag);
    put_data(outstr, EnergyTag, RealType, etot, 3, 0);
    put_data(outstr, "nsteps", IntType, &DIAGnsteps(diags), 0);
    put_data(outstr, "cputime", RealType, &cputime_helper, 0);
    put_tes(outstr, DiagnosticsTag);
    }

/* endof: binaryout.c */
