/*
 * NBODYIO.C   C-version I/O for Fortran-C interface of Aarseth's nbody0.f
 *             Has a simple F2C interface  
 *		22-nov-91  fixed aix interface
 *		 6-jan-00  converted to F77_FUNC macros
 *		21-jan-00  added reset option 
 *		 9-apr-01  fixed some -DSINGLEPREC errors
 *              21-feb-04  prototypes to keep -Wall silent
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>

#include <filestruct.h>
#include <history.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>
      
local double l_eta;                  /* local copies of things for Aarseth */
local double l_deltat;               /* now double, used to be float */
local double l_tcrit;
local double l_eps2;

local int   nbody;                  /* actual number of bodies */
local int   bits;                   /* bitmask what's in snapshot */
local bool  Qstep;
local bool  Qreset;
local bool  Qf3dot;

local stream instr;                    /* - pointer to input file */
local stream outstr;                   /* - pointer to output file */

local Body *btab = NULL;               /* - pointer to body structure */

extern double cputime(void);

#include "proto.h"			/* fortran to c prototypes */

/*============================================================================*/

/*
 *  pars_2_aarseth: kludge to prepare transport parameters to fortran
 *		and start reading the input snapshot file
 *
 */

void pars_2_aarseth(void)
{
    real time;
    
    l_eta = (double) getdparam("eta");
    l_deltat = (double) getdparam("deltat");
    l_tcrit = (double) getdparam("tcrit");
    l_eps2 = (double) sqr(getdparam("eps"));

    Qreset = getbparam("reset");
    Qf3dot = getbparam("f3dot");

    if (hasvalue("options"))
    	Qstep = TRUE;
    else
    	Qstep = FALSE;

    
    if (NDIM != 3)
        error("%d: Code not supported for non 3D",NDIM);
    instr = stropen(getparam("in"),"r");     /* open input file */
    get_history(instr);             /* get history from input file */
    get_snap (instr, &btab, &nbody, &time, &bits);  /* read snapshot in */

    if ((bits & TimeBit)==0)        /* set time to zero if not present */
        time=0.0;
    if ((bits & MassBit)==0)            /* check if essentials are present */
        error ("No masses in input snapshot");
    if ((bits & PhaseSpaceBit)==0)
        error ("No phase space coordinates in input snapshot");

    if (time >= (double) l_tcrit)
       error("Input snapshot has larger timestamp (%g) than required tcrit(%g)",
			time, l_tcrit);

    outstr = stropen(getparam("out"),"w");      /* open output file */
    put_history(outstr);            /* copy history to output file */
    bits = (TimeBit | MassBit | PhaseSpaceBit);     /* set things to copy */
    if (Qstep) bits |= AuxBit;
}

/* 
 *  call_2_aarseth:	the actual routine which does the fortran calling
 */
void call_2_aarseth(void)
{
    nbody0();           /* call fortran hard core worker routine */
}

/* 
 * (fortran) INPARS: read global integration parameters
 */

void inpars (nmax, n, eta, deltat, tcrit, eps2, reset, use3dot)
     int *nmax, *n, *reset, *use3dot;
     double *eta, *deltat, *tcrit, *eps2;
{
    if (nbody > *nmax)
        error("Too many particles (%d) in inputfile, recompile (%d) program",
			nbody, *nmax);
    *n = nbody;
    *eta    = l_eta;
    *deltat = l_deltat;
    *tcrit  = l_tcrit;
    *eps2   = l_eps2;
    *reset  = (Qreset ? 1 : 0);
    *use3dot= (Qf3dot ? 1 : 0);
}

/*
 *  (fortran) INBODS: read particle masses and phases into array
 */

void inbods (n, m, x, v)
int   *n;
double *m, *x, *v;
{
    Body *bp;
    int i;
			/* WARNING: what about double and DOUBLE here ? */
    for (bp=btab, i=0; i<*n; bp++, i++) {
        *m++ = Mass(bp);            /* copy mass & increase pointer */
        SETV(x,Pos(bp));            /* position */
        SETV(v,Vel(bp));            /* and velocity */
        x += NDIM;  v += NDIM;      /* increase array pointers */
    }
}

/*
 *  (fortran) OUTBODS: outputs characteristics of particle i  (i=1,...nbody)
 *        
 *      In C version just prepares data for flushing later on
 */

void outbods (m, x, v, step, i)
double *m, *x, *v, *step;
int *i;
{
    Body *bp;

    if (*i < 0 || *i > nbody)
        error("OUTBODS: counting error particle\n");

    bp = (btab + (*i - 1));   /* point to proper position in Body */

    Mass(bp) = *m;
    SETV(Pos(bp), x);
    SETV(Vel(bp), v);
    if (Qstep) Aux(bp) = *step;
}

/*
 * (fortran) OUTENE: output energy
 *
 *   In C: also flushes particles buffer to file
 */

void outene (tnext, nsteps, e)
double *tnext, *e;
int   *nsteps;
{
    real time;

    dprintf(0,"time = %g   steps = %d   energy = %g cpu = %10.3g min\n",
            *tnext, *nsteps, *e, cputime());

    time = (real) *tnext;
    put_snap (outstr, &btab, &nbody, &time, &bits);
}

