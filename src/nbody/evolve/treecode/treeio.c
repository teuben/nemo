/*
 * TREEIO.C   C-version I/O for Fortran-C interface of Hernquist' treecode
 *          Has a crude F2C interface 
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>

#include <filestruct.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

local int   nbody;                  /* actual number of bodies */
local int   bits;                   /* bitmask what's in snapshot */

local stream instr;                    /* - pointer to input file */
local stream outstr;                   /* - pointer to output file */

local Body *btab = NULL;               /* - pointer to body structure */
local char outbuff[256];                /* local buffer to print messages */

/*============================================================================*/
/*                  Fortran-to-C interface stuff:
 *  Since we only handle simple parameters here (i.e. no character
 *  variables), the only thing to be done is to properly define the name
 *  of the C routine, because it will be called as a fortran subroutine. 
 *  General rule: avoid character variables and life is fairly simple in F2C.
 *
 *  In VMS the C names stay the same
 *  in UNICOS the C name is to be in UPPER CASE.
 *  In most (BSD) unix versions C name gets an underscore (_) appended
 */

#if !defined(VMS)
#if defined(unicos)
#define treecode TREECODE
/*	will not work, due to character I/O */
#else
#define treecode treecode_
#define second   second_
#define inparams inparams_
#define inbods   inbods_
#define outbods  outbods_
#define outenrgy outenrgy_
#define startout startout_
#define outcpu   outcpu_
#define stopout  stopout_
#define outterm  outterm_
#define outerror outerror_
#define outlog   outlog_
#define outbox   outbox_
#define cellcom  cellcom_
#define paramcom paramcom_
#define bodycell bodycell_
#define timecom	 timecom_
#define cpucom   cpucom_
#endif
#endif

#include "nmax.h"
#include "common.h"

extern struct fcb1 paramcom;		/* access common blocks from fortran */
extern struct fcb3 cellcom;
extern struct fcb5 bodycell;
extern struct fcb8 timecom;
extern struct fcb9 cpucom;
/*============================================================================*/

/*
 *  pars_2_lars: kludge to prepare transport parameters to fortran
 *		and start reading the input snapshot file
 *
 */

pars_2_lars()
{
    double time;
    
    timecom.dtime = (float) getdparam("dtime");
    timecom.nsteps = (float) getdparam("nsteps");
    timecom.nout = (float) getdparam("nout");
    paramcom.tol = (float) getdparam("tol");
    paramcom.eps = (float) getdparam("eps");
    paramcom.usequad = (int) getbparam("usequad");

    instr = stropen(getparam("in"),"r");     /* open input file */
    get_history(instr);             /* get history from input file */
    get_snap (instr, &btab, &nbody, &time, &bits);  /* read snapshot in */
    paramcom.nbodies = nbody;           /* copy to fortran common blocks */
    if ((bits & TimeBit)==0)        /* set time to zero if not present */
        time=0.0;
    if ((bits & MassBit)==0)            /* check if essentials present */
        error ("No masses in input snapshot\n");
    if ((bits & PhaseSpaceBit)==0)
        error ("No phases in input snapshot\n");

    outstr = stropen(getparam("out"),"w");      /* open output file */
    put_history(outstr);            /* copy history to output file */
    bits = (TimeBit | MassBit | PhaseSpaceBit);     /* set things to copy */
}

/* 
 *  call_2_lars:	the actual routine which does the fortran calling
 */
call_2_lars()
{
    treecode();           /* call fortran hard core worker routine */
}

/* 
 * (fortran) INPARAMS: dummy routine for lars
 */

inparams ()
{
}

/*
 *  (fortran) INBODS: read particle masses and phases into common array's
 */
inbods ()
{
    Body *bp;
    int i;

    dprintf(0,"INBODS: moving %d particles\n",nbody);
    for (bp=btab, i=0; i<nbody; bp++, i++) {
	bodycell.mass[i] = (float) Mass(bp);
	bodycell.pos[i]            = (float) Pos(bp)[0];
	bodycell.pos[i+NBODCELL]   = (float) Pos(bp)[1];
	bodycell.pos[i+2*NBODCELL] = (float) Pos(bp)[2];
	bodycell.vel[i]         = (float) Vel(bp)[0];
	bodycell.vel[i+NMAX]    = (float) Vel(bp)[1];
	bodycell.vel[i+2*NMAX]  = (float) Vel(bp)[2];
    }
}

startout()		/* dummy for lars */
{
}


/*
 *  (fortran) OUTBODS: outputs characteristics of particle 
 *        
 */

outbods ()
{
    Body *bp;
    int i;

    dprintf(0,"outputting particles:\n");
    for (i=0; i<nbody; i++)
      dprintf(0,"1: %g %g %g %g %g %g\n",
	bodycell.pos[i], bodycell.pos[i+NBODCELL],bodycell.pos[i+2*NBODCELL],
	bodycell.vel[i], bodycell.vel[i+NMAX],    bodycell.vel[i+2*NMAX]);


}

outbox()
{
    dprintf(0,"expanding the box %g %g %g %g\n",
	cellcom.rmin[0], cellcom.rmin[1], cellcom.rmin[2], cellcom.rsize);
}
 
outenrgy ( mtot, e, ek, ep, am, cmpos, cmvel)
float *mtot, *e, *ek, *ep, *am, *cmpos, *cmvel;
{
    dprintf(0,"outenrgy: \n");

    printf ("mtot = %12.4g\n", *mtot);
    printf ("e, ek, ep = %12.4g %12.4g %12.4g\n",*e, *ek, *ep);
    printf ("amx, amy, amz = %12.4g  %12.4g  %12.4g \n",am[0], am[1], am[2]);
    printf ("cmpos = %12.4g %12.4g %12.4g \n",cmpos[0], cmpos[1], cmpos[2]);
    printf ("cmvel = %12.4g %12.4g %12.4g \n",cmvel[0], cmvel[1], cmvel[2]);
}

outlog (istep)
int *istep;
{
    dprintf(0,"outlog: %d\n",*istep);
}

outhead(outunit)
int *outunit;
{
    dprintf(0,"outhead: %d\n",*outunit);
}

outterm(message, n, mlen)           /* SHIT F2C interface */
char *message;
int *n, mlen;
{
    dprintf(0,"outterm: \n");

    strncpy(outbuff,message,mlen);
    outbuff[mlen] = '\0';
    dprintf(0,"%s %d\n",outbuff,*n);
}

outerror(message,mlen)          /* f2c */
char *message;
int mlen;
{
    dprintf(0,"outerror: \n");

    strncpy(outbuff,message,mlen);
    outbuff[mlen] = '\0';
    dprintf(0,"ERROR: %s\n",outbuff);
}

outcpu()
{
    dprintf(0,"outcpu: \n");

    printf("\n\nTotal CPU time uses (seconds): %f\n",
                        cpucom.cputime1 - cpucom.cputime0);
}

stopout()           /* dummy for lars */
{
}


/* 
 * FORTRAN second:  return CPU so far in seconds 
 */
second ( sec )
float *sec;
{
    double cputime();

    *sec = (float) (60.0*cputime());
    return ;
}

#ifdef sun3
		/* fool the bad SUN OS 4.0.1 compiler  */
units()
{
	error("ERROR: units called\n");
}
#endif
