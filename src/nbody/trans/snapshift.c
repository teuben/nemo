/*
 *  SNAPSHIFT:	shift an nbody system
 *
 *	24-Apr-87	V1.0 created      				PJT
 *       5-May-87       v1.1 more  intelligent read_ and write_snap	PJT
 *	 6-May-87	V1.2 treat headline in appending mode		PJT
 *	 9-Mar-88	V1.3 history(3) used, some dprintf() cleanup	PJT
 *	 6-Jun-88	V2.0 new filestruct				PJT
 *	13-mar-89	V2.1 new shapshot macros			PJT
 *	29-mar-89	V2.2 drange() is now nemoinpd()			PJT
 *      25-may-90       V2.3 added mode parameter                       PJT
 *	nov-93		V3.0 fixed missing data bug
 *	15-mar-95	V3.1 nemo_main, minor formatting		pjt
 *	 4-apr-97	V3.1b SINGLEPREC				pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {
    "in=???\n               Input file (snapshot)",
    "out=???\n              Output file (snapshot)",
    "rshift=0,0,0\n         Shift in position",
    "vshift=0,0,0\n         Shift in velocity",
    "times=all\n            Which times to work on",
    "mode=add\n             Add or subtract the shift? [add|sub]",
    "VERSION=3.1b\n	    4-apr-97 pjt",
    NULL,
};

string usage="shift an N-body system";

#define TIMEFUZZ	0.0001	/* tolerance in time comparisons */

void nemo_main(void)
{
    stream instr, outstr;
    real   tsnap, fact;
    vector rshift, vshift;
    string headline = NULL;
    string times;
    Body *btab = NULL, *bp;
    int i, nbody, bits, nrshift, nvshift;
    char *mode;

    nrshift = nemoinpr(getparam("rshift"),rshift,NDIM);
    if (nrshift!=NDIM)
    	error("rshift= needs %d values, got %d", NDIM, nrshift);

    nvshift = nemoinpr(getparam("vshift"),vshift,NDIM);
    if (nvshift!=NDIM)
    	error("vshift- needs %d values, got %d", NDIM, nvshift);
			
    dprintf(1,"(rshift,vshift): %g %g %g %g %g %g\n",
	rshift[0],rshift[1],rshift[2],vshift[0],vshift[1],vshift[2]);

    mode = getparam("mode");
    switch (*mode) {
      case 's':
      case '-':
                fact = -1.0;    break;
                
      case 'a':
      case '+':
                fact = 1.0;     break;
      default:
                dprintf(0,"Unknown mode %s, add assumed\n",mode);
                fact = 1.0;	break;
    }
    for (i=0; i<NDIM; i++) {
        rshift[i] *= fact;
        vshift[i] *= fact;
    }

    times = getparam("times");
					/* open files */
    instr = stropen(getparam("in"), "r");
    outstr = stropen(getparam("out"), "w");

    get_history(instr);		/* get data history from input file */
    put_history(outstr);	/* write "    "     to output file */
    for(;;) {			/* do the work in an infinite loop */
    	get_history(instr);     /* skip any unneeded history/headline stuff */
        if (!get_tag_ok(instr, SnapShotTag))
		break;			/* done with work in loop */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
#if 0
        if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0) {
            unlink(getparam("out"));
            error("%s: essential data missing\tbits = %x\n", getargv0(), bits);
        }
#else
        if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0) continue;
#endif        	
        if ((bits & TimeBit) == 0)
        	tsnap = 0.0;
        else if (!streq(times,"all") && !within(tsnap, times, TIMEFUZZ))
        	continue;		/* however skip this snapshot */
        dprintf (1,"Snapshot time=%f shifting\n",tsnap);
        for (bp = btab; bp < btab+nbody; bp++) {
            for (i=0; i<NDIM; i++)
               Pos(bp)[i] += rshift[i];
            for (i=0; i<NDIM; i++)
               Vel(bp)[i] += vshift[i];
        }
        put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    }
    strclose(instr);
    strclose(outstr);
}


