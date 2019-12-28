/*
 * SNAPCOPY.C: copy an N-body snapshot.
 *	
 *	11-Apr-89	V1.0	created 		PJT
 *	26-oct-90          a    ??
 *	11-jan-97	   b    report # copied, NEMOV2 style
 *	25-mar-97	   c    minor fix SINGLEPREC
 *       1-apr-01          d    compiler warning
 *      13-feb-04       V1.1f    silenced more compiler warnings (shetty bug?)
 *      15-nov-06        1.2    set time to 0 if it was absent     PJT/AP
 *    27-dec-2019        1.3    special case body= selection
 *
 *	BUG: should optionally copy other sets within the snapshot
 *	     set, e.g. diagnostics and story
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
				/* new filestruct */
#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

#include <bodytransc.h>

string defv[] = {
    "in=???\n           Input snapshot",
    "out=???\n          Output snapshot",
    "select=1\n         Selection expression",
    "times=all\n        Times to select",
    "precision=double\n Precision of results to store (double/single) [unused]",
    "keep=all\n         Items to copy in snapshot",
    "ibody=-1\n         One body to select (not implemented yet)",
    "VERSION=1.3\n      27-dec-2019 PJT",
    NULL,
};

string usage="copy an N-body snapshot";

string cvsid="$Id$";


#define TIMEFUZZ	0.0001	/* tolerance in time comparisons */

void nemo_main(void)
{
    stream instr, outstr;
    real   tsnap;
    string times, precision, keep;
    Body   *btab = NULL, *bpi, *bpo;
    int    i, nbody, nout, nreject, bitsi, bitso, vis, visnow, vismax;
    bool   Qall;
    int    ibody = getiparam("ibody");
    iproc_body sfunc;

    times = getparam("times");
    sfunc = btitrans(getparam("select"));
    instr = stropen(getparam("in"), "r");
    outstr = stropen(getparam("out"), "w");
    precision = getparam("precision");	/* unused */
    if (!streq(precision,"double")) {
        warning("Precision \"%s\" not supported yet, use csf convert=d2f/f2d",
            precision);
    }
    keep = getparam("keep");
    bitso = 0;
    if (scanopt(keep,"all")) {
        Qall = TRUE;
    } else {
        Qall = FALSE;
        if (scanopt(keep,"mass"))
            bitso |= MassBit;
        if (scanopt(keep,"phase"))
            bitso |= PhaseSpaceBit;
        if (scanopt(keep,"phi"))
            bitso |= PotentialBit;
        if (scanopt(keep,"acc"))
            bitso |= AccelerationBit;
        if (scanopt(keep,"aux"))
            bitso |= AuxBit;
        if (scanopt(keep,"key"))
            bitso |= KeyBit;
        if (scanopt(keep,"time"))
            bitso |= TimeBit;
    }


    get_history(instr);
    put_history(outstr);		
    for (;;) {
    	get_history(instr);		/* skip over stuff we can forget */
        if (!get_tag_ok(instr, SnapShotTag))
		break;			/* done with work in loop */
        get_snap(instr, &btab, &nbody, &tsnap, &bitsi);
        if ((bitsi & MassBit) == 0 && (bitsi & PhaseSpaceBit) == 0) {
	    continue;       /* just skip */
        }
        				/* do some initial output */
        if ((bitsi & TimeBit) == 0) {	/* if not time flag */
        	tsnap = 0.0;		/* set time to zero */
		bitsi |= TimeBit;
		warning("time reset to 0");
        } else if (!streq(times,"all") && !within(tsnap, times, TIMEFUZZ))
        	continue;		/* now skip this snapshot */
        if (bitsi & KeyBit)
            dprintf(0,"Warning: Keyfield reinitialized\n");
        for (bpi = btab; bpi < btab+nbody; bpi++)
            Key(bpi) = 0;                    /* set to false */
        visnow = vismax = nreject = 0;
        do {				/* loop through all particles */
            visnow++;
            for (bpi = btab, i=0; i<nbody; bpi++,i++) {
                vis = (*sfunc)(bpi, tsnap, i);
		dprintf(2,"sfunc [%d] = %d\n",i,vis);
                vismax = MAX(vismax,vis);
                if (vis==visnow)
                    Key(bpi) = 1;
            }
        } while (visnow < vismax);          /* do all layers */
        nreject = 0;
        for (bpi = btab, bpo = btab, i=0; i<nbody; bpi++,i++) {
            if (!Key(bpi)) { 
               nreject++;       /* count number of stars rejected */
               continue;
            }
            if (nreject)        /* only copy while out of sync */
                bcopy (bpi, bpo, sizeof(Body));
            Key(bpo) = i;       /* counting from zero */
            bpo++;
            
        }
        nout = nbody - nreject;
        if (nout) {
            if (Qall) {             /* if all old things selected */
                if (nreject)        /* and some stars were rejected */
                    bitsi |= KeyBit;    /* then explicitly add Key field */
            } else
                bitsi = bitso;      
            put_snap(outstr, &btab, &nout, &tsnap, &bitsi);
            dprintf (1,"Snapshot time=%f copied %d particles\n",
				tsnap,nout);
        } else
           dprintf(0,"No particles to copy at tsnap=%f\n",tsnap);
    }
}
