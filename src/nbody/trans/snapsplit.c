/*
 * SNAPCUT.C: cut an N-body snapshot into pieces - useful for low
 *		memory machines which can process snapshots in series
 *		The reverse would be snapmerge ?
 *
 *      27-jan-94  V1.0 Created             peter teuben
 *	28-jan-94  V1.1 allow nbody= to have more entries
 *	22-dec-95  released 'as is'
 *      14-sep-02  V2.0 use mstropen() to optionally split into different files
 *      13-jun-07  V2.1 using within() from stdinc.h   WD
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
    "in=???\n           Input snapshot",
    "out=???\n          Output snapshot (multiple allowed via %)",
    "nbody=\n           Size of one (or more) snapshot(s)",
    "nsnap=\n           Number of pieces to cut a snapshot into",
    "times=all\n        Series of times-ranges to select",
    "VERSION=2.1\n	13-jun-07 PJT/WD",
    NULL,
};

string usage="cut an N-body snapshot in pieces for serial processing";

#define TIMEFUZZ	0.0001	/* tolerance in time comparisons */

#define MAXSNAP 1024

/* 
 * commented out as it is already in stdinc.h, causing compilation error WD
 *
   extern bool within(real, string, real);
*/

nemo_main()
{
    stream instr, outstr;
    mstr *mp;
    real   tsnap;
    string times, precision, keep;
    Body   *btab = NULL, *bp;
    int    i, n1, nc, nsnap, nin, nout, bits;
    int    nbody, nbodys[MAXSNAP];

    times = getparam("times");

    if (hasvalue("nbody")) {
        nsnap = nemoinpi(getparam("nbody"),nbodys,MAXSNAP);
        if (nsnap<1) error("parsing nbody=");
        nbody = nsnap==1 ? nbodys[0] : 0;
        
    } else if (hasvalue("nsnap")) {
        nsnap = getiparam("nsnap");
        /* */
    } else
        error("Need either nbody= or nsnap= to cut into pieces");

    instr = stropen(getparam("in"), "r");
    mp = mstr_init(getparam("out"));

    for (;;) {
    	get_history(instr);		/* skip over stuff we can forget */
        if (!get_tag_ok(instr, SnapShotTag))
		break;			/* done with work in loop */

        get_snap(instr, &btab, &nin, &tsnap, &bits);
        if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0) {
	    continue;       /* just skip */
        }

        if (!streq(times,"all") && !within(tsnap, times, TIMEFUZZ))
        	continue;		/* now skip this snapshot */
        dprintf (0,"Snapshot time=%f to be copied\n",tsnap);

        if (nbody > 0)
            nout = nbody;
        else
            nout = (nin+nsnap-1)/nsnap;

        for (bp=btab, i=0, nc=0; i<nin; bp+=nout, i+=nout, nc++) {
            if (nin-i > nout)            
                n1 = nout;
            else
                n1 = nin - i;

	    outstr = mstr_open(mp, "w");
	    dprintf(0,"mstropen: count=%d\n",mstr_count(mp));
	    if (mstr_count(mp)==1 || mstr_multi(mp))
	      put_history(outstr);
            put_snap(outstr, &bp, &n1, &tsnap, &bits);
        }
        if (n1 != nout)
            dprintf(0,"Wrote %d %d-body snapshots + 1 %d-body snapshot\n",
                nc-1,nout,n1);
        else
            dprintf(0,"Wrote %d %d-body snapshots\n",nc,nout);
    }
    strclose(instr);
    mstr_close(mp);
}
