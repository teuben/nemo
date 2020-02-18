/*
 *  SNAPKEY:  set keys in a snapshot - see also snapmask
 *
 *      17-feb-2020     0.1    quick hack from snapmass
 *
 */
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <image.h>
#include <loadobj.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

#include <bodytrans.h>


string defv[] = {
    "in=???\n		      input (snapshot) file",
    "out=???\n                output (snapshot) file",
    "keys=\n                  Manually supply all keys",
    "keyfile=\n		      ascii table with keys",
    "key=\n                   bodytrans expression",
    "VERSION=0.3\n            17-feb-2020 PJT",
    NULL,
};

string usage="(re)assign keys to a snapshot";

string cvsid="$Id$";

#define TIMEFUZZ	0.0001	/* tolerance in time comparisons */
#define MAXNBODY        100000

extern btiproc btitrans(string);


void nemo_main(void)
{
    stream instr, outstr, tabstr;
    real   tsnap;
    Body  *btab = NULL, *bp;
    int i, n, nbody, bits;
    int keys[MAXNBODY], nkeys;
    bool  first = TRUE;
    btiproc   kfunc = NULL;
 
    instr = stropen(getparam("in"), "r");
    outstr = stropen(getparam("out"), "w");
    nkeys = nemoinpi(getparam("keys"),keys,MAXNBODY);
    if (hasvalue("key"))
      kfunc = btitrans(getparam("key"));

    for(;;) {
					/* input main data */
    	get_history(instr);                    /* skip history & comments */
        if (!get_tag_ok(instr, SnapShotTag)) {
            break;
            error("Snapmass (in): Need a snapshot");
        }
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (bits&KeyBit)
            dprintf(0,"Warning: existing keys overwritten\n");
        if (!(bits&TimeBit)) {
            dprintf(1,"Warning: time=0.0 assumed\n");
            tsnap = 0.0;
        }
        
        if (first) {
            put_history(outstr);
            first = FALSE;
        }
	if (nkeys > 0) {
	  for (bp=btab, i=0; bp<btab+nbody; bp++, i++)
	    Key(bp) = keys[i];
        } else if (kfunc != NULL) {
	  for (bp=btab, i=0; bp<btab+nbody; bp++, i++)
	    Key(bp) = (*kfunc)(bp,tsnap,i);
	} else
            error("no keys assigned");

	bits |= KeyBit;    /* turn mass bit on anyways */
        put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    }

    strclose(instr);
    strclose(outstr);
}
