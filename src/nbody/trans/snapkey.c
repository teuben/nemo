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
    "VERSION=0.2\n            17-feb-2020 PJT",
    NULL,
};

string usage="(re)assign keys to a snapshot";

string cvsid="$Id$";

#define TIMEFUZZ	0.0001	/* tolerance in time comparisons */
#define MAXNBODY        100000

extern real_proc getrfunc(string , string , string , int *);

void nemo_main(void)
{
    stream instr, inmassstr, outstr, ccdstr;
    real   tsnap, tsnapmass, mass_star, mtot, mrange[2], norm;
    real   xpos, ypos, xmin, ymin, idx, idy;
    Body  *btab = NULL, *bp;
    Body  *bmasstab = NULL, *bmassp;
    imageptr iptr = NULL;
    int i, n, nbody, nbodymass, bits, bitsmass, seed, ix, iy;
    rproc_body  bfunc;
    real_proc   mfunc;
    int keys[MAXNBODY], nkeys;
    bool  Qnorm, first = TRUE;

    instr = stropen(getparam("in"), "r");
    outstr = stropen(getparam("out"), "w");
    nkeys = nemoinpi(getparam("keys"),keys,MAXNBODY);

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
        } else             
            error("bad flow logic");

	bits |= KeyBit;    /* turn mass bit on anyways */
        put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    }

    strclose(instr);
    strclose(outstr);
}
