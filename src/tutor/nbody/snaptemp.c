/*
 *  SNAPTEMP:  template to operate on snapshot files
 *
 *  With the $NEMOBIN/cc script it should compile as follows:
 *             cc -o snaptemp snaptemp.c -lnemo -lm
 *
 *      26-jan-94	V1.0    Created                 Peter Teuben
 *	 4-nov-94       renamed from snaptrans to snaptemp	pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {   /* Nemo_Keys: "key=val\n    help",       <--- required */
    "in=???\n			input (snapshot) file",
    "out=???\n                  output (snapshot) file",
    "times=all\n                times to process",
    "VERSION=1.0a\n		4-nov-94 PJT",
    NULL,
};

string usage = "template for snapshot operations";	/*   <--- required */

void nemo_main(void) /* this replaces main(argc,argv)        <--- required */
{
    stream instr, outstr;
    string times = getparam("times");
    real   tsnap;
    Body  *btab = NULL, *bp;
    int    i, nbody, bits;
    bool   first = TRUE;

    instr = stropen(getparam("in"), "r");
    outstr = stropen(getparam("out"), "w");

    for (;;) {    /* infinite loop, broken only when ran out of snapshots */
    	get_history(instr);                    /* read history */
        if (!get_tag_ok(instr, SnapShotTag)) break; /* check if done */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);  /* get next */
        
        /* Operate on the snapshot here : */
        for (bp=btab, i=0; i<nbody; bp++, i++) {
            /* all the work goes here */
        }
        if (first) {
            put_history(outstr);
            first = FALSE;
        }
        put_snap(outstr, &btab, &nbody, &tsnap, &bits);     /* output */
    }
    strclose(instr);
    if (first) {
    	warning("No snapshots processed");
  	strdelete(outstr,TRUE);
    } else
        strclose(outstr);
}
