/*
 * SNAPADDBH.C:    Add a black hole to the center of a snapshot, 
 *                 and modify the velocities of the stars accordingly
 *
 *      22-mar-04:  after Nelly Mouawad's ideas originally in a modified mkplummer
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <snapshot/snapshot.h>
#include <snapshot/barebody.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {	                 	/* DEFAULT INPUT PARAMETERS */
    "in=???\n	    Input snapshot file",
    "out=???\n      Output snapshot file",
    "bh=1\n         Mass of the BH to add",
    "VERSION=1.0\n  22-mar-04 PJT",
    NULL,
};

string usage = "add a black hole to a snapshot";


nemo_main()
{
    stream instr, outstr;
    Body *btab = NULL;
    int i, nbody, bits;
    real tsnap, bh;

    instr = stropen(getparam("in"), "r");           /* open input file */
    outstr = stropen(getparam("out"), "w");

    bh = getdparam("bh");

    get_history(instr);
    while (get_tag_ok(instr, SnapShotTag)) {
	get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (bits & PhaseSpaceBit) {
            dprintf(1,"Processing time %g bits=0x%x\n",tsnap,bits);
	    add_bh(btab, nbody, bh);
	    put_snap(outstr, &btab, &nbody, &tsnap, &bits);
	}
    }
    strclose(instr);
    strclose(outstr);
}

add_bh(Body *btab, int nbody, real bh)
{
  warning("this doesn't do a thing yet");
}
