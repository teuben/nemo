/*
 * SNAPADDBH.C:    Add a black hole to the center of a snapshot, 
 *                 and modify the velocities of the stars accordingly
 *
 * mode 0:  the sqrt(M/r) is in the same plane as (r,v) but pure centrifugal
 * mode 1:  sqrt(M/r) is added to the original velocity the star had
 *
 *
 * 
 *
 *      22-mar-04:  after Nelly Mouawad's ideas originally in a modified mkplummer
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

string defv[] = {	                 	/* DEFAULT INPUT PARAMETERS */
    "in=???\n	    Input snapshot file",
    "out=???\n      Output snapshot file",
    "bh=1\n         Mass of the BH to add",
    "mode=0\n       Mode how to add",
    "fstream=0\n    Amount of streaming added (not implemented)",
    "VERSION=1.0\n  29-mar-04 PJT",
    NULL,
};

string usage = "add a black hole to a snapshot";


local void add_bh(Body *btab, int nbody, real bh, int mode);

void nemo_main(void)
{
    stream instr, outstr;
    Body *btab = NULL;
    int nbody, bits,mode;
    real tsnap, bh;

    instr = stropen(getparam("in"), "r");           /* open input file */
    outstr = stropen(getparam("out"), "w");

    bh = getdparam("bh");
    mode = getiparam("mode");

    get_history(instr);
    while (get_tag_ok(instr, SnapShotTag)) {
	get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (bits & PhaseSpaceBit) {
            dprintf(1,"Processing time %g bits=0x%x\n",tsnap,bits);
	    add_bh(btab, nbody, bh, mode);
	    put_snap(outstr, &btab, &nbody, &tsnap, &bits);
	}
    }
    strclose(instr);
    strclose(outstr);
}

local void add_bh(Body *btab, int nbody, real bh, int mode)
{
  vector r,v,j,w;
  real vlen, rlen, factor;
  Body *bp;

  for (bp=btab; bp<btab+nbody; bp++) {
    SETV(r,Pos(bp));
    SETV(v,Vel(bp));
    CROSSVP(j,r,v);
    CROSSVP(w,j,r);
    vlen = absv(v);
    rlen = absv(r);
    factor = sqrt(bh/rlen) / (rlen*rlen*vlen);
    SMULVS(w,factor);
    if (mode==0) {
      SETV(Vel(bp),w);
    } else if (mode == 1) {
      SADDV(Vel(bp),w);
    } else {
      error("Unknown mode %d",mode);
    }
  }
}
