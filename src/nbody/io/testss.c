/*
 * TESTSS.C: test routine for get_snap  and put_snap methods.
 *	     (formarly called SNAPTEST)
 *
 *  xx-xxx-8x   Created in the dark ages                            JEB
 *  12-sep-90   Allow a very primitive edit (nlist,i,option,data)   PJT
 *  20-feb-92   fixed a few portability problems		    PJT
 *   7-mar-92   happy GCC2.0 - renamed it from SNAPTEST to TESTSS   pjt
 */

#include <stdinc.h>
#include <vectmath.h>
#include <getparam.h>
#include <filestruct.h>

string defv[] = {
    "in=???\n       Input file",
    "out=???\n      Output file",
    "nlist=2\n	    Maximum number of bodies to print",
    "i=-1\n         Star number to process (>=0)",
    "option=\n      Option [story]",
    "data=\n        Data supplied with option",
    "times=\n       Times to process; defaults to all",
    "VERSION=1.1a\n 7-mar-92  PJT",
    NULL,
};

string usage = "test routine for <get_snap>  and <put_snap> methods";

/*
 * Include files from local snapshot directory.
 */

#include "snapshot.h"
#include "mybody.h"
#include "myget_snap.c"
#include "myput_snap.c"

void
nemo_main()
{
    stream instr, outstr;
    Body *btab = NULL, *bp;
    char *data, *option, *cp, *allocate();
    int nbody, bits, i, nlist;
    real tsnap;

    i = getiparam("i");
    nlist = getiparam("nlist");
    option = getparam("option");
    data = getparam("data");

    instr = stropen(getparam("in"), "r");               /* INPUT */
    get_history(instr);
    if (hasvalue("times"))
	get_snap_by_t(instr, &btab, &nbody, &tsnap, &bits, getparam("times"));
    else
	get_snap(instr, &btab, &nbody, &tsnap, &bits);
    strclose(instr);

    printf("get_snap returned with bits= 0x%x\n", bits);       /* PROCESS */
    printf("  nbody = %d\n", nbody);
    if (bits & TimeBit)
	printf("  tsnap = %f\n", tsnap);
    for (bp = btab; bp < btab + MIN(nbody,nlist); bp++) {
	if (bits & MassBit)
	    printf("  body%3d: mass = %f\n", bp - btab, Mass(bp));
	if (bits & PhaseSpaceBit) {
	    printf("           pos  = [%f,%f,%f]\n",
		   Phase(bp)[0][0], Phase(bp)[0][1], Phase(bp)[0][2]);
	    printf("           vel  = [%f,%f,%f]\n",
		   Phase(bp)[1][0], Phase(bp)[1][1], Phase(bp)[1][2]);
	}
    }
    if (i>=0 && *option) {
        if (streq(option,"story")) {
            if (*data) {
                if (Story(btab+i)) {
                    cp = allocate(strlen(data)+strlen(Story(btab+i))+1);
                    strcpy(cp,Story(btab+i));
                    strcat(cp,"\n");
                    strcat(cp,data);
                    data = cp;
                }
                Story(btab+i) = data;
            } else
                error("No story supplied! (data=)\n");
        } else
            error("Not a valid option option={story,}\n");
    }

    outstr = stropen(getparam("out"), "w!");         /* OUTPUT overwrite */
    put_history(outstr);
    put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    strclose(outstr);
} /* nemo_main */
