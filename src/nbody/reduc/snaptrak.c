/*
 * SNAPTRAK.C: track centroids of user-specified particle groups.
 *
 *     12-feb-89    V1    JEB
 *     19-nov-93    V1.1  NEMO V2.x    PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <snapshot/body.h>
#include <snapshot/snapshot.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			input file name",
    "out=???\n			output file name",
    "group=???\n		specify group of each particle",
    "times=all\n		range of times to process",
    "headline=\n		random mumble for humans",
    "VERSION=1.1b\n		16-feb-97 PJT",
    NULL,
};

string usage = "track centroids of user-specified particle groups";

Body *bodytab = NULL, *grouptab = NULL;

int nbody, ngroup;

real tsnap;

extern iproc btitrans(string);

nemo_main()
{
    stream instr, outstr;
    string times;
    iproc group;
    int bits, groupbits = TimeBit | PhaseSpaceBit | KeyBit;

    instr = stropen(getparam("in"), "r");
    get_history(instr);
    if (! streq(getparam("headline"), ""))
	set_headline(getparam("headline"));
    group = btitrans(getparam("group"));
    times = getparam("times");
    outstr = stropen(getparam("out"), "w");
    put_history(outstr);
    do {
	get_snap_by_t(instr, &bodytab, &nbody, &tsnap, &bits, times);
	groupbits = groupbits | (bits & MassBit);
	if (bits & PhaseSpaceBit) {
	    snaptrak(group);
	    put_snap(outstr, &grouptab, &ngroup, &tsnap, &groupbits);
	}
    } while (bits != 0);
    strclose(outstr);
}

snaptrak(group)
iproc group;
{
    int i, ig;
    Body *b, *g;

    if (grouptab == NULL) {
	ngroup = 0;
	for (i = 0, b = bodytab; i < nbody; i++, b++) {
	    ig = (group)(b, tsnap, i);
	    ngroup = MAX(ngroup, ig);
	}
	fprintf(stderr, "[snaptrak: allocating %d groups]\n", ngroup);
	grouptab = (Body *) allocate(ngroup * sizeof(Body));
    }
    for (i = 0, g = grouptab; i < ngroup; i++, g++) {
	Mass(g) = 0.0;
	CLRV(Pos(g));
	CLRV(Vel(g));
	Key(g) = 0;
    }
    for (i = 0, b = bodytab; i < nbody; i++, b++) {
	ig = (group)(b, tsnap, i);
	if (ig > ngroup)
	    error("snaptrak: cant expand group array\n");
	if (ig > 0) {
	    g = &grouptab[ig - 1];
	    Mass(g) = Mass(g) + Mass(b);
	    ADDV(Pos(g), Pos(g), Pos(b));
	    ADDV(Vel(g), Vel(g), Vel(b));
	    Key(g) = Key(g) + 1;
	}
    }
    for (i = 0, g = grouptab; i < ngroup; i++, g++) {
	if (Key(g) == 0)
	    error("snaptrak: group %d empty\n", i + 1);
	DIVVS(Pos(g), Pos(g), Key(g));
	DIVVS(Vel(g), Vel(g), Key(g));
    }
}
