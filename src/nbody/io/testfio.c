/*
 * TESTFIO: NEMO fortran I/O (example)
 *	Not meant for export into $NEMOBIN
 *	
 *      16-Nov-88       V1.0 written	                    PJT
 *	15-sep-90	V1.1 CPU tester for Fio vs. Cio     PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

/*#include "fbodyio.c"	*/

string defv[] = {
    "in=???\n			input file name",
    "out=???\n			output file name",
    "VERSION=1.0",		15-sep-90  PJT",
    NULL,
};

#define TIMEFUZZ	0.0001	/* tolerance in time comparisons */

main(argc, argv)
int argc;
string argv[];
{
    char *malloc();
    float *pos, *vel, *mass;
    stream instr, outstr;
    real   mscale, pscale, tsnap;
    vector rscale, vscale;
    string headline = NULL;
    string times;
    Body *btab = NULL, *bp;
    int i, nbody, bits, nrscale, nvscale, ndim;
					/* init NEMO environment */
    initparam(argv, defv);
    ndim = NDIM;
					/* open files */
    instr = stropen(getparam("in"), "r");
    outstr = stropen(getparam("out"), "w");

    /* do only one snapshot */
    	get_history(instr);
        if (!get_tag_ok(instr, SnapShotTag))
		error("No snapshot");
        GetSnap(instr, &btab, &nbody, &tsnap, &bits);
	mass = (float *) malloc(nbody*sizeof(float));
	pos = (float *) malloc(nbody*NDIM*sizeof(float));
	vel = (float *) malloc(nbody*NDIM*sizeof(float));
	if (vel==NULL)
	    error("Cannot allocate data for fortran I/O");

        fio(WRITE_FLAG,nbody,btab,mass,pos,vel);
        fwork_ (&nbody,&ndim,mass,pos,vel);           /* fortran worker */
        fio(READ_FLAG,nbody,btab,mass,pos,vel);

        put_history(outstr);		
        PutSnap(outstr, &btab, &nbody, &tsnap, &bits);
}
