/*
 *  SNAPRC: convert a snapshot to Carlberg's binary 'RV' format
 *
 *	4-May-94	V1.0 created      		TRQ
 */

#include <stdinc.h>
#include <getparam.h>
#include <math.h>
#include <stdlib.h>
#include <vectmath.h>		/* otherwise NDIM undefined */
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>


string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			Input file (snapshot)",
    "out=???\n			Output file (RV format)",
    "times=all\n		Times to select snapshot",
    "VERSION=1.0\n		4-may-94",
    NULL,
};

string usage = "convert a snapshot to Carlberg's binary 'RV' format";

void nemo_main(void)
{
    stream instr, outstr;
    real   tsnap;
    string times;
    Body *btab = NULL, *bp;
    float mass, pos[NDIM], vel[NDIM], t;
    int j, nbody, bits;

            
    instr = stropen(getparam("in"), "r");	/* open input file */
    outstr = stropen(getparam("out"), "w");	/* open output file */
    times = getparam("times");

    get_history(instr);                 /* read history */
    for(;;) {                /* repeating until first or all times are read */
	get_history(instr);
        if (!get_tag_ok(instr, SnapShotTag))
            break;                                  /* done with work */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (!streq(times,"all") && !within(tsnap,times,0.0001))
            continue;                   /* skip work on this snapshot */
        if ( (bits & PhaseSpaceBit) == 0)
            continue;                   /* skip work, only diagnostics here */

        dprintf(0,"Converting time=%f to RV format",tsnap);
        
	fwrite(&nbody,sizeof(int),1,outstr);
        for (bp = btab; bp < btab+nbody; bp++) {
	    mass = Mass(bp);
            fwrite(&mass,sizeof(float),1,outstr);
        }
        t = tsnap;
        fwrite(&t,sizeof(float),1,outstr);
        for (bp = btab; bp < btab+nbody; bp++) {
	    for(j = 0; j < NDIM; j++) {
		pos[j] = Pos(bp)[j];
		vel[j] = Vel(bp)[j];
	    }
	    fwrite(pos,sizeof(float),NDIM,outstr);
	    fwrite(vel,sizeof(float),NDIM,outstr);
        }
    }
    strclose(instr);
    strclose(outstr);
}
