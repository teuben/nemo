/*
 * TESTIO: Test various NEMO I/O and CPU 
 *	Not meant for export into $NEMOBIN
 *	
 *      16-Nov-88       V1.0 written	                    PJT
 *	15-sep-90	V1.1 CPU tester for Fio vs. Cio     PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <math.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

#include "get_snapshot.c"
#include "put_snapshot.c"

/*#include "fbodyio.c"	OLD OLD OLD */

string defv[] = {
    "in=???\n		 input file name",
    "out=\n		 output file name (optional)",
    "io=get_snapshot\n   mode of I/O {get_snapshot, get_snap}",
    "work=bodywork\n     which work routine?",
    "VERSION=1.1\n       16-sep-90  PJT",
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
    string mode_io, mode_work;
    real   mscale, pscale, tsnap;
    vector rscale, vscale;
    string outname, headline = NULL;
    string times;
    Body *btab = NULL, *bp;
    int i, nbody, bits, nrscale, nvscale, ndim, nr, mr=32;
    real retval[32];
    SS ss;
					/* init NEMO environment */
    initparam(argv, defv);
    ndim = NDIM;
					/* open files */
    instr = stropen(getparam("in"), "r");
    outname = getparam("out");
    if (*outname)
        outstr = stropen(outname,"w");
    else
        outstr = NULL;
    mode_io = getparam("io");
    mode_work = getparam("work");

  if (streq("get_snapshot",mode_io)) {
 
    ini_snapshot(&ss);
    get_history(instr);
    if (outstr) put_history(outstr);
    while (get_snapshot(instr,&ss)) {
        dprintf(0,"get_snapshot: Read snapshot at time=%f\n",ss.time);
        if (*mode_work) {
            bodywork_(&ss.nbody,&ndim,ss.mass,ss.phase,retval,&mr,&nr);
            printf("[work(%d)] ",nr);
            for (i=0; i<nr; i++) printf("%f ",retval[i]);
            printf("\n");
	}
        if (outstr) put_snapshot(outstr,&ss);
        get_history(instr); /* catch any sandwiched ones */
    }

  } else {

    get_history(instr);
    if (outstr) put_history(outstr);
    while(get_snap(instr, &btab, &nbody, &tsnap, &bits)) {
        dprintf(0,"get_snap: read snapshot at time=%f\n",tsnap);
	if (*mode_work) {
            bodywork(nbody,btab);
	}
        if (outstr) put_snap(outstr, &btab, &nbody, &tsnap, &bits);
        get_history(instr); /* catch any sandwiched ones */
    }
  }
}



ini_snapshot(ssp)
SS *ssp;
{
    ssp->nbody = 0;
    ssp->bits = 0;
    ssp->time = 0.0;
    ssp->mass = NULL;
    ssp->phase = NULL;
}

/*	only in THREEDIM */

#define MULVV(w,v,u)            /* MULtiply Vector by Vector */         \
{                                                                       \
    register real *_wp = (w), *_vp = (v), *_up = (u);                   \
    *_wp++ = (*_vp++) * (*_up++);                                       \
    *_wp++ = (*_vp++) * (*_up++);                                       \
    *_wp   = (*_vp  ) * (*_up);                                         \
}
		
		

bodywork(nbody,btab)
int nbody;
Body *btab;
{
    Body *bp;
    int i;
    vector sumpos, sumvel, tmp;
    real   summass;

    CLRV(sumpos);
    CLRV(sumvel);

    for (bp=btab; bp<btab+nbody; bp++) {
        SETV(tmp,Pos(bp));
        MULVS(tmp,tmp,Mass(bp));
	MULVV(tmp,tmp,tmp);
        INCADDV(sumpos,tmp);

        SETV(tmp,Vel(bp));
        MULVS(tmp,tmp,Mass(bp));
	MULVV(tmp,tmp,tmp);
        INCADDV(sumvel,tmp);
	
        summass += Mass(bp);
    }
    for (i=0; i<NDIM; i++) {
        sumpos[i] = sqrt(sumpos[i])/summass;
        sumvel[i] = sqrt(sumvel[i])/summass;
    }
    printf("%f %f %f %f %f %f\n",
        sumpos[0],sumpos[1],sumpos[2],
        sumvel[0],sumvel[1],sumvel[2]);
}
