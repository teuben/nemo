/*
 *  SNAPTIPSY: convert a snapshot to tipsy binary format
 *
 *	14-Apr-94	V1.0 created      		TRQ
 *	23-dec-94	V1.1 added out=- and mode=	PJT
    -- not completed yet -- use only mode=d
 *	17-aug-00       V1.2 padding byte ....		PJT
 *      16-jun-11       V1.3 add warning - padding was turned off???   PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <math.h>
#include <stdlib.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
/*  #include <bodytrans.h>  CANNOT USE THIS */
#include "tipsydefs.h"

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			Input file (snapshot)",
    "out=-\n                    Output file (tipsy binary format)",
    "times=all\n		Times to select snapshot",
    "mode=dark\n                Output mode (dark|gas|star)",
    "swap=f\n                   Swap bytes on output?",
    "VERSION=1.3\n		16-jun-2011",
    NULL,
};

string usage="convert snapshot to tipsy binary format";

#define MAXOPT    50

#define MODE_DARK  1
#define MODE_GAS   2
#define MODE_STAR  3

extern string *burststring(string,string);
extern rproc btrtrans(string expr);    /* bodytrans.h : conflict with tipsy */

void nemo_main()
{
    stream instr, outstr;
    real   tsnap;
    string times, mode;
    bool Qswap;
    Body *btab = NULL, *bp;
    int i, ndim, nbody, bits, ParticlesBit, omode;
    struct dump header;
    struct dark_particle *dark;

    ParticlesBit = (MassBit | PhaseSpaceBit | PotentialBit | AccelerationBit |
            AuxBit | KeyBit);
    instr = stropen(getparam("in"), "r");	/* open input file */
    outstr = stropen(getparam("out"),"w");

    times = getparam("times");
    Qswap = getbparam("swap");
    mode = getparam("mode");
    switch (*mode) {
        case 'd':   omode = MODE_DARK; break;
        case 'g':   omode = MODE_GAS;  break;
        case 's':   omode = MODE_STAR; break;
        default:    error("mode=%s not supported",mode);
    }
    get_history(instr);                 /* read history */

    if (sizeof(header) != 32) 
      warning("%d: TIPSY header not 32 in size",sizeof(header));

    for(;;) {                /* repeating until first or all times are read */
	get_history(instr);
        if (!get_tag_ok(instr, SnapShotTag))
            break;                                  /* done with work */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (!streq(times,"all") && !within(tsnap,times,0.0001))
            continue;                   /* skip work on this snapshot */
        if ( (bits & ParticlesBit) == 0)
            continue;                   /* skip work, only diagnostics here */
	ndim = header.ndim = 3;
	header.nbodies = nbody;
	header.nsph = header.nstar = header.ndark = 0;
#ifdef TIPSY_NEEDPAD
	header.junk = 1;		/* padding byte */
#endif
	switch (omode) {
	    case MODE_DARK: 	header.ndark = nbody;  break;
	    case MODE_GAS: 	header.nsph  = nbody;  break;
	    case MODE_STAR: 	header.nstar = nbody;  break;
            default:    error("Badly parsed mode %s",mode);
        }
	header.time = tsnap;
        if (Qswap) {
            bswap(&header.time,    sizeof(double), 1);
            bswap(&header.nbodies, sizeof(int),    1);
            bswap(&header.ndim,    sizeof(int),    1);
            bswap(&header.nsph,    sizeof(int),    1);
            bswap(&header.ndark,   sizeof(int),    1);
            bswap(&header.nstar,   sizeof(int),    1);
        }
	fwrite((char *)&header,sizeof(header),1,outstr) ;
	dark = malloc(nbody*sizeof(*dark));	
	if(dark == NULL)
	   abort();
        for (bp = btab, i=0; bp < btab+nbody; bp++, i++) {
	    int j;

	    dark[i].mass = Mass(bp);
	    dark[i].phi = Phi(bp);
	    for(j = 0; j < ndim; j++) {
		dark[i].pos[j] = Pos(bp)[j];
		dark[i].vel[j] = Vel(bp)[j];
	    }
        }
        if (Qswap)
            bswap((char *)dark, sizeof(float), nbody*sizeof(*dark));
	fwrite((char *)dark, sizeof(*dark), nbody, outstr);
    }
    strclose(instr);
    strclose(outstr);
}
