/*
 *  SNAPATLAS: tabulate a snapshot in atlas(5NEMO) format
 *
 *       8-apr-03       V1.0 created                    PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <math.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <bodytransc.h>

string defv[] = {
    "in=???\n			Input file (snapshot)",
    "times=all\n		Times to select snapshot",
    "options=i,m,x,y,z,vx,vy,vz,phi,aux,key\n	Things to output",
    "format=%.16e\n		Format used to output floating point numbers",
    "header=f\n			Add header to output?",
    "VERSION=1.0\n		31-dec-02 PJT",
    NULL,
};

string usage="tabulate a snapshot in ASCII atlas format";

#define MAXOPT    50

extern string *burststring(string,string);

void nemo_main()
{
    stream instr, tabstr = stdout;
    real   tsnap, dr, aux;
    string times;
    Body *btab = NULL, *bp, *bq;
    bool   Qhead, Qfirst = TRUE;
    int i, n, nbody, bits, nsep, isep, nopt, ParticlesBit;
    char fmt[20],*pfmt;
    string *opt;
    rproc_body fopt[MAXOPT];

#if defined(SINGLEPREC)
    warning("SINGLEPREC mode: data will not be in full precision");
#endif

    ParticlesBit = (MassBit | PhaseSpaceBit | PotentialBit | AccelerationBit |
            AuxBit | KeyBit | DensBit | EpsBit);
    instr = stropen(getparam("in"), "r");	/* open input file */
    Qhead = getbparam("header");

    opt = burststring(getparam("options"),", ");
    nopt = 0;					/* count options */
    while (opt[nopt]) {				/* scan through options */
        fopt[nopt] = btrtrans(opt[nopt]);
        nopt++;
        if (nopt==MAXOPT) {
            dprintf(0,"\n\nMaximum number of options = %d exhausted\n",MAXOPT);
            break;
        }
    }
    if (Qhead) {
	fprintf(tabstr,"# ");
	for (i=0; i<nopt; i++)
	  fprintf(tabstr,"%s ",opt[i]);
	fprintf(tabstr,"\n");
    }

    times = getparam("times");
    pfmt = getparam("format");
    strcpy (fmt,pfmt);
    if (strchr(fmt,' ')==NULL && strchr(fmt,',')==NULL)
        strcat (fmt," ");       /* append blank if user did not specify sep */

    get_history(instr);                 /* read history */
    for(;;) {                /* repeating until first or all times are read */
	get_history(instr);
        if (!get_tag_ok(instr, SnapShotTag))
            break;                                  /* done with work */
#if 0
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (!streq(times,"all") && !within(tsnap,times,0.0001))
            continue;                   /* skip work on this snapshot */
#else
        get_snap_by_t(instr, &btab, &nbody, &tsnap, &bits, times);	
#endif
        if ( (bits & ParticlesBit) == 0)
            continue;                   /* skip work, only diagnostics here */

	if (Qfirst)                     /* add blank line between snapshots */
	  Qfirst = FALSE;
	else
	  fprintf(tabstr," \n");
	fprintf(tabstr,"%d\n",nbody);
	fprintf(tabstr,fmt,tsnap);
	fprintf(tabstr,"\n");
	for (bp = btab, i=0; bp < btab+nbody; bp++, i++) {
	  for (n=0; n<nopt; n++) {
	    aux = fopt[n](bp,tsnap,i);
	    fprintf(tabstr,fmt,aux);
	  }
	  fprintf(tabstr,"\n");
	}
    }
    strclose(instr);
}
