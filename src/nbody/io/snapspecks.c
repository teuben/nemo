/*
 *  SNAPSPECKS: tabulate a snapshot in specks format
 *
 *	7-nov-00    PJT     Created, Q&D
 *     11-mar-01    PJT     time floating
 *     19-may-01	    integer again?
 *     22-jun-07    PJT     moved into NEMO
 */

#include <nemo.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			Input file (snapshot)",
    "options=x,y,z,m\n            Output variables",
    "format=%g\n                Format for output",
    "times=all\n		Times to select snapshot",
    "VERSION=1.1\n	 	22-jun-07 PJT",
    NULL,
};

string usage="tabulate a snapshot in specks format for PartiView";


extern string *burststring(string,string);


#define MAXOPT    50

void nemo_main()
{
    stream instr, tabstr;
    real   tsnap, dr, aux;
    string times;
    Body *btab = NULL, *bp, *bq;
    bool   Qsepar, Qhead;
    int i, n, nbody, bits, nsep, isep, nopt, ParticlesBit, datatime = 0;
    char fmt[20],*pfmt;
    string *opt;
    rproc btrtrans(), fopt[MAXOPT];

    ParticlesBit = (MassBit | PhaseSpaceBit | PotentialBit | AccelerationBit |
            AuxBit | KeyBit);
    instr = stropen(getparam("in"), "r");	/* open input file */

    opt = burststring(getparam("options"),", ");
    nopt = 0;					/* count options */
    while (opt[nopt]) {				/* scan through options */
        fopt[nopt] = btrtrans(opt[nopt]);
        nopt++;
        if (nopt==MAXOPT) {
            warning("Maximum number of options = %d exhausted\n",MAXOPT);
            break;
        }
    }
    printf("##  file created by snapspecks\n");
    for (i=0; i<nopt; i++)
        printf("# %s ",opt[i]);
    printf("\n");

    times = getparam("times");
    pfmt = getparam("format");
    strcpy (fmt,pfmt);
    if (strchr(fmt,' ')==NULL && strchr(fmt,',')==NULL)
        strcat (fmt," ");       /* append blank if user did not specify sep */


    printf("datavar 0 lum\n");

    get_history(instr);                 /* read history */

    for(;;) {                /* repeating until first or all times are read */
	get_history(instr);
        if (!get_tag_ok(instr, SnapShotTag))
            break;                        
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (!streq(times,"all") && !within(tsnap,times,0.0001))
            continue;                   /* skip work on this snapshot */
        if ( (bits & ParticlesBit) == 0)
            continue;                   /* skip work, no data here */

        printf("# time=%g\n",tsnap);
        printf("datatime %d\n",datatime++);
        for (bp = btab, i=0; bp < btab+nbody; bp++, i++) {
                for (n=0; n<nopt; n++) {
                    aux = fopt[n](bp,tsnap,i);
                    printf(fmt,aux);
                }
                printf("\n");        
        }
    }
    strclose(instr);
}
