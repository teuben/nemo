/*
 *  SNAPPRINT: tabulate a snapshot
 *
 *	19-Mar-88	V1.0 created      		PJT
 *	 6-jun-88	V1.2 new filestruct		PJT
 *	27-sep-88	V1.3 using 'options' keyword	PJT
 *      27-oct-88       V1.4 added r,vr,vt              PJT
 *			    a      aux
 *	 7-feb-89	V1.5 multiple snapshots if times= given	PJT  (bug)
 *	10-feb-89	V1.6 bodytrans() now allowed	PJT
 *	12-apr-89	V1.7 removed times= bug		PJT
 *	25-may-90	V1.8 tidiedup keywords,added tab	PJT
 *	15-nov-90	V1.9 got rid of old 1.5 code	PJT
 *	20-jan-94          a extra decl for solaris2    pjt
 *	24-mar-94	   b fix gcc warnings
 *	 3-apr-96	   c added usage line		pjt
 *	20-nov-96	V2.0 added header=		pjt
 *				for nbody,time
 *       9-oct-01          a  toying with time selection pjt
 *      31-dec-02       V2.1 gcc3/SINGLEPREC             pjt
 *       4-sep-03       V2.2 allow CSV output based      pjt
 *      24-feb-04       V2.4 add newline=t               pjt
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
    "options=x,y,z,vx,vy,vz\n	Things to output",
    "format=%g\n		Format used to output numbers",
    "separ=0\n			Special table of interparticle distances",
    "times=all\n		Times to select snapshot",
    "tab=\n			Standard output or table file?",
    "header=f\n			Add header (nbody,time)to output?",
    "newline=f\n                add newline in the header?",
    "csv=f\n                    Use Comma Separated Values format",
    "comment=f\n                Add table columns as common, instead of debug",
    "VERSION=2.4\n		24-feb-04 PJT",
    NULL,
};

string usage="tabulate a snapshot";

#define MAXOPT    50

extern string *burststring(string,string);

void nemo_main()
{
    stream instr, tabstr;
    real   tsnap, dr, aux;
    string times;
    Body *btab = NULL, *bp, *bq;
    bool   Qsepar, Qhead = getbparam("header");
    bool   Qcsv = getbparam("csv");
    bool   Qcomment = getbparam("comment");
    bool   Qnewline = getbparam("newline");
    int i, n, nbody, bits, nsep, isep, nopt, ParticlesBit;
    char fmt[20],*pfmt;
    string *opt;
    rproc_body fopt[MAXOPT];

    ParticlesBit = (MassBit | PhaseSpaceBit | PotentialBit | AccelerationBit |
            AuxBit | KeyBit | DensBit | EpsBit);
    instr = stropen(getparam("in"), "r");	/* open input file */

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
    if (Qcomment) {
      printf("# ");
      for (i=0; i<nopt; i++)
        printf("\t%s",opt[i]);
      printf("\n");
    } else {
      for (i=0; i<nopt; i++)
        dprintf(0,"%s ",opt[i]);
      dprintf(0,"\n");
    }

    if (hasvalue("tab")) {
	pfmt = getparam("tab");
	dprintf(0,"Saving table in %s\n",pfmt);
        tabstr = stropen(pfmt,"w");
    } else
        tabstr = stdout;
    times = getparam("times");
    pfmt = getparam("format");
    strcpy (fmt,pfmt);
    if (!Qcsv && strchr(fmt,' ')==NULL && strchr(fmt,',')==NULL)
        strcat (fmt," ");       /* append blank if user did not specify sep */

    nsep = getiparam("separ");
    if (nsep) {
      dprintf(1,"Printing log10 of every %d-th PP distance\n",nsep);
      Qsepar=TRUE;
    } else
      Qsepar=FALSE;


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
	if (!Qsepar) {				/* printf options */
	    if (Qhead) {
	      fprintf(tabstr,"%d ",nbody);
	      if (Qnewline) fprintf(tabstr,"\n");
	      fprintf(tabstr,"%g\n",tsnap);
	    }
            for (bp = btab, i=0; bp < btab+nbody; bp++, i++) {
                for (n=0; n<nopt; n++) {
                    aux = fopt[n](bp,tsnap,i);
		    if (Qcsv && n>0) fprintf(tabstr,",");
                    fprintf(tabstr,fmt,aux);
                }
                fprintf(tabstr,"\n");
            }
        } else {
            isep=nsep;
            for (bp=btab+1; bp<btab+nbody; bp++)
                for (bq=btab; bq<bp; bq++) {
                    if (isep--)
			continue;
                    else
                        isep=nsep;
                    dr = 0;
                    for (i=0; i<NDIM; i++)
                    	dr += sqr(Pos(bp)[i]-Pos(bq)[i]);
                    printf("%f\n",0.5*log10(dr));
                }
        }
    
    }
    strclose(instr);
}
