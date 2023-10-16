/*
 *  SNAPMNMX: find min and/or max of snapshot bodyvars'c
 *
 *	16-Apr-91	V1.0 created      		PJT
 *      13-may-91       V1.1 added time to list of options  PJT
 *	 6-nov-93	V1.2 moment, NEMO V2.			pjt
 *      17-aug-2022     V1.3 Allow out=
 */

#include <stdinc.h>
#include <getparam.h>
#include <moment.h>
#include <vectmath.h>		/* otherwise NDIM undefined */
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			Input file (snapshot)",
    "var=x,y,z\n                Variables to find min/max of",
    "mode=min,max\n             Modes: {time,min,max,mean,sigma}",
    "times=all\n                Times of snapshot",
    "format=%g\n                Format to print with",
    "out=\n                     If min or max given, output in a 1-body snapshot",
    "VERSION=1.3\n		17-aug-2022 pjt",
    NULL,
};

string usage = "find min and/or max of snapshot variables";

#define MAXOPT    50

void nemo_main()
{
    stream instr, tabstr, outstr;
    real   tsnap, ekin, etot, dr, r, rv, v, vr, vt, aux;
    real   varmin[MAXOPT], varmax[MAXOPT];
    real   var0[MAXOPT], var1[MAXOPT], var2[MAXOPT];
    Moment var[MAXOPT];
    string headline=NULL, options, times, mnmxmode;
    Body *btab = NULL, *bp, *bq;
    bool   Qmin, Qmax, Qmean, Qsig, Qtime, scanopt();
    int i, n, nbody, bits, nsep, isep, nopt, ParticlesBit, iout, imin, imax;
    char fmt[20],*pfmt;
    string *burststring(), *opt;
    rproc btrtrans(), fopt[MAXOPT], faux;
    bool Qout;

    ParticlesBit = (MassBit | PhaseSpaceBit | PotentialBit | AccelerationBit |
            AuxBit | KeyBit);
    instr = stropen(getparam("in"), "r");	/* open input file */
    mnmxmode= getparam("mode");

    opt = burststring(getparam("var"),", ");
    nopt = 0;					/* count options */
    while (opt[nopt]) {				/* scan through options */
        fopt[nopt] = btrtrans(opt[nopt]);
        nopt++;
        if (nopt==MAXOPT) {
            dprintf(1,"\n\nMaximum number of var's = %d exhausted\n",MAXOPT);
            break;
        }
    }
    dprintf(1,"var: \n");
    for (i=0; i<nopt; i++)
        dprintf(1,"%s ",opt[i]);
    dprintf(1,"\n");
    dprintf(1,"mode: %s\n",mnmxmode);

    if (hasvalue("out") && nopt==1) {
      Qout = TRUE;
      outstr = stropen(getparam("out"),"w");
    } else
      Qout = FALSE;

    Qmin = scanopt(mnmxmode,"min");
    Qmax = scanopt(mnmxmode,"max");
    Qmean = scanopt(mnmxmode,"mean");
    Qsig = scanopt(mnmxmode,"sigma");
    Qtime = scanopt(mnmxmode,"time");
    if (!Qmin && !Qmax && !Qmean && !Qsig && !Qtime) 
        error("No mode selected");
#if 0
    pfmt = getparam("tab");
    if (pfmt!=NULL && *pfmt!=NULL) {
	dprintf(1,"Saving table in %s\n",pfmt);
        tabstr = stropen(pfmt,"w");
    } else
#endif
      if (Qout)
	tabstr = stderr;
      else
        tabstr = stdout;

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
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (!streq(times,"all") && !within(tsnap,times,0.0001))
            continue;                   /* skip work on this snapshot */
        if ( (bits & ParticlesBit) == 0)
            continue;                   /* skip work, only diagnostics here */

         iout = imin = imax = -1;
	 for (bp = btab, i=0; bp < btab+nbody; bp++, i++) {
	   for (n=0; n<nopt; n++) {
	     aux = fopt[n](bp,tsnap,i);
	     if (i==0) ini_moment(&var[n],2,0);
	     accum_moment(&var[n], aux, 1.0);
	     if (Qout) {
	       if (aux==min_moment(&var[0])) imin=i;
	       if (aux==max_moment(&var[0])) imax=i;
	     }
	   }
	 }
	    
	 if (Qtime)
	   fprintf(tabstr,fmt,tsnap);
	 if (Qmin) {
	   if (Qout) iout = imin;
	   for (n=0; n<nopt; n++)
	     fprintf(tabstr,fmt,min_moment(&var[n]));
	 }
	 if (Qmax) {
	   if (Qout) iout = imax;
	   for (n=0; n<nopt; n++)
	     fprintf(tabstr,fmt,max_moment(&var[n]));
	 }
	 if (Qmean) {
	   for (n=0; n<nopt; n++)
	     fprintf(tabstr,fmt,mean_moment(&var[n]));
	 }
	 if (Qsig) {
	   for (n=0; n<nopt; n++)
	     fprintf(tabstr,fmt,sigma_moment(&var[n]));
	 }
	 if (imin < 0)
	   fprintf(tabstr,"\n");
	 else
	   fprintf(tabstr,"%d\n",iout);	      
	 fprintf(tabstr,"\n");
	 if (iout>=0) {
	   int ibody = 1;
	   Body *b1 = &btab[iout];
	   put_snap(outstr, &b1, &ibody, &tsnap, &bits);	      	      
	 }
	 
    } // for()
    strclose(instr);
}
