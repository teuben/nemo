/*
 * SNAPTRIM.C: cut a snapshot file down to size.
 *
 *	V1.3: 12-apr-89		allow early retirement (checkall=false)	PJT
 *	V1.4: 12-may-89		interactive
 *	V1.5  10-dec-90		fixed up for new release
                    ============ HAS TO BE FIXED UP WHEN SETPARAM DONE ===
 *	   a  31-oct-91    	set TIMEFUZZ to 0.00001
 *	   b   7-oct-92		warning if nothing ever output
 *	   c   4-feb-93		usage
 *	   d   9-aug-95         fixed beyond() for new # syntax in times=
 *	   e  16-feb-97		SINGLEPREC support
 *      V1.6   5-mar-98         supporting time=first and time=last
 *	      15-jun-02         debug output
 */

/* #define INTERACT */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <snapshot/snapshot.h>

string defv[] = {               /* DEFAULT INPUT PARAMETERS		    */
    "in=???\n                     input file name",
    "out=???\n			  output file name",
    "times=all\n		  time range to scan through",
    "partcyc=1\n		  keep one in partcyc particle frames",
    "diagcyc=1\n		  keep one in diagcyc diagnostics",
    "amnesia=false\n		  if true, do not output history, etc",
    "checkall=false\n             must it check all snapshots",
#if defined(INTERACT)
    "more=y\n                     needs interactive SETPARAM part",
#endif
    "VERSION=1.6\n		  5-mar-98 PJT",
    NULL,
};

string usage="cut a snapshot file down to size";

#define TimeFuzz  0.00001	  /* slop tolerated in time comparisons */

extern bool within(real, string, real);
extern bool beyond(real, string, real);

nemo_main()
{
    stream instr, outstr;
    string times;
    int partcyc, diagcyc, npart, ndiag;
    bool pramflag, timeflag, partflag, diagflag;
    bool checkall, more, first=TRUE, something=FALSE;
    bool Qfirst=FALSE, Qlast=FALSE, Qdone=FALSE;
    real time;
    char c;

    instr = stropen(getparam("in"), "r");
    get_history(instr);
    outstr = stropen(getparam("out"), "w");
    if (! getbparam("amnesia"))
	put_history(outstr);
    times = getparam("times");
        Qfirst = streq(times,"first");
        Qlast = streq(times,"last");
    partcyc = getiparam("partcyc");
    diagcyc = getiparam("diagcyc");
    checkall = getbparam("checkall");
#if defined(INTERACT)
    more=getbparam("more");
#endif
    npart = ndiag = 0;
    for (;;) {
        get_history(instr);                                 /* skip junk */
        if (!get_tag_ok(instr, SnapShotTag))
            break;                                /* done with snapshots */
	get_set(instr, SnapShotTag);
	pramflag = timeflag = partflag = diagflag = FALSE;
	if (get_tag_ok(instr, ParametersTag)) {
	    pramflag = TRUE;
	    get_set(instr, ParametersTag);
	    if (get_tag_ok(instr, TimeTag)) {
		timeflag = TRUE;
		get_data_coerced(instr, TimeTag, RealType, &time, 0);
		dprintf(1,"Found Time=%g\n",time);
	    }
	    get_tes(instr, ParametersTag);
	}
	if (! timeflag ||
	      (streq(times, "all") || within(time, times, TimeFuzz)) ||
	      ( (first && Qfirst) || Qlast) ) {
	    if (get_tag_ok(instr, ParticlesTag)) {
		partflag = (partcyc != 0 && npart % partcyc == 0);
		npart++;
	    }
	    if (get_tag_ok(instr, DiagnosticsTag)) {
		diagflag = (diagcyc != 0 && ndiag % diagcyc == 0);
		ndiag++;
	    }
	} else if (timeflag && beyond(time,times,TimeFuzz))  /* all done? */
            if (!checkall) {
                 get_tes(instr,SnapShotTag);
                 break;            /* done with all - break reading loop */
            }
#if defined(INTERACT)
        printf("Time = %f ; partflag = %d diagflag = %d, read it? (y/n): ",
                time, partflag, diagflag);
#if defined(SETPARAM)
        if (Qinter) {
        } else {
        }
#else
        do {
	   fflush(stdin);
           c = (char )getchar();
        } while (c!='y' && c!='n');
        if (c=='n')
            partflag = diagflag = FALSE;
#endif
#endif
	if (partflag || diagflag) {
            something = TRUE;
            if (!Qlast)
	        dprintf(0,"time =%8.3f\tnpart =%4d\tndiag =%4d\toutputing %s\n",
		   timeflag ? time : 0.0, npart, ndiag,
		   partflag ? "particles" : "diagnostics");
	    if (Qlast) {
                rewind(outstr);
	        put_history(outstr);
            }
	    put_set(outstr, SnapShotTag);
	    if (pramflag)
		copy_item(outstr, instr, ParametersTag);
	    if (partflag)
		copy_item(outstr, instr, ParticlesTag);
	    if (diagflag)
		copy_item(outstr, instr, DiagnosticsTag);
	    put_tes(outstr, SnapShotTag);
            if (Qfirst) Qdone=TRUE;
	}
	get_tes(instr, SnapShotTag);
	fflush(outstr);
	if (Qdone) break;
    }
    strclose(instr);
    if (Qlast && something)
        dprintf(0,"time =%8.3f\tnpart =%4d\tndiag =%4d\toutputing %s\n",
                    timeflag ? time : 0.0, npart, ndiag,
                    partflag ? "particles" : "diagnostics");
    if (!something) warning("Nothing was ever written out");
    strclose(outstr);
}

/*
 * BEYOND.C: determine if a floating point number is beyond specified range,
 * represented as a string of the form "<subrange1>,<subrange2>,..." where
 * each <range> is either a single floating point number, or a pair of
 * numbers seperated by a ":".  To allow for small uncertainties in the
 * values tested, floating-point comparison is done with a specified
 * fuzzyness parameter.
 */


bool beyond(real val, string range, real fuzz)
{
    char *endptr, *subptr, *sepptr, *colptr;
    real sublow, subhi;

    if (*range == '#') return FALSE;

    endptr = range + strlen(range);		/* point to term. NULL */
    for (subptr = range; subptr != endptr; ) {	/* for each subrange */
        sepptr = strchr(subptr, ',');		/*   pnt to subrange end */
	if (sepptr == NULL)			/*   last subrange listed? */
	    sepptr = endptr;			/*     fix up subend ptr */
	colptr = strchr(subptr, ':');		/*   scan subrange for : */
	if (colptr > sepptr)			/*   in another subrange? */
	    colptr = NULL;			/*     then dont use it */
	sublow = atof(subptr) - fuzz/2.0;	/*   set low end of range */
	if (colptr != NULL)			/*   high end specified? */
	    subhi = atof(colptr+1) + fuzz/2.0;	/*     set high end */
	else
	    subhi = sublow + fuzz;		/*     just use low end */
        if (val > subhi)                        /*   beyond subrange ?  */
	    return (TRUE);
	subptr = sepptr;			/*   advance subrange ptr */
	if (*subptr == ',')			/*   more ranges to do? */
	    subptr++;				/*     move on to next */
    }
    return (FALSE);
}

