/*
 *  SNAPMASK:  mask out particles while copying snapshot data 
 *
 *	 18-Jun-87	created			                    PJT
 *	  9-Mar-88	V1.2 data-history added	                    PJT
 *	  7-Jun-88	V1.3 new filestruct			    PJT
 *       19-Aug-88      V1.4 experimenting with reentry get_snap    PJT
 *	 24-oct-88	V1.5 also copy or init Key-field; file=	    PJT
 *	  4-mar-89          a  :  new getsnap/putsnap macros	    PJT
 *	  6-apr-89	V1.6 skip snapshots without particles	    PJT
 *	 13-apr-89	    a  : calling nemoinpi to parse 'select' PJT
 *       18-nov-90      V1.7 helpvec                                PJT
 *       20-may-94      V1.8 - nemo V2.x ; also warn only once      pjt
 *	  8-aug-96      local var                                   pjt
 *       20-may-01      warning->dprintf
 *       10-jun-01      added \n
 *       25-may-02      V1.9 don't mess up when upper bound > nbody pjt
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>	
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#define REALLOC
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>


string defv[] = {		
    "in=???\n			  ascii input file name ",
    "out=???\n			  snapshot output file name ",
    "select=all\n		  select string ",
    "times=all\n		  times select string ",
    "keyfile=\n			  file with Key field to select from ",
    "keyoffset=0\n                offsets to be applied extra to outkey ",
    "VERSION=1.9\n		  25-may-02 PJT",
    NULL,
};

string usage = "mask out particles while copying snapshot data";

local string select_str, times;

local string headline;		       /* random text message */
local stream instr,outstr;	       /* file streams */

local int    nbody;                    /* number of particles in input snapshot */
local double tsnap;                    /* current time in snapshot */
local Body   *btab=NULL;	       /* pointer to input snapshot */

local bool   *Qsel=NULL;              /* pointer to select */
local bool   *Qself=NULL;             /* pointer to select from file */
local int    *sel=NULL;

#define TIMEFUZZ	0.0001	      /* tolerance in time comparisons */


nemo_main()
{
    char *fname;
    int  i, nbody_out, nbody_old, nbody_file, nbody_max, bits, offset, nret;
    double ekin;
    Body  *bp, *bq;
    stream fstr;
    bool first = TRUE;

    select_str=getparam("select");
    times=getparam("times");
    fname=getparam("keyfile");
    offset=getiparam("keyoffset");
    instr = stropen (getparam("in"), "r");
    outstr = stropen (getparam("out"), "w");

    if (*fname != 0) {    /* get initial selection from file */
        fstr = stropen(fname,"r");
        get_history(fstr);                      /* get history */
        while (get_tag_ok(instr, HeadlineTag))  /* skip headlines */
            headline = get_string(instr, HeadlineTag);
        if (!get_tag_ok(instr, SnapShotTag))    /* must be a snapshot */
            error("No snapshots in key file=");
        get_snap(fstr, &btab, &nbody_file, &tsnap, &bits);
        if ((bits & KeyBit) == 0)
            error("Key file= must contain key field");
        for (i=0, bp=btab, nbody_max = -1; i<nbody_file; i++, bp++)
	    if (Key(bp) > nbody_max)
                nbody_max = Key(bp);
        dprintf (1,"File %s has maximum key %d\n",fname,nbody_max);
        nbody_max++;
        Qself = (bool *) allocate(nbody_max * sizeof(bool));
        for (i=0; i<nbody_max; i++)
            Qself[i] = FALSE;
        for (i=0, bp=btab; i<nbody_file; i++, bp++)
            Qself[Key(bp)] = TRUE;
        strclose(fstr);
    }

    get_history(instr);
    put_history(outstr);
    for(;;) {			/* loop until done reading snapshots */
    	get_history(instr);
        while (get_tag_ok(instr, HeadlineTag))
        	headline = get_string(instr, HeadlineTag);
        if (!get_tag_ok(instr, SnapShotTag))
		break;			/* done with work in loop */

        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (bits == 0 || bits == 1)
            continue;                   /* no information to copy */

        if (Qsel==NULL || nbody > nbody_old) {       /* allocate & init space */
            if (Qsel==NULL) {
                sel = (int *) allocate(nbody*sizeof(int));
                Qsel = (bool *) allocate(nbody*sizeof(bool));
            } else {
                sel = (int *) reallocate(sel, nbody*sizeof(int));
                Qsel = (bool *) reallocate(Qsel, nbody*sizeof(bool));
            }
            if (Qself != NULL && nbody > nbody_max) {
                Qself = (bool *) reallocate(Qself, nbody*sizeof(bool));
                for (i=nbody_max; i<nbody; i++)
                    Qself[i] = FALSE;               /* init tail to false */
                nbody_max = nbody;
            }
            if (!streq("all",select_str)) {
                for (i=0; i<nbody; i++) {
                    Qsel[i] = FALSE;		/* default: don't use them all */
        	    sel[i]= -1;
                }
		nret = nemoinpi(select_str,sel,nbody);
                for (i=0; i<nret; i++)
                    Qsel[sel[i]]=TRUE;
            } else				/* all are selected */
                for (i=0; i<nbody; i++)
                    Qsel[i] = TRUE;
            if (Qself != NULL)            /* also mask with keys from file */
                for (i=0; i<nbody; i++)
                    Qsel[i] &= Qself[i];                
         }
        if (headline != NULL) {      /* see if special output to be written */
            put_string(outstr, HeadlineTag, headline);
            headline = NULL;
        }
        if ((bits & TimeBit) == 0)
        	tsnap = 0.0;
        else if (!streq(times,"all") && !within(tsnap, times, TIMEFUZZ))
        	continue;		/* however skip this snapshot */
        dprintf (1,"Snapshot time=%f masked\n",tsnap);
	if ((bits & KeyBit) == 0) {		/* if Key not set */
	    bits |= KeyBit;			/* set it for output */
            if (first) dprintf(0,"Keys field initialized to order in file\n");
	    for (i=0, bp=btab; i<nbody; i++, bp++)
                Key(bp) = i;			/* and initialize Key's */
        }
        for (i=0,bp=btab,bq=btab; bp < btab+nbody; i++,bp++) {
            Key(bp) += offset;  /* apply extra offset */
	    if (!Qsel[i])
	        continue;       /* no copy */
	    if (bp==bq) {
	        bq++;
	        continue;       /* no need to copy yet, still in sync */
	    }
            bcopy(bp, bq, sizeof(Body));
            bq++;
        }
        nbody_out = bq - btab;
/*      bits = bits & (TimeBit | MassBit | PhaseSpaceBit | KeyBit);	 */
        put_snap(outstr, &btab, &nbody_out, &tsnap, &bits);
        nbody_old = nbody;   /* remember last allocated space */
        first = FALSE;
    }
    
    strclose (instr);
    strclose (outstr);
}
