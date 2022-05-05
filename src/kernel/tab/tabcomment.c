/*
 * 	Add comments to a table, or comment certain lines, or show just comments
 *
 *	1-aug-92  written		Peter Teuben
 *	6-aug-92  V1.1 can also delete comment lines	PJT
 *     18-jun-98  V1.1b   increased debug output level by 1 	PJT
 *      8-dec-01      c   MAX_LINELEN
 *      9-apr-22  V2.1 using new tablev2, no more MAX_LINELEN
 *
 *  Modes:
 *      - show table and commented comments (default)
 *      - show only table (delete=t)
 *      - show only raw comments (raw=t)
 */

#include <stdinc.h>
#include <getparam.h>
#include <ctype.h>
#include <table.h>

string defv[] = {
	"in=???\n		  input ascii table",
	"out=-n                   output ascii table",
#if 0
	"lines=\n		  Lines to comment (1..Nlines)",
#endif
	"alpha=t\n		  Comment lines that start with alpha?",
	"blank=t\n		  Comment blank lines?",
	"punct=t\n		  Comment lines seem punctiation?",
	"delete=f\n		  Delete the comment lines?",
	"raw=f\n                  Show only the raw comments?",
	"comment=#\n              The comment character!",
	"VERSION=2.1\n		  22-apr-2022 PJT",
	NULL,
};

string usage = "Add comments to a table, or comments certain lines, or show just comments";


void nemo_main()
{
    stream instr, outstr;
    table  *tptr;
    char   *cp;
    bool   Qalpha, Qblank, Qpunct, Qkeep, Qraw;
    int    nlines=0, nreal=0;
    string comment;

    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");
    tptr = table_open(instr,1);
    
    Qalpha = getbparam("alpha");
    Qblank = getbparam("blank");
    Qpunct = getbparam("punct");
    Qkeep = !getbparam("delete");
    Qraw   = getbparam("raw");
    comment = getparam("comment");
    if (Qraw) {
      Qkeep=FALSE;
      warning("new raw mode");
    }
    dprintf(0,"start ncols=%d nrows=%d\n",table_ncols(tptr),table_nrows(tptr));
    while ( (cp=table_line(tptr)) ) {
        nlines++;

        while (isspace(*cp))            /* skip whitespace */
            cp++;

        if (*cp == '\0' && Qblank) {        /* if blank, comment */
            if (Qkeep) fprintf(outstr,"%c %s",*comment,cp);
	    if (Qraw)  fprintf(outstr,"%s",cp);
            continue;
        }

        if (isalpha(*cp) && Qalpha) {        /* if alpha, comment */
            if (Qkeep) fprintf(outstr,"%c %s",*comment,cp);
	    if (Qraw)  fprintf(outstr,"%s",cp);
            continue;
        }

        if (ispunct(*cp) && Qpunct) {        /* if punct, comment */
            if (Qkeep) fprintf(outstr,"%c %s",*comment,cp);
	    if (Qraw)  fprintf(outstr,"%s",cp);
            continue;
        }

        nreal++;
	if (!Qraw)
	  fprintf(outstr,"%s",cp);          /* else, regular output */
    }
    dprintf(1,"%s: commented %d/%d lines\n",getparam("in"),nlines-nreal,nlines);
    table_close(tptr);
    strclose(instr);
    strclose(outstr);
}
