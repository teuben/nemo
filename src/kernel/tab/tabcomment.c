/*
 * 	Add comments to a table, or comment certain lines
 *
 *	1-aug-92  written		Peter Teuben
 *	6-aug-92  V1.1 can also delete comment lines	PJT
 *     18-jun-98  V1.1b   increased debug output level by 1 	PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <ctype.h>

string defv[] = {
	"in=???\n		  input ascii table",
	"out=???\n		  output ascii table",
#if 0
	"lines=\n		  Lines to comment (1..Nlines)",
#endif
	"alpha=t\n		  Comment lines that start with alpha?",
	"blank=t\n		  Comment blank lines?",
	"punct=t\n		  Comment lines seem punctiation?",
	"delete=f\n		  Delete those comment lines?",
	"VERSION=1.1a\n		  18-jun-98 PJT",
	NULL,
};

string usage = "Add comments to a table, or comments certain lines";


nemo_main()
{
    stream instr, outstr;
    char   line[2048], *cp;
    bool   Qalpha, Qblank, Qpunct, Qkeep;
    int    nlines=0, nreal=0;

    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");
    Qalpha = getbparam("alpha");
    Qblank = getbparam("blank");
    Qpunct = getbparam("punct");
    Qkeep = !getbparam("delete");
    while (fgets(line,2048,instr) != NULL) {    /* loop all lines */
        nlines++;
	cp = line;

        while (isspace(*cp))            /* skip whitespace */
            cp++;

        if (*cp == '\0' && Qblank) {        /* if blank, comment */
            if (Qkeep) fprintf(outstr,"# %s",line);
            continue;
        }

        if (isalpha(*cp) && Qalpha) {        /* if alpha, comment */
            if (Qkeep) fprintf(outstr,"# %s",line);
            continue;
        }

        if (ispunct(*cp) && Qpunct) {        /* if punct, comment */
            if (Qkeep) fprintf(outstr,"# %s",line);
            continue;
        }

        nreal++;
        fprintf(outstr,"%s",line);          /* else, regular output */
    }
    dprintf(1,"%s: commented %d/%d lines\n",getparam("in"),nlines-nreal,nlines);
}
