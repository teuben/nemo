/*
 *    unwrap a table coming out of Paxton's Fortran-90 code
 *
 *     21-oct-04  written for ADASS poster		Peter Teuben
 */

#include <stdinc.h>
#include <getparam.h>
#include <ctype.h>

string defv[] = {
	"in=???\n		  input ascii table",
	"out=???\n		  output ascii table",
	"VERSION=1.0\n		  21-oct-04 PJT",
	NULL,
};

string usage = "Unwrap a f90 table";

#ifndef MAX_LINELEN
#define MAX_LINELEN  128
#endif

nemo_main()
{
    stream instr, outstr;
    char   line[MAX_LINELEN], *cp;
    bool   Qalpha, Qblank, Qpunct, Qkeep;
    int    len, nlines=0, nreal=0;
    string comment;

    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");
    while (fgets(line,MAX_LINELEN,instr) != NULL) {    /* loop all lines */
        nlines++;
	cp = line;
	len = strlen(line);
	line[len-2] = 0;     /* remove newline */
	cp++;                /* skip first space */
        fprintf(outstr,"%s",cp);           /* buffer out */
	if (len < 81) fprintf(outstr,"\n");   /* and add newline when this line was not full */
    }
}
