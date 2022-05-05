/*
 *    unwrap a table coming out of Paxton's Fortran-90 code
 *
 *     21-oct-04    written for ADASS poster		Peter Teuben
 *      4-mar-2022  converted for new tableV2 - but what was this for?
 */

#include <stdinc.h>
#include <getparam.h>
#include <ctype.h>
#include <table.h>

string defv[] = {
	"in=???\n		  input ascii table",
	"out=???\n		  output ascii table",
	"VERSION=1.1\n		  3-mar-2022 PJT",
	NULL,
};

string usage = "Unwrap a f90 table";

#ifndef MAX_LINELEN
#define MAX_LINELEN  128
#endif

void nemo_main()
{
    stream instr, outstr;
    table  *tp;
    string s;
    char   *cp;
    int    len, nlines=0;
    
    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");
    tp = table_open(instr,0);
    
    while ( (s=table_line(tp)) != NULL) {    /* loop all lines */
      nlines++;
      cp = s;
      len = strlen(s);
      s[len-2] = 0;     /* remove newline ??? */
      cp++;             /* skip first space */
      fprintf(outstr,"%s",cp);           /* buffer out */
      if (len < 81) fprintf(outstr,"\n");   /* and add newline when this line was not full */
    }
    table_close(tp);
}
