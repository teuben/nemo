/*
 * TABDATE: convert a date/time column: wrapping strftime() and strptime()
 *
 *      1-feb-05    0.1 Created
 */

#include <stdinc.h>
#include <getparam.h>
#include <ctype.h>
#include <time.h>

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n           input file name(s)",
    "out=???\n          output file name",
    "format1=???\n      format to read with",
    "format2=???\n      format to write with",
    "VERSION=0.1\n      1-feb-05 PJT",
    NULL
};

string usage = "Convert a date/time column";

string cvsid = "$Id$";



stream instr, outstr;			/* file streams */

#ifndef MAX_LINELEN
#define MAX_LINELEN  2048
#endif

string format1;
string format2;


nemo_main()
{
    format1 = getparam("format1");
    format2 = getparam("format2");

    instr = stropen(getparam("in"),"r");
    outstr = stropen (getparam("out"),"w");

    convert(instr,outstr);
}
convert(stream instr, stream outstr)
{
  char   line[MAX_LINELEN];          /* input linelength */
  int  nlines;
  struct tm tm;
        
  nlines=0;               /* count lines read so far */
  for(;;) {                       /* loop over all lines in file(s) */
    if (get_line(instr, line) < 0)
      return 0;
    dprintf(3,"LINE: (%s)\n",line);
    if(line[0] == '#') continue;		/* don't use comment lines */
    nlines++;
    strptime(line,format1,&tm);
    strftime(line,MAX_LINELEN,format2,&tm);
    fputs(line,outstr);
    fputs("\n",outstr);
  } /* for(;;) */
}
