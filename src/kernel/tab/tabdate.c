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
    "tcol=\n            Columns for date (default: all)",
    "time0=\n           If given, reference input times w.r.t this time (format2=%s)",
    "VERSION=0.2\n      9-may-05 PJT",
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

string time0 = 0;
bool need_time0 = FALSE;

nemo_main()
{
  format1 = getparam("format1");
  if (hasvalue("time0"))  {
    if (streq(getparam("time0"),"0"))
      need_time0 = TRUE;
    else 
      time0 = getparam("time0");
    format2 = strdup("%s");
  } else
    format2 = getparam("format2");

  
  instr = stropen(getparam("in"),"r");
  outstr = stropen (getparam("out"),"w");
  
  convert(instr,outstr);
}

convert(stream instr, stream outstr)
{
  char   line[MAX_LINELEN];          /* input linelength */
  int  nlines;
  long long nt0, nt;
  struct tm tm, tm0;

  if (time0) {
    strptime(time0,format1,&tm0);
    strftime(line,MAX_LINELEN,"%s",&tm0);
    nt0 = atoll(line);
    dprintf(0,"Command line uses: time0=%lld sec since 1970.0\n",nt0);
  }

  nlines=0;               /* count lines read so far */
  for(;;) {                       /* loop over all lines in file(s) */
    if (get_line(instr, line) < 0)
      return 0;
    dprintf(3,"LINE: (%s)\n",line);
    if(line[0] == '#') continue;		/* don't use comment lines */
    nlines++;
    strptime(line,format1,&tm);
    if (need_time0) {
      dprintf(0,"First line: %s\n",line);
      strftime(line,MAX_LINELEN,"%s",&tm);
      time0 = strdup(line);
      nt0 = atoll(time0);
      dprintf(0,"First line uses: time0=%lld sec since 1970.0\n",nt0);
      need_time0 = FALSE;
    }
    strftime(line,MAX_LINELEN,format2,&tm);
    if (time0) {
      nt = atoll(line);
      sprintf(line,"%lld",nt-nt0);
    }
    fputs(line,outstr);
    fputs("\n",outstr);
  } /* for(;;) */
  dprintf(1,"Processed %d lines\n",nlines);
}
