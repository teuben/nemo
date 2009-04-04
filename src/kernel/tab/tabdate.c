/*
 * TABDATE: convert a date/time column: wrapping strftime() and strptime()
 *
 *      1-feb-05    0.1 Created
 *      9-may-05    0.2 added time0= (tcol= still defaulted to 'all')
 *      3-apr-09    0.3 initialize the 'struct tm' structures better, using dcol= now
 */

#include <stdinc.h>
#include <getparam.h>
#include <extstring.h>
#include <ctype.h>
/* #define _XOPEN_SOURCE /* glibc2 needs this */
#include <time.h>

extern string *burststring(string, string);

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n           input file name(s)",
    "out=???\n          output file name",
    "format1=???\n      format to read with",
    "format2=???\n      format to write with",
    "dcol=0\n           Column for date conversion (0=whole line)",
    "time0=\n           If given, reference input times w.r.t this time (format2=%s)",
    "scale=1\n          Scale factor output seconds are divided by when time0 used" 
    "VERSION=0.3\n      4-apr-09 PJT",
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
int dcol;
real scale;

string time0 = 0;
bool need_time0 = FALSE;
bool use_scale = FALSE;

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
  scale = getrparam("scale");
  use_scale = (scale != 1.0);
  dcol = getiparam("dcol");
  
  instr = stropen(getparam("in"),"r");
  outstr = stropen (getparam("out"),"w");
  
  convert(instr,outstr);
}

convert(stream instr, stream outstr)
{
  char   line[MAX_LINELEN];          /* input linelength */
  int  nlines, nwords, nskip;
  long long nt0, nt;
  real dt;
  struct tm tm, tm0;
  string *bp;


  if (time0) {
    tm0.tm_sec= tm0.tm_min = tm0.tm_hour = 0;
    strptime(time0,format1,&tm0);
    strftime(line,MAX_LINELEN,"%s",&tm0);
    nt0 = atoll(line);
    dprintf(0,"Using time0=%lld sec since 1970.0 using %s on %s \n",
	    nt0,format1,time0);
  }

  nlines=0;               /* count lines read so far */
  nskip=0;
  for(;;) {                       /* loop over all lines in file(s) */
    if (get_line(instr, line) < 0)
      return 0;
    dprintf(3,"LINE: (%s)\n",line);
    if(line[0] == '#') continue;		/* don't use comment lines */
    nlines++;
    tm.tm_sec= tm.tm_min = tm.tm_hour = 0;

    if (dcol>0) {
      bp = burststring(line,", ");                /* split line in words */
      nwords = xstrlen(bp,sizeof(string))-1;
      if (nwords<dcol) {
	nskip++;
	continue;
      }
      strcpy(line,bp[dcol-1]);                    /* only pick out 1 column now */
      freestrings(bp);
    }

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
      if (use_scale) {
	dt = (nt-nt0)/scale;
	sprintf(line,"%g",dt);
      } else
	sprintf(line,"%lld",nt-nt0);
    }
    fputs(line,outstr);
    fputs("\n",outstr);
  } /* for(;;) */
  dprintf(1,"Processed %d lines\n",nlines);
  if (nskip) warning("%d/%d lines skipped because dcol too large",nskip,nlines);
}
