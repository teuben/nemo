/*
 * TABDMS: convert column from/to HMS/DMS
 *
 *      24-jul-94  V0.1 created	PJT
 *	28-aug-98   0.2 ?
 *      13-sep-99   0.3 added hms PJT  - but not finished
 *      24-jan-00   0.4 added hms/dms with sign for decimal stuff
 *       8-dec-01   0.6 std MAX_ macros
 *      10-jan-03     a SGN -> SIGN
 *      31-jan-04   0.7 added fromhms=
 *
  
  TODO:  shuld use strptime() for parsing 

   %a, %A     Sun/Mon/...  or Sunday/...
   %b, %B     Jan/...	      January/...
   %d         day of month (01..31)
   %D         mm/dd/yy
   %e         day of month, blank padded (1..31)
   %H         hour (00..23)
   %j         day of year (001..366)
   %k         hour (0..23)
   %m         month (01..12)
   %M         minute (00..59)
   %S         second (00..60)
   %T         time (hh:mm:ss)
   %y         last two digits of year (00..99)
   %Y         year (1970..)
 */

#include <stdinc.h>
#include <getparam.h>
#include <ctype.h>

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n           input file name(s)",
    "out=???\n          output file name",
    "todms=\n           list of columns (1..) to convert to dms (0=none)",
    "tohms=\n           list of columns (1..) to convert to hms (0=none)",
    "fromhms=\n         list of columns (1..) to convert from hms (0=none)",
    "fromdms=\n         list of columns (1..) to convert from dms (0=none)",
    "separator=:\n      separator between output D-M-S.S",
    "format=%g\n        output format",
    "scale=1\n          scale applied to degrees in the fromhms=/fromdms=",
    "VERSION=0.8\n      25-sep-2013 PJT",
    NULL
};

string usage = "Convert to HMS/DMS tables";

string cvsid="$Id$";


stream instr, outstr;			/* file streams */

string sep;                             /* separator */
string fmt;                             /* output format */
real   scale;                           /* output scale factor */

#ifndef MAX_COL
#define MAX_COL 256
#endif

#ifndef MAX_LINELEN
#define MAX_LINELEN  2048
#endif

bool   keepc[MAX_COL+1];                 /* columns to keep (t/f) in old fmt */
int    colmode[MAX_COL+1];		 /* 0, +/-24, +/- 360 */
int    fromhms[MAX_COL+1];               /* convert HMS to decimal (seconds)? */
int    fromdms[MAX_COL+1];               /* convert DMS to decimal (seconds)? */

extern string *burststring(string, string);

nemo_main()
{
    setparams();

    instr = stropen(getparam("in"),"r");
    outstr = stropen (getparam("out"),"w");

    convert(instr,outstr);
}

setparams()
{
    int    col[MAX_COL];                    /* columns to skip for output */
    int    ncol, i, j;

    for (i=0; i<MAX_COL; i++) {
        colmode[i+1] = 0;
	fromhms[i+1] = 0;
	fromdms[i+1] = 0;
    }

    if (hasvalue("todms")) {
    	ncol = nemoinpi(getparam("todms"),col,MAX_COL);
        if (ncol < 0) error("Error parsing todms=%s",getparam("todms"));
        for (i=0; i<ncol; i++) {
            j = ABS(col[i]);
            if (j>MAX_COL) error("Column %d too large, MAX_COL=%d",j,MAX_COL);
            if (colmode[j] || fromhms[j]) error("Column %d already specified",j);
            colmode[j] = ( col[i] > 0 ? 360 : -360);
        }
    }
    if (hasvalue("tohms")) {
    	ncol = nemoinpi(getparam("tohms"),col,MAX_COL);
        if (ncol < 0) error("Error parsing tohms=%s",getparam("tohms"));
        for (i=0; i<ncol; i++) {
            j = ABS(col[i]);
            if (j>MAX_COL) error("Column %d too large, MAX_COL=%d",j,MAX_COL);
            if (colmode[j] || fromhms[j]) error("Column %d already specified",j);
            colmode[j] = ( col[i] > 0 ? 24 : -24);
        }
    }
    if (hasvalue("fromhms")) {
    	ncol = nemoinpi(getparam("fromhms"),col,MAX_COL);
        if (ncol < 0) error("Error parsing fromhms=%s",getparam("fromhms"));
        for (i=0; i<ncol; i++) {
            j = ABS(col[i]);
            if (j>MAX_COL) error("Column %d too large, MAX_COL=%d",j,MAX_COL);
            if (colmode[j] || fromhms[j]) error("Column %d already specified",j);
	    fromhms[j] = 1;
        }
    }
    if (hasvalue("fromdms")) {
    	ncol = nemoinpi(getparam("fromdms"),col,MAX_COL);
        if (ncol < 0) error("Error parsing fromdms=%s",getparam("fromdms"));
        for (i=0; i<ncol; i++) {
            j = ABS(col[i]);
            if (j>MAX_COL) error("Column %d too large, MAX_COL=%d",j,MAX_COL);
            if (colmode[j] || fromdms[j]) error("Column %d already specified",j);
	    fromdms[j] = 1;
        }
    }
    sep = getparam("separator");
    fmt = getparam("format");
    scale = getrparam("scale");
}



convert(stream instr, stream outstr)
{
    char   line[MAX_LINELEN];          /* input linelength */
    double dval;        
    int    nval, i, nlines, sign, decimalval,nhms,ndms;
    string *outv, *hms, *dms;       /* pointer to vector of strings to write */
    string seps=", \t";             /* column separators  */
    real   dd, mm, ss;
        
    nlines=0;               /* count lines read so far */
    for(;;) {                       /* loop over all lines in file(s) */

        if (get_line(instr, line) < 0)
            return 0;
        dprintf(3,"LINE: (%s)\n",line);
	if(line[0] == '#') continue;		/* don't use comment lines */
        nlines++;


        outv = burststring(line,seps);
        i=0;
        while (outv[i]) {
            if (colmode[i+1] == 0) {        /* no special mode: just copy */
	      if (fromhms[i+1]) {
		hms = burststring(outv[i],sep);
		nhms = xstrlen(hms,sizeof(string))-1;
		dval = atof(hms[0])*3600;
		if (nhms > 1) dval += atof(hms[1])*60;
		if (nhms > 2) dval += atof(hms[2]);
		if (nhms > 3) error("HMS string %s too many %s",outv[i],sep);
                fprintf(outstr,fmt,dval/3600.0);
                fputs(" ",outstr);
		freestrings(hms);
	      } else if (fromdms[i+1]) {
		dms = burststring(outv[i],sep);
		ndms = xstrlen(dms,sizeof(string))-1;
		dval = atof(dms[0])*3600;
		if (ndms > 1) dval += atof(dms[1])*60;
		if (ndms > 2) dval += atof(dms[2]);
		if (ndms > 3) error("DMS string %s too many %s",outv[i],sep);
                fprintf(outstr,fmt,dval*15/3600.0);
                fputs(" ",outstr);
		freestrings(dms);
	      } else {
                fputs(outv[i],outstr);
                fputs(" ",outstr);
	      }
            } else {
                if (nemoinpd(outv[i],&dval,1)<0) 
                    error("syntax error decoding %s",outv[i]);
		if (colmode[i+1] < 0) {     
                    sign = SGN(dval);
		    dval = ABS(dval);       /* do we allow < 0 ??? */
                    decimalval = (int)floor(dval);
                    fprintf(outstr,"%d",decimalval*sign);
                    fprintf(outstr,"%s",sep);
                    dval -= decimalval;
                    dval *= ABS(colmode[i+1]);
		}
                sign = SGN(dval);
                dval = ABS(dval);
                dd = floor(dval);
                dval = (dval-dd)*60.0;
                mm = floor(dval);
                ss = (dval-mm)*60.0;
                fprintf(outstr,"%02d",(int)dd*sign);
                fprintf(outstr,"%s",sep);
                fprintf(outstr,"%02d",(int)mm);
                fprintf(outstr,"%s",sep);
                fprintf(outstr,"%06.3f ",ss);

            }
            i++;
        }
        fputs("\n",outstr);
	freestrings(outv);
    } /* for(;;) */
}
