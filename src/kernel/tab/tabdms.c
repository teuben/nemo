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
 *      30-dec-21   0.9 use the new atox(), major fixes
 *      31-dec-21   1.0 implemented center=
 *
  
  TODO:  should use strptime() for parsing 

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


   TODO:  allow a special mode with a center, and all coordinates are then relative to
          the center value. Of course output can be in degrees, minutes are arcsec
 */

#include <stdinc.h>
#include <extstring.h>
#include <getparam.h>
#include <table.h>
#include <ctype.h>

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n           input file name",
    "out=???\n          output file name",
    "todms=\n           list of columns (1..) to convert to dms (0=none)",
    "tohms=\n           list of columns (1..) to convert to hms (0=none)",
    "fromdms=\n         list of columns (1..) to convert from dms (0=none)",
    "fromhms=\n         list of columns (1..) to convert from hms (0=none)",
    "separator=:\n      separator between output D M and S.S",
    "format=%g\n        output format",
    "scale=1\n          scale applied to degrees in the fromhms=/fromdms= and center=",
    "center=\n          If coordinates are to reference to a center",
    "VERSION=1.1\n      7-may-2022 PJT",
    NULL
};

string usage = "Convert to HMS/DMS";



stream instr, outstr;			/* file streams */
table *tptr;                            /* table */

string sep;                             /* separator */
string fmt;                             /* output format */
real   scale;                           /* output scale factor */
bool Qcenter;
real x_cen, y_cen, x_fac, cosdec;       /* center parameters */

#ifndef MAX_COL
#define MAX_COL 256
#endif

#ifndef MAX_LINELEN
#define MAX_LINELEN  2048
#endif

int    colmode[MAX_COL+1];		 /* 0, +/-24, +/- 360 for conversion to HMS/DMS  */
int    fromhms[MAX_COL+1];               /* convert HMS to decimal (seconds)? */
int    fromdms[MAX_COL+1];               /* convert DMS to decimal (seconds)? */
int    ncol_xy, col_x, col_y;


// local
void setparams();
void convert(stream instr, stream outstr);
void diff(stream instr, stream outstr);
void center(string cen);
  
void nemo_main()
{
    setparams();

    instr = stropen(getparam("in"),"r");
    outstr = stropen (getparam("out"),"w");

    if (Qcenter)
      diff(instr,outstr);
    else
      convert(instr,outstr);
}

void setparams()
{
    int    col[MAX_COL];                    /* columns to skip for output */
    int    ncol, i, j;

    Qcenter = hasvalue("center");
    ncol_xy = 0;
    x_fac = 1.0;

    for (i=0; i<MAX_COL; i++) {
        colmode[i+1] = 0;         // will be 0, 24 or 360
	fromhms[i+1] = 0;         // columns start at 1
	fromdms[i+1] = 0;         // 0 means not used here
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
	    if (Qcenter) {
	      if (ncol > 1) error("Cannot deal with more than one HMS column");
	      col_x = j;
	      ncol_xy++;
	      x_fac = 15.0;
	    }
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
	    if (Qcenter) {
	      if (ncol_xy == 2 && ncol > 1) error("Cannot deal with multiple DMS columns");
	      if (ncol_xy == 0)
		col_x = j;
	      else
		col_y = j;
	      ncol_xy++;
	    }
        }
    }

    if (Qcenter) center(getparam("center"));

    if (Qcenter && ncol_xy != 2) error("No fromhms= and fromdms= column selected");
    dprintf(1,"center cols: %d %d\n",col_x, col_y);
    
    // @todo   H:M:S is normal,
    //     NED  example 00h07m15.84s, +27d42m29.1s
    //     CASA example 12:22:55.212, +15.49.15.500    (who came up with THAT?)

    sep = getparam("separator");   
    fmt = getparam("format");
    scale = getrparam("scale");

}

void center(string scen)
{
  int ncen;
  double xycen[2];

  ncen = nemoinpx(scen, xycen, 2);
  if (ncen != 2) error("Error parsing %s", scen);
  x_cen = xycen[0] * x_fac;
  y_cen = xycen[1];
  cosdec = cos(y_cen * PI/180.0);
  dprintf(1,"center:  %g %g\n", x_cen, y_cen);
}

void diff(stream instr, stream outstr)
{
  //char   line[MAX_LINELEN];          /* input linelength */
  char *line;
  double xval, yval;
  int    nval, i, nlines;
  string *outv;                   /* pointer to vector of strings to write */
  string seps=", \t";       /* column separators  */

  tptr = table_open(instr, 1);
  
  nlines=0;               /* count lines read so far */
  for(;;) {                       /* loop over all lines in file(s) */
    line = table_line(tptr);
    if (line == NULL) return;
    dprintf(3,"LINE: (%s)\n",line);
    if(line[0] == '#') continue;		/* don't use comment lines */
    nlines++;

    outv = burststring(line,seps);
    xval = atox(outv[col_x-1]) * x_fac;
    yval = atox(outv[col_y-1]);
    printf("%s %s   %g %g\n",
	   outv[col_x-1], outv[col_y-1],
	   (xval-x_cen)*cosdec*scale, (yval-y_cen)*scale);
  } /* for(;;) */
}



void convert(stream instr, stream outstr)
{
  //char   line[MAX_LINELEN];          /* input linelength */
  char *line;
  double dval;        
  int    nval, i, nlines, sign, decimalval,nhms,ndms;
  string *outv, *hms, *dms;       /* pointer to vector of strings to write */
  string seps=", \t";             /* column separators  */
  real   dd, mm, ss;

  tptr = table_open(instr, 1);
  
  nlines=0;               /* count lines read so far */
  for(;;) {                       /* loop over all lines in file(s) */
    line = table_line(tptr);
    if (line==NULL) return;
    dprintf(3,"LINE: (%s)\n",line);
    if(line[0] == '#') continue;		/* don't use comment lines */
    nlines++;

    outv = burststring(line,seps);
    i=0; // column counter, internally 0-based, but 1-based where user interface was involved
    while (outv[i]) {
      if (colmode[i+1] == 0) {          // convert from sexagesimal, or pass through (fromhms/fromdms)
	if (fromhms[i+1]) {
	  dval = atox(outv[i]) * 15;
	  fprintf(outstr,fmt,dval);
	  fputs(" ",outstr);
	} else if (fromdms[i+1]) {
	  dval = atox(outv[i]);
	  fprintf(outstr,fmt,dval);
	  fputs(" ",outstr);
	} else {                        // pass through
	  fputs(outv[i],outstr);
	  fputs(" ",outstr);
	}
      } else {                          // todms/tohms:   convert degree to sexagesimal
	if (nemoinpd(outv[i],&dval,1)<0) 
 	  error("syntax error decoding %s",outv[i]);
	if (colmode[i+1] < 0) {
	  dprintf(1,"colmode %d: %d\n",i,colmode[i+1]);	     
	  sign = SGN(dval);
	  dval = ABS(dval);       /* do we allow < 0 ??? */
	  decimalval = (int)floor(dval);
	  fprintf(outstr,"%d",decimalval*sign);
	  fprintf(outstr,"%s",sep);
	  dval -= decimalval;
	  dval *= ABS(colmode[i+1]);
	}
	if (colmode[i+1] == 24)
	  scale = 15;
	else
	  scale = 1;
	sign = SGN(dval);
	dval = ABS(dval/scale);
	dd = floor(dval);
	dval = (dval-dd)*60.0;
	mm = floor(dval);
	ss = (dval-mm)*60.0;
	fprintf(outstr,"%02d",(int)dd*sign);
	fprintf(outstr,"%s",sep);
	fprintf(outstr,"%02d",(int)mm);
	fprintf(outstr,"%s",sep);
	fprintf(outstr,"%06.3f ",ss);
	// @todo     if the final string is 60.000 it should become 0.000
	//           and mm and dd might need to be incremented
	
      }
      i++;
    } // while
    fputs("\n",outstr);
    freestrings(outv);
  } /* for(;;) */
}
