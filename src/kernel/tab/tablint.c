/*
 * TABLINT:	table syntax checker
 *
 *      13-jan-95  V0.1 prototype - Arie's initial idea
 *      16-feb-95  V0.1 added col*= keywords, out=, 
 *       8-dec-01  V0.2b   std. MAX_  
 *      18-may-05  V0.3 various to get it running for NEMO V3.x
 */
 
/**************** INCLUDE FILES ********************************/ 


#include <stdinc.h>
#include <getparam.h>
#include <extstring.h>
#include <strlib.h>
#include <ctype.h>

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {
    "in=???\n           input ASCII table filename",
    "lint=\n            lint specification file",
    "maxcol=\n          override max number of columns to check",
    "coltype=\n         * set/override type of columns",
    "colrange=\n        * set/override ranges in numeric columns",
    "colblank=\n        * set/override blank value for columns",
    "colmatch=\n        * set/override matches of ascii columns",
    "collength=\n       * set/override length of columns",
    "linelen=\n         override max linelength of a row",
    "out=\n             * optional output file of selected rows",
    "VERSION=0.3\n      18-may-05 PJT",
    NULL
};

string usage = "table syntax checker";



/**************** GLOBAL VARIABLES *******************************/

string infile;
stream instr, linstr;			/* file streams */
stream outstr = NULL;


#ifndef MAX_COL
#define MAX_COL          256             /* MAXIMUM number of columns */
#endif

#ifndef MAX_LINELEN
#define MAX_LINELEN  2048
#endif

#define MNEWDAT          80		/* space needed for one number */

#define COL_NONE    0           /* valid column types */
#define COL_CHAR    1
#define COL_INT     2
#define COL_REAL    3
#define COL_ALPHA  11
#define COL_UPPER  12
#define COL_LOWER  13
#define COL_DIGIT  14
#define COL_ALNUM  15

#define CHECK_LEN   0x01        /* checker options applied to a column */
#define CHECK_MIN   0x02
#define CHECK_MAX   0x04
#define CHECK_BLN   0x08


typedef struct col {            /* definition of a column */

    int nr;                     /* column nr (1..maxcol) */
    string name;                /* column name (optional) */
    int type;                   /* 'enumerated' type of the column */
    int check;                  /* bitmap of checking options */

    int length[2];              /* min and max length of column */

    union range_tag {           /* allowed min and max value */
        int irange[2];          /* separately checked with */
        float frange[2];        /* CHECK_MIN and CHECK_MAX */
    } range;

    union blank_tag {           /* column value that represents a blank */
        char *cblank;           /* checked if CHECK_BLN is set */
        int iblank;
        float fblank;
    } blank;

    string match;               /* match string */

} col, *colptr;


int maxcol, ncol=0;
col *cols = NULL;
int nlines = 0;
int linelen = MAX_LINELEN;

extern string *burststring(string, string);    

local char *skipdigit(char *);
local int bad_int(string), bad_real(string);
local int check_col(col *, string , int, string);
local col *get_col(int);
/****************************** START OF PROGRAM **********************/

void nemo_main()
{
    int linelen;

    maxcol = (hasvalue("maxcol") ? getiparam("maxcol") : -1);

    if (hasvalue("lint")) {
        linstr = stropen(getparam("lint"),"r"); 
        read_lintfile(linstr);
    }
    read_col_param("type");
    read_col_param("range");
    read_col_param("blank");
    read_col_param("match");
    read_col_param("length");

    if (hasvalue("out")) outstr = stropen(getparam("out"),"w");
    if (hasvalue("linelen")) linelen = getiparam("linelen");
    else linelen = MAX_LINELEN;

    infile = getparam("in");
    instr  = stropen(infile,"r");               
    lint(instr,linelen,outstr);
}


read_lintfile(stream instr)
{
    char   line[MAX_LINELEN];          /* input linelength */
    double dval[MAX_COL];            /* number of items (values on line) */
    char   newdat[MNEWDAT];         /* to store new column in ascii */
    int    nval, n, i, icol;
    string *sp;
        
    nlines=0;                       /* count lines read so far */

    for(;;) {                       /* loop over all lines in lintfile */

        if(fgets(line,MAX_LINELEN,instr)==NULL)
            break;
        nlines++;
        dprintf(3,"LINE: (%s)",line);
	if(line[0] == '#') continue;		/* don't use comment lines */

        sp = burststring(line," \t\n");
        n = xstrlen(sp, sizeof(string)) - 1;
        if (n==0) continue;
        if (streq(sp[0],"maxcol")) {            /* set max # columns allowed */
            if (maxcol>0) 
                warning("Cannot reset maxcol=%d, set at %d",
                        atoi(sp[1]),maxcol);
            else
                maxcol = atoi(sp[1]);
        } else if (streq(sp[0],"col")) {        /* parse a column descriptor */
            icol = atoi(sp[1]);
            if (icol<1 || (maxcol>0 && icol>maxcol))
                error("Illegal icol=%d for maxcol=%d",icol,maxcol);
            read_col(instr,icol);
        } else if (streq(sp[0],"row")) {        /* parse a row descriptor */
            read_row(instr);
        } else
            warning("Ignoring lint directive %s",sp[0]);
    } /* for(;;) */
    dprintf(0,"Read %d lines from lint file\n",nlines);
}

read_col(stream instr, int icol)
{
    char   line[MAX_LINELEN];          /* input linelength */
    string *sp;
    int n;
    col *cp;

    cp = get_col(icol);

    for(;;) {                       /* loop over all lines in file(s) */

        if(fgets(line,MAX_LINELEN,instr)==NULL) {
            error("Bad end reading col %d",icol);
            break;
        }
            
        dprintf(3,"LINE: (%s)",line);
	if(line[0] == '#') continue;		/* don't use comment lines */
        nlines++;

        sp = burststring(line," \t\n");
        n = xstrlen(sp, sizeof(string)) - 1;
        if (n==0) continue;
        if (streq(sp[0],"}"))
            break;
        else if (n<2)
            error("not a valid column command");

        set_column(cp,sp[0],&sp[1]);
    }
}

read_col_param(string cmd)
{
    char key[32];
    string *sp;
    int i, n, icol;
    col *cp;

    sprintf(key,"col%s",cmd);
    if (hasvalue(key)) {
        sp = burststring(getparam(key),",");
        n = xstrlen(sp,sizeof(string))-1;
        if (n%2) error("Wrong argument count in %s=%s",key,getparam(key));
        for (i=0; i<n; i += 2) {
            if (bad_int(sp[i])) 
                error("(%s) syntax error column nr=%s",key,sp[i]);
            icol = atoi(sp[i]);
            cp = get_col(icol);
            set_column(cp, cmd, &sp[1]);
        }
    }
}


col *get_col(int icol)
{
    col *cp;
    int i;

    if (cols==NULL) {
        if (maxcol<0) maxcol = MAX_COL;
        dprintf(1,"Allocating space for %d columns\n",maxcol);
        cols = (col *) allocate(maxcol*sizeof(col));
        for (i=0, cp=cols; i<maxcol; i++, cp++) {
            cp->nr = icol;
            cp->type = COL_NONE;
            cp->name = NULL;
        }
    }

    if (icol < 0 || icol > maxcol) error("Cannot handle column %d",icol);
    cp = &cols[icol-1];
    return cp;
}    


set_column(col *cp, string cmd, string *sp)
{
    if (streq(cmd,"name"))
        cp->name = scopy(sp[0]);
    else if (streq(cmd,"type")) {
        if (streq(sp[0],"int"))
            cp->type = COL_INT;
        else if (streq(sp[0],"real"))
            cp->type = COL_REAL;
        else if (streq(sp[0],"char"))
            cp->type = COL_CHAR;
        else if (streq(sp[0],"alpha"))
            cp->type = COL_ALPHA;
        else if (streq(sp[0],"upper"))
            cp->type = COL_UPPER;
        else if (streq(sp[0],"lower"))
            cp->type = COL_LOWER;
        else if (streq(sp[0],"digit"))
            cp->type = COL_DIGIT;
        else if (streq(sp[0],"alnum"))
            cp->type = COL_ALNUM;
        else 
            warning("Column (%d) type \"%s\" not understood",
                    cp->nr,sp[0]);
    } else if (streq(cmd,"range"))
        parse_range(cp,sp[0]);
    else if (streq(cmd,"blank"))
        parse_blank(cp,sp[0]);
    else if (streq(cmd,"length"))
        parse_length(cp,sp[0]);
    else if (streq(cmd,"match"))
        cp->match = scopy(sp[0]);
    else
        warning ("Ignoring column (%d) directive \"%s\"",
                    cp->nr,cmd);
}

read_row(stream instr)
{
    char   line[MAX_LINELEN];          /* input linelength */
    string *sp;
    int n;
    
    for(;;) {                       /* loop over all lines in file(s) */

        if(fgets(line,MAX_LINELEN,instr)==NULL) {
            warning("Bad end reading row");
            break;
        }
            
        dprintf(3,"LINE: (%s)",line);
	if(line[0] == '#') continue;		/* don't use comment lines */
        nlines++;

        sp = burststring(line," \n");
        n = xstrlen(sp, sizeof(string)) - 1;
        if (n==0) continue;
        if (streq(sp[0],"}"))
            break;
        warning ("Ignoring row directive: %s",line);
    }
}

parse_range(col *cp, string s)
{
    char *c = strchr(s,':');

    if (c==0) {
        warning("Must use ':' character in range command (column %d)",
            cp->nr);
        return;
    }
    if (!(cp->type == COL_INT || cp->type == COL_REAL)) {
        warning("Column (%d) type must be set before range defined",
            cp->nr);
        return;
    }

    cp->check = 0;

    c = s;
    if (*c!=':') {          /* minimum defined */
        cp->check |= CHECK_MIN;
        if (cp->type == COL_INT)
            cp->range.irange[0] = atoi(c);
        else if (cp->type == COL_REAL)
            cp->range.frange[0] = atof(c);
        c = strchr(c,':');
    } 
    c++;    
    if (*c) {               /* maximum defined */
        cp->check |= CHECK_MAX;
        if (cp->type == COL_INT)
            cp->range.irange[1] = atoi(c);
        else if (cp->type == COL_REAL)
            cp->range.frange[1] = atof(c);
    }

}

parse_blank(col *cp, string s)
{
    char *c = s;
}

parse_length(col *cp, string s)
{
    char *c = s;
}

/*
 * LINT: read lines from a file; each line can have up to maxline
 *       characters (excluding the newline?)
 *       This routine could also be used to check ascii fits tables
 *       that do not have newlines but fixed record length, in which
 *       maxline must be set to NAXIS1 of that table.
 */

lint( stream instr, int maxline, stream outstr)
{
  char   *line, msg[64];
  int    n, i, nlines, nilines, nolines;
  string *sp;
  col    *cp;    
  bool   doout, first = TRUE;
  
  line = (char *) allocate((maxline+1)*sizeof(char));
  
  nlines=0;               /* count lines read so far */
  nilines = nolines = 0;
  
  for(;;) {                   /* loop over all lines in file(s) */
    if (fgets(line,maxline+1,instr)==NULL) break;
    nlines++;
    dprintf(3,"LINE: %s",line);
    if(line[0] == '#' || line[0] == '\\' || line[0] == '|') {
      /* if (Qcomment && Qoutput) .... output .... */
      continue;		/* skip comment lines */
    }
    nilines++;
    sp = burststring(line," ,\t\n");
    n = xstrlen(sp,sizeof(string))-1;
    if (n==0) continue;                 /* blank line */
    
    if (first) {
      first = FALSE;
      if (cols==NULL) {
	if (ncol == 0) ncol = n;
      } else {
	if (ncol == 0) ncol = maxcol;       
      }
      dprintf(1,"Found %d columns\n",ncol);
    }
    if (cols) {
      doout = TRUE;            
      for (i=0, cp=cols; i<n; i++, cp++) {       /* check all columns */
	if (!check_col(cp, sp[i], nlines, msg)) {
	  dprintf(1,"# %s\n",msg);
	  doout=FALSE;
	}
      }
      if (doout) nolines++;
      if (outstr && doout) fputs(line,outstr);
    }
    if (n != ncol)
      warning("%s: line %d found %d/%d columns",infile,nlines,n,ncol);
  } /* for(;;) */
  dprintf(1,"Read %d lines from datafile\n",nlines);
  dprintf(1,"Read %d data lines from datafile\n",nilines);
  if (outstr)
    dprintf(1,"Wrote %d data lines to output\n",nolines);
  else
    dprintf(1,"Selected %d data lines\n",nolines);
}

int check_col(col *cp, string c, int linenr, char *msg)
{
    int ival;
    float fval;
    char *p = c;

    if (cp==0) return 0;

    msg[0] = 0;
    if (cp->type == COL_INT) {
        if (bad_int(c)) {
            strcat(msg,"bad syntax for integer");
            return 0;
        }
        if (cp->check) {
            ival = atoi(c);
            dprintf(1,"0x%x : checking %d to be in range %d : %d\n",
                cp->check, ival, cp->range.irange[0], cp->range.irange[1]);
            if (cp->check & CHECK_MIN && ival < cp->range.irange[0]) {
                strcat(msg,"below minimum");
                return 0;
            }
            if (cp->check & CHECK_MAX && ival > cp->range.irange[1]) {
                strcat(msg,"above maximum");
                return 0;
            }
        }
    } else if (cp->type == COL_REAL) {
        if (bad_real(c)) {
            strcat(msg,"bad syntax for real");
            return 0;
        }
        if (cp->check) {
            fval = atof(c);
            dprintf(1,"0x%x : checking %g to be in range %g : %g\n",
                cp->check, fval, cp->range.frange[0], cp->range.frange[1]);
            if (cp->check & CHECK_MIN && fval < cp->range.frange[0]) {
                strcat(msg,"below minimum");
                return 0;
            }
            if (cp->check & CHECK_MAX && fval > cp->range.frange[1]) {
                strcat(msg,"above maximum");
                return 0;
            }
        }
    } else if (cp->type == COL_ALPHA) {
        while (*p) if (!isalpha(*p++)) return 0;
    } else if (cp->type == COL_UPPER) {
        while (*p) if (!isupper(*p++)) return 0;
    } else if (cp->type == COL_LOWER) {
        while (*p) if (!islower(*p++)) return 0;
    } else if (cp->type == COL_DIGIT) {
        while (*p) if (!isdigit(*p++)) return 0;
    } else if (cp->type == COL_ALNUM) {
        while (*p) if (!isalnum(*p++)) return 0;
    }
    return 1;
}


/*
 *  Check if a string is a valid integer
 */


int bad_int(string word)
{
    char *cp = word;

    if (cp==NULL) return 1;     /* non-valid and null strings */
    if (*cp==0) return 1;       /* are illegal */

    if (*cp == '+' || *cp == '-') cp++;     /* sign */
    cp = skipdigit(cp);                     /* number */

    if (*cp) {
        dprintf(1,"bad_int: %s\n",cp);
        return 1;
    }
    return 0;
}

int bad_real(string word)
{
    char *cp = word;

    if (cp==NULL) return 1;     /* non-valid and null strings */
    if (*cp==0) return 1;       /* are illegal */

    if (*cp == '+' || *cp == '-') cp++;     /* sign */
    cp = skipdigit(cp);                     /* digits */
    if (*cp == '.') cp++;                   /* dot */
    cp = skipdigit(cp);                     /* digits */
    if (*cp && strchr("dDeE",*cp)) cp++;    /* exponent */
    if (*cp == '+' || *cp == '-') cp++;     /* sign */
    cp = skipdigit(cp);                     /* digits */

    if (*cp) {
        dprintf(1,"bad_real: %s\n",cp);
        return 1;
    } 
    return 0;
}


char *skipdigit(char *cp)
{
    while (*cp) {
        if (isdigit(*cp))
            cp++;
        else
            return cp;
    }
    return cp;

}
