/*
 * TSF: type the contents of a structured file in human-readable
 *      and possibly abbreviated form.
 * V 1.0: Joshua Barnes 1986	original program development
 *        Peter Teuben 7/87	casting pointers fixed for Ultrix
 * V 1.1  Joshua Barnes 7/87    max-precision option
 * V 1.X  Lyman Hurd 8/87	using deferred input
 * V 2.0  Joshua Barnes 4/88    new filestruct package
 *        Peter Teuben 3/90     made GCC happy by supplying silly return 
 *                                     codes in find_name, find_fmt]
 * V2.1   PJT 12-sep-90		added helpvec to defv
 * V2.2   PJT 5-oct-90		added dprintf() of some system parameters
 * V2.3   PJT 8-dec-90		minor improvements 
 * V2.4   PJT 10-feb-91         added item= keyword
 * V2.5   PJT 12-jun-91		nemo_main, to make it return proper exit()
 * V2.6   pjt 4-mar-94          ansi
 * V2.7   pjt 21-mar-01         XML output option, written in 
 *                              "Coffee Plantation", Tucson, AZ.
 * V3.0   pjt 14-jun-02         integer output not in octal anymore
 * V3.1   pjt 29-aug-02         confirmed bug item=XXX doesn't work anymore
 * V3.1a  pjt 25-nov-03         fix minor (but fatal) xml syntax error
 * V3.1e  pjt 10-aug-09         fix bug handled items > 2GB
 * V3.2   wd  08-sep-11         margin automatic adapts to terminal width
 *                              NOTE:  in a pipe or redirect this is not correct
 * V3.3   pjt 18-sep-2013       fixed the terminal width problem on pipes?
 *            23-oct-2013       another fix, using both stdin and stdout
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#ifdef unix
# include <sys/ioctl.h>
# include <stdio.h>
# include <unistd.h>
#endif

string defv[] = {
    "in=???\n                     input file name ",
    "maxprec=false\n		  print nums with max precision ",
    "maxline=4\n                  max lines per item ",
    "allline=false\n              print all lines (overrides maxline) ",
    "indent=2\n			  indentation of compound items ",
#ifdef unix
    "margin=\n                    righthand margin (default: window width)",
#else
    "margin=72\n		  righthand margin ",
#endif
    "item=\n                      Select specific item",
    "xml=f\n                      output data in XML format? (experimental)",
    "octal=f\n                    Force integer output in octal again?",
    "VERSIONf=3.3b\n		  24-oct-2013 PJT ",
    NULL,
};

string usage="type contents of a (binary) structured file";

string cvsid="$Id$";

local stream instr;			/* input stream from struct. file   */
local bool maxprec;			/* if true, use max precision forms */
local bool allline;			/* if true, disable maxline limit   */
local int maxline, indent, margin;      /* formating parameters             */
local string testtag;                   /* if non-zero, print only this tag */
local bool xml;                         /* output data in XML format        */
local bool octal;                       /* integer output in octal ?        */

/* local functions */
void   print_item   (string);
void   print_set    (string);
void   print_tes    (string);
void   print_header (string, string, int *);
void   print_data   (string, string, int *);
string find_name    (string);
string find_fmt     (string);
bool   outstr       (string);
void   end_line     (void);
void   change_indent(int);

#define BUFLEN   512


void nemo_main()
{
    string *tags;

    dprintf(2,"TSF: MaxSetLen = %d\n",MaxSetLen);   /* DEBUG - PJT */
    instr = stropen(getparam("in"), "r");
    maxprec = getbparam("maxprec");
    maxline = getiparam("maxline");
    indent = getiparam("indent");
#ifdef unix
    if(!hasvalue("margin")) {
      struct winsize w1,w2;
      ioctl(STDIN_FILENO, TIOCGWINSZ, &w1);
      ioctl(STDOUT_FILENO, TIOCGWINSZ, &w2);
      dprintf(1,"w_in=%d w_out=%d\n",w1.ws_col,w2.ws_col);
      if (w1.ws_col==0)
	margin = w2.ws_col;
      else if (w2.ws_col==0)
	margin = w1.ws_col;
      else
	margin = MIN(w1.ws_col,w2.ws_col);
    } else {
      margin = getiparam("margin");
    }
#else
    margin = getiparam("margin");
#endif
    if (margin<10) {
      warning("awfully small margin=%d",margin);
      margin=80;
    }
    if (margin>BUFLEN) {
      warning("BUFLEN %d might be too small (margin=%d)",BUFLEN,margin);
      margin=80;
    }
    dprintf(2,"margin=%d\n",margin);

    allline = getbparam("allline");
    if (hasvalue("item")) warning("item= is broken");
    testtag = getparam("item");
    xml = getbparam("xml");
    if (xml) {     /* if XML output used, it should print all items */
      allline = TRUE;
      /* @todo : should this be added:
         <?xml version='1.0' encoding='UTF-8'?> 
       */
      printf("<nemo>\n");
    }
    octal = getbparam("octal");
    while ((tags = list_tags(instr)) != NULL) {
        print_item(*tags);
	free(*tags);
	free((char *)tags);
    }
    if (xml) printf("</nemo>\n");
}

void print_item(string tag)
{
    string type, *tags, *tp;
    int *dims;

    type = get_type(instr, tag);
    if (streq(type, SetType)) {
	get_set(instr, tag);
	print_set(tag);
	tags = list_tags(instr);
	for (tp = tags; *tp != NULL; tp++)
	    print_item(*tp);
	get_tes(instr, tag);
	print_tes(tag);
	for (tp = tags; *tp != NULL; tp++)
	    free(*tp);
	free((char *)tags);
    } else {
	dims = get_dims(instr, tag);
	print_header(tag, type, dims);
	(void) outstr(" ");
	print_data(tag, type, dims);
	end_line();
	if (dims != NULL)
	    free((char *)dims);
    }
    free(type);
}

void print_set(string tag)
{
    char buf[BUFLEN];  /* danger: buffer overflow */

    if (xml)
      sprintf(buf, "<%s>", tag);
    else
      sprintf(buf, "set %s", tag);
    (void) outstr(buf);
    end_line();
    change_indent(indent);
}
    
void print_tes(string tag)
{
    char buf[BUFLEN];  /* danger: buffer overflow */

    change_indent(- indent);
    if (xml)
      sprintf(buf, "</%s>",tag);
    else
      sprintf(buf, "tes");
    (void) outstr(buf);
    end_line();
}

void print_header(string tag, string type, int *dims)
{
    char buf[BUFLEN];  /* danger: buffer overflow */
    int *dp;
    
    if (xml)
      sprintf(buf,"<%s type=\"%s\" ",tag,find_name(type));
    else
      sprintf(buf, "%s %s", find_name(type), tag);
    (void) outstr(buf);
    if (dims != NULL) {				/* is this a plural item?   */
      if (xml) {
	sprintf(buf,"dim=\"");
        outstr(buf);
      }
      for (dp = dims; *dp != 0; dp++) {	/*   loop over dimensions   */
	sprintf(buf, "[%d]", *dp);		/*     format a dimension   */
	(void) outstr(buf);                 /*     and print it out     */
      }
      if (xml) {
	sprintf(buf,"\"");
        outstr(buf);
      }
    }
    if (xml) {
      sprintf(buf,">");
      outstr(buf);
    }
    
}

void print_data(string tag, string type, int *dims)
{
    string fmt;
    size_t dlen;
    char buf[BUFLEN];  /* danger: buffer overflow */
    byte *dat, *dp;
    bool Qprint;

    Qprint = !(*testtag && streq(tag,testtag)==0);  /* if only print testtag... */

    fmt = find_fmt(type);
    dlen = get_dlen(instr, tag);
    dat = (byte*) allocate(dlen+1);
    get_data_sub(instr, tag, type, dat, dims, FALSE);
    if (streq(type, CharType) && Qprint)
	(void) outstr(dims != NULL ? "\"" : "\'");
    for (dp = dat; dp < dat + dlen; ) {		/* loop over data array     */
	if (streq(type, AnyType)) {		/*   output generic data?   */
	    sprintf(buf, fmt, *((byte *) dp));
	    dp += sizeof(byte);
	} else if (streq(type, CharType)) {	/*   output readable chars? */
	    sprintf(buf, fmt, *((char *) dp));
	    dp += sizeof(char);
	} else if (streq(type, ByteType)) {	/*   output bytes of data?  */
	    sprintf(buf, fmt, *((byte *) dp));
	    dp += sizeof(byte);
	} else if (streq(type, ShortType)) {	/*   output short integers? */
	    sprintf(buf, fmt, *((short *) dp));
	    dp += sizeof(short);
	} else if (streq(type, IntType)) {	/*   output standard ints?  */
	    sprintf(buf, fmt, *((int *) dp));
	    dp += sizeof(int);
	} else if (streq(type, LongType)) {	/*   output long integers?  */
	    sprintf(buf, fmt, *((long *) dp));
	    dp += sizeof(long);
	} else if (streq(type, HalfpType)) {	/*   output short int? */
	    sprintf(buf, fmt, *((short *) dp));
	    dp += sizeof(short);
	} else if (streq(type, FloatType)) {	/*   output floating point? */
	    sprintf(buf, fmt, *((float *) dp));
	    dp += sizeof(float);
	} else if (streq(type, DoubleType)) {	/*   output double numbers? */
	    sprintf(buf, fmt, *((double *) dp));
	    dp += sizeof(double);
	} else
	    error("print_data: type %s unknown\n", type);
	if (Qprint)
	  if (! outstr(buf))
	    break;
    }
    if (streq(type, CharType) && Qprint)
	(void) outstr(dims != NULL ? "\"" : "\'");
    if (xml && Qprint) {
      sprintf(buf," </%s>",tag);
      outstr(buf);
    }
    free((char *)dat);
}

/*
 * TYPEINFO: struct for info about data types.
 */

struct typeinfo {
    string titype;		/* type code used by filestruct package     */
    string tiname;		/* C-style declarator for this type         */
    string tilpfmt;		/* printf format for low-precision output   */
    string tihpfmt;		/* printf format for max-precision output   */
    string tixpfmt;             /* printf format for XML output */
};

/*
 * TYPETAB: table of info about data types.
 */

struct typeinfo typetab_o[] = {
    { AnyType,    "any",      "%#o ",     "%#o ",     "%#o", },
    { CharType,   "char",     "%c",       "%c",       "%c",  },
    { ByteType,   "byte",     "%#o ",     "%#o ",     "%#d ", },
    { ShortType,  "short",    "%#o ",     "%#o ",     "%#d ", },
    { IntType,    "int",      "%#o ",     "%#o ",     "%#d ", },
    { LongType,   "long",     "%#lo ",    "%#lo ",    "%#d ", },
    { HalfpType,  "halfp",    "%#o ",     "%#o ",     "%#d ", },
    { FloatType,  "float",    "%#g ",     "%15.8e ",  "%#g ", },
    { DoubleType, "double",   "%#g ",     "%23.16e ", "%#g ", },
    { SetType,    "set",      NULL,       NULL,        NULL, },
    { TesType,    "tes",      NULL,       NULL,        NULL, },
    { NULL, },
};


struct typeinfo typetab_d[] = {
    { AnyType,    "any",      "%#o ",     "%#o ",     "%#o", },
    { CharType,   "char",     "%c",       "%c",       "%c",  },
    { ByteType,   "byte",     "%#d ",     "%#d ",     "%#d ", },
    { ShortType,  "short",    "%#d ",     "%#d ",     "%#d ", },
    { IntType,    "int",      "%#d ",     "%#d ",     "%#d ", },
    { LongType,   "long",     "%#ld ",    "%#ld ",    "%#d ", },
    { HalfpType,  "halfp",    "%#d ",     "%#d ",     "%#d ", },
    { FloatType,  "float",    "%#g ",     "%15.8e ",  "%#g ", },
    { DoubleType, "double",   "%#g ",     "%23.16e ", "%#g ", },
    { SetType,    "set",      NULL,       NULL,        NULL, },
    { TesType,    "tes",      NULL,       NULL,        NULL, },
    { NULL, },
};

string find_name(string type)
{
    struct typeinfo *tip;

    for (tip = typetab_o; tip->titype != NULL; tip++)
	if (streq(tip->titype, type))
	    return (tip->tiname);
    error("find_name: type %s unknown\n", type);
    return NULL;   /* make stringent compilers happy */
}

string find_fmt(string type)
{
    struct typeinfo *tip;

    for (tip = octal ? typetab_o : typetab_d; tip->titype != NULL; tip++)
	if (streq(tip->titype, type))
	  return (xml ? tip->tixpfmt :
		  maxprec ? tip->tihpfmt : tip->tilpfmt);
    error("find_fmt: type %s unknown\n", type);
    return NULL;   /* make stringent compilers happy */
}

/*
 * State variables for formatting functions.
 */

int curcoll = 0;	/* current cursor column, numberd from 0 */
int curline = 0;	/* current line of item, numbered from 0 */

int baseindent = 0;	/* indentation of item */
int actindent;		/* actual indentation of line */

bool outstr(string str)
{
    if (curcoll + (int)strlen(str) >= margin) { /* string too long? */
	printf("\n");				/*   begin next line */
	curcoll = 0;
	curline++;
	if ((!allline) && (curline == maxline))	/*   too much output? */
	    str = ". . .";
	actindent = baseindent + indent;	/*   indent text of item */
    }
    while (curcoll < actindent) {		/* space over to start */
	printf(" ");
	curcoll++;
    }
    printf("%s", str);                          /* output string */
    curcoll = curcoll + strlen(str);            /* and count length */
    return (allline || (curline < maxline));
}

void end_line()
{
    printf("\n");
    curcoll = curline = 0;
    actindent = baseindent;
}

void change_indent(int delta)
{
    baseindent = baseindent + delta;
    actindent = baseindent;
}
