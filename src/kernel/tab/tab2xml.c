/*
 * TAB2XML: table into XML, mostly to produce a VOTable
 *
 *      13-feb-04   created
 *      15-feb-04   V1.1  - decided to keep the traditional out=-
 *
 * TODO:
 *    - "TABLE name=" is wrong, it's the input file
 */

#include <stdinc.h>
#include <getparam.h>
#include <table.h>
#include <extstring.h>
#include <ctype.h>

string defv[] = {
  "in=???\n           input ascii table",
  "out=???\n          output VOTable",
  "col=\n             column numbers to select",
  "colnames=\n        Column names",
  "VERSION=1.1\n      15-feb-04 PJT",
  NULL
};

string usage = "Select columns from a table and convert to VOTable/XML";


string input, output;			/* file names */
stream instr, outstr;			/* file streams */
int   ninput;				/* number of input files */

#ifndef MAX_COL
#define MAX_COL 256
#endif

#ifndef MAX_LINELEN
#define MAX_LINELEN  2048
#endif


#define MNEWDAT          80		/* space needed for one number */

int    keep[MAX_COL+1];                 /* column numbers to keep */
int    nkeep;                           /* actual number of skip columns */
int    maxcol = 0;                      /* largest column number */
char   colsep;                          /* character to separate columns */
string *colnames;
bool   Qall;

local void setparams(void);
local void convert(stream , stream);
local void tab2space(char *);

extern  string *burststring(string, string);

local string vot_header=
"<?xml version='1.0'?>\n"
"<!DOCTYPE VOTABLE SYSTEM 'http://us-vo.org/xml/VOTable.dtd'>\n"
"<VOTABLE version='1.0'>\n"
"<!--\n"
" !  VOTable written by NEMO's tab2xml program\n"
" !-->\n"
"<RESOURCE>";

local string vot_trailer=
"</TABLEDATA>\n"
"</DATA>\n"
"</TABLE>\n"
"</RESOURCE>\n"
"</VOTABLE>";


void nemo_main()
{
  setparams();
  convert (instr,outstr);
}

local void setparams(void)
{
  int i;
  
  input = getparam("in");
  output = getparam("out");
  instr  = stropen(input,"r");
  outstr  = stropen(output,"w");
  
  if (hasvalue("col")) {
    Qall = FALSE;
    nkeep = nemoinpi(getparam("col"),keep,MAX_COL);
    if (nkeep<=0 || nkeep>MAX_COL)
      error("Too many columns given (%d) to keep (MAX_COL %d)",
	    nkeep,MAX_COL);
    for (i=0; i<nkeep; i++){
      if (keep[i]<=0) error("Illegal column number %d",keep[i]);
      maxcol = MAX(maxcol,keep[i]);
    }
  } else {
    Qall = TRUE;
    nkeep = -1;   /* trigger that nkeep needs to be set */
    for (i=0; i<=MAX_COL; i++) keep[i] = i;
  }
  if (hasvalue("colnames")) {
    colnames = burststring(getparam("colnames"),",");
  } else
    colnames = NULL;
}


local void convert(stream instr, stream outstr)
{
  char   line[MAX_LINELEN];          /* input linelength */
  int    i, nlines, noutv;
  string *outv;                   /* pointer to vector of strings to write */
  char   *seps=", \t";       /* column separators  */
  static bool first = TRUE;
  
  nlines=0;               /* count lines read so far */
  
  fprintf(outstr,"%s\n",vot_header);
  fprintf(outstr,"<TABLE name=\"%s\">\n",input);
 
  for (;;) {
    if (get_line(instr, line) < 0)      /* EOF */
      break;
    
    dprintf(3,"LINE: (%s)\n",line);
    if (iscomment(line)) continue;
    
    nlines++;
    tab2space(line);                  /* work around a Gipsy (?) problem */
    
    outv = burststring(line,seps);
    noutv = xstrlen(outv,sizeof(string)) - 1;
    if (noutv < maxcol)
      error("Too few columns in input file (%d < %d)",noutv,maxcol);
    if (Qall) { 
      if (first) {
	for (i=0; i<noutv; i++)
	  if (colnames)
	    fprintf(outstr,"<FIELD datatype=\"double\" name=\"%s\"/>\n",colnames[i]);
	  else
	    fprintf(outstr,"<FIELD datatype=\"double\" name=\"col%d\"/>\n",i+1);
	fprintf(outstr,"<DATA>\n");
	fprintf(outstr,"<TABLEDATA>\n");
	first = FALSE;
      }
      nkeep = noutv;
    } else if (first) {
      for (i=0; i<nkeep; i++)
	  if (colnames)
	    fprintf(outstr,"<FIELD datatype=\"double\" name=\"%s\"/>\n",colnames[i]);
          else
	    fprintf(outstr,"<FIELD datatype=\"double\" name=\"col%d\"/>\n",i+1);
      fprintf(outstr,"<DATA>\n");
      fprintf(outstr,"<TABLEDATA>\n");
      first = FALSE;
    }
    fprintf(outstr,"  <TR>\n");	
    for (i=0; i<nkeep; i++)
      fprintf(outstr,"    <TD>%s</TD>\n",outv[keep[i]]);
    fprintf(outstr,"  </TR>\n");	
  }
  fprintf(outstr,"%s\n",vot_trailer);
}
/*
 * small helper function, replaces tabs by spaces before processing.
 * this prevents me from diving into gipsy parsing routines and  fix
 * the problem there 
 * PJT - June 1998.
 */

local void tab2space(char *cp)
{
    while (*cp) {
        if (*cp == '\t') *cp = ' ';
        cp++;
    }
}

