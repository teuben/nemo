/*
 * IDFIO : Input Directive File  I/O
 *
 *    ...many years in the making...
 *
 *	12-aug-2014	Created - static version only, no library yet
 *      13-aug-2014     static, but with one open ended (last) array
 *
 *  WARNING:  this program has memory leaks like a bucket shot with a bbgun
 */

#include <stdinc.h>

#define IDF_COMMENT "c"
#define IDF_REAL    "r"
#define IDF_INT     "i"
#define IDF_BOOL    "b"
#define IDF_STRING  "s"
#define IDF_QSTRING "qs"

local string idf_types[] = { "c", "r", "i", "b", "s", "qs"};

  /* debugging output ? */
local bool Q = FALSE;

typedef struct cstring {
  string          val;
  struct cstring *nxt;
} cstring;

/*  'type:key'  or 'type:key[]'   or   'type:key[len]' (not yet) */
typedef struct my_idf {
  int itype;     /* one if IDF_xxx types (1,2,3....) */
  string type;   /* same one, but the string value */
  string key;    /* keyword name (part of -- getargv) */
  string len;    /* reference to a length (optional) */
  string val;    /* string value of key */
  int row;       /* 1=1st row    */
  int col;       /* 1=1st column */
  int nvals;     /* -1=scalar  0=to end of file >1: taken */
  string raw;    /* orginal line (raw entry) */
  string out;    /* new value to be output */
  cstring cout;  /* chained string of outvalues until 0 */
} IDF;

typedef struct my_line {
  int row;       /* 1=1st line */
  string data;   /* point to raw data */
} LINE;


IDF    *idf_open_string(string *idfs);
IDF    *idf_open_file(string file);
void    idf_write(int nidf, IDF *idf, string file);
string *line_open_file(string file);


void strip_newline(char *line) {
  int n;
  if (line == 0) return;
  if (*line == 0) return;
  n = strlen(line);
  if (line[n-1] == '\n') line[n-1] = 0;
}


string *line_open_file(string fname) {
  stream istr;
  char *cp, line[MAX_LINELEN];
  int i, n;
  string *lines;

  n = nemo_file_lines(fname,0);
  lines = (string *) allocate((n+1)*sizeof(string));

  istr = stropen(fname,"r");
  for (i=0; i<n; i++) {
    if (fgets(line,MAX_LINELEN,istr) == NULL) error("short file");
    strip_newline(line);
    lines[i] = strdup(line);
  }
  lines[n] = 0;
  strclose(istr);
  return lines;
}


#ifdef TOOLBOX


#include <getparam.h>

string defv[] = {
  "idf=???\n      Input Directive File (use blank for self-generated)",
  "par=\n         Input parameter file",
  "out=\n         Output parameter file",
  "lineno=f\n     Add linenumbers to output?",
  "checktype=f\n  Type checking on parameters?",
  "report=f\n     Report all final RDF key=val pairs",
  "VERSION=1.1a\n 14-aug-2014 PJT", 
  NULL,
};

string testidf[] = {
  "# a comment that is ignored by the parser",
  "c:comment1  # a true comment line that will be skipped during parsing",
  "r:a         # the real parameter 'a'",
  "i:n         # an integer value",
  "b:q         # a boolean value",
  "s:s1        # a direct string",
  "qs:s2       # a quoted string",
  "i:n1 i:n2   # two integers",
  "c:comment2  # another true parting comment",
  "# a parting comment",
  NULL,
};

string usage = "IDF I/O routine tester";

nemo_main()
{
  stream istr;
  float fval;
  double dval;
  int i, k, l, n, n2, iw, i1, i2, lineno, ival;
  string *sp, *idf0, *pars;
  char *cp, *cp1, line[MAX_LINELEN];
  bool Qline = getbparam("lineno");
  bool Qshow_idf;
  bool Qtype = getbparam("checktype");
  bool Qreport = getbparam("report");
  int argc, nidf, nidf2, nw, nopen, nrow;
  IDF *idf;
  string *argv, *av, *w;
  cstring *csp, *csn;
  
  if (hasvalue("idf")) {
    n = nemo_file_lines(getparam("idf"),0);
    dprintf(1,"%s : nlines=%d\n",getparam("idf"),n);
    idf0 = (string *) allocate((n+1)*sizeof(string));
    istr = stropen(getparam("idf"),"r");
    for (i=0; i<n; i++) {
      if (fgets(line,MAX_LINELEN,istr) == NULL) error("short file");
      strip_newline(line);
      idf0[i] = strdup(line);
    }
    idf0[n] = 0;
    strclose(istr);
  } else {
    warning("Testing IDF");
    idf0 = testidf;
  }

  /* report, and pre-parse and count the true IDF's */

  dprintf(1,"report IDF\n");
  Qshow_idf = !hasvalue("out");
  nidf = 0;
  for (sp = idf0, lineno=0; *sp; sp++) {
    if (*sp[0] == '#') continue;
    lineno++;
    if (Qshow_idf) {
      if (Qline)
	printf("%d: %s\n",lineno,*sp);	
      else
	printf("%s\n",*sp);	
    }
    w = burststring(*sp," \t");
    nw = xstrlen(w,sizeof(string))-1;
    for (i=0; i<nw; i++) {
      if (*w[i] == '#') break;
      cp = strchr(w[i],':');
      if (cp==NULL) error("Missing : on %s",w[i]);
      nidf++;
    }
  } 
  dprintf(0,"Found %d IDF_parameters in %d lines in idf file\n",nidf,lineno);

  /* now fully parse IDF */

  idf = (IDF *) allocate(nidf*sizeof(IDF));
  nidf2 = 0;
  nopen = 0;
  nrow = 0;
  for (sp = idf0; *sp; sp++) {
    if (*sp[0] == '#') continue;
    nrow++;
    if (nopen) error("Cannot handle any parameters after an open ended array");
    w = burststring(*sp," \t");
    nw = xstrlen(w,sizeof(string))-1;
    for (i=0; i<nw; i++) {
      if (*w[i] == '#') break;
      cp = strchr(w[i],':');
      if (cp==NULL) error("Missing : on %s",w[i]);
      *cp++ = 0;
      /* now w[i] points to type; cp to keyword, possibly with a [] */
      idf[nidf2].type = strdup(w[i]);
      idf[nidf2].key  = strdup(cp);
      idf[nidf2].row  = nrow;
      idf[nidf2].col  = i+1;
      cp1 = strchr(cp,'[');
      if (cp1) {
	*cp1++ = 0;
	idf[nidf2].key  = strdup(cp);
	if (*cp1 == ']') {
	  if (nopen) error("Cannot handle more than one open ended array");
	  idf[nidf2].nvals = 0;
	  nopen++;
	} else
	  error("fixed dimensioned arrays not allowed yet : %s",cp);
      } else {
	idf[nidf2].key   = strdup(cp);
	idf[nidf2].nvals = -1;
      }
      nidf2++;
    }
  }
  if (nidf2 != nidf) error("nidf2=%d != nidf=%d\n",nidf2,nidf);


  /* report the full IDF */

  if (Q) {
      for (i=0; i<nidf; i++) {
      dprintf(1,"###: [%d,%d] %s %s\n", idf[i].row, idf[i].col, idf[i].type, idf[i].key);
    }
  }

  if (hasvalue("par")) {
    pars = line_open_file(getparam("par"));
    n2 = xstrlen(pars,sizeof(string))-1;
    dprintf(0,"Found %d lines in par file\n",n2);
    if (Q)
      if (n2 != lineno) warning("par file not same as idf");

    for (l=0, i=0; l<n2; l++) {    /* loop over all lines :   l counts lines, i counts idf's */
      /*  idf[i] is the current IDF */
      /*  pars[l] is the current line */
      /*  idf[i].row should match the current line */


      /* a comment spans (by definition) the whole line */
      if (i<nidf && streq(idf[i].type,IDF_COMMENT)) {
	  idf[i].out = strdup(pars[l]);
	  idf[i].cout.val = strdup(pars[l]);
	  idf[i].cout.nxt = 0;
	  i++;
	  continue;
      } 
      /* for now any types are not comments, and so treated same way */
      w = burststring(pars[l]," \t");
      nw = xstrlen(w,sizeof(string))-1;

#if 1
      /* just do it, no checking */
      for (iw=0; iw<nw; iw++) {
	if (*w[iw] == '#') break;
	if (i < nidf) {
	  idf[i].out = strdup(w[iw]);
	  idf[i].cout.val = strdup(w[iw]);
	  idf[i].cout.nxt = 0;
	} else {
	  if (nopen==0) error("Too many values, and no open ended array");

	  csn = (cstring *) allocate(sizeof(cstring));
	  csn->val = strdup(w[iw]);
	  csn->nxt = 0;

	  csp = &idf[nidf-1].cout;
	  while (csp->nxt)
	    csp = csp->nxt;
	  csp->nxt = csn;
	}
	i++;
      }
      if (Q) {
	if (i==nidf) {
	  warning("end of idf %d %d",l,n2);
	  if (l<n2-1) warning("not exhausting lines in par file");
	}
      }
#else
      /* messy double checking */
      for (i1=i; i<nidf && idf[i1].row==idf[i].row; i1++) {
	dprintf(0,"i=%d i1=%d row=%d row1=%d\n",i,i1,idf[i].row,idf[i1].row);
      }
      n2 = i1-i;
      if (nw != n2) error("nw=%d != n2=%d i=%d l=%d (%s)",nw,n2,i,l,pars[l]);
      for (i1=i, iw=0; i<nidf && idf[i1].row==idf[i].row; i1++, iw++) {
	idf[i1].out = strdup(w[iw]);
      }
#endif
    }

    if (Q) {
      for (i=0; i<nidf; i++) {
	printf("###: [%d,%d] %-3s %-10s = %s\n", idf[i].row, idf[i].col, idf[i].type, idf[i].key, idf[i].out);
      }
    }
  }

  /* report the command line tail after the -- */

  dprintf(1,"getargv\n");
  argv = getargv(&argc);
  dprintf(1,"  argc=%d\n",argc);
  if (argc>0) {
    argv++; /* skip the -- */
    if (Qtype) warning("Type checking not implemented yet");
    for (av = argv; *av; av++) {
      if (Q) printf("arg: %s\n",*av);
      cp = strchr(*av,'=');
      if (cp) {
	*cp++ = 0;
	for (i=0; i<nidf; i++) {
	  if (streq(*av,idf[i].key)) {
	    dprintf(1,"Patching key=%s with val=%s\n",*av,cp);
	    if (idf[i].nvals < 0) {
	      idf[i].out      = strdup(cp);
	      idf[i].cout.val = strdup(cp);
	      idf[i].cout.nxt = 0;
	    } else {
	      w = burststring(cp," ,");
	      nw = xstrlen(w,sizeof(string))-1;
	      idf[i].out      = strdup(w[0]);
	      idf[i].cout.val = strdup(cp);
	      idf[i].cout.nxt = 0;
	      csn = &idf[i].cout;
	      for (iw=0; iw<nw; iw++) {
		csn->val = strdup(w[iw]);
		if (iw<nw-1) {
  		  csn->nxt = (cstring *) allocate(sizeof(cstring));
		  csn = csn->nxt;
		} else
		  csn->nxt = 0;
	      }
	    }
	    break;
	  }
	}
	if (i==nidf) error("%d: key=%s not registered in IDF",i,*av);
      } else {
	error("cannot handle arugments yet that are not key=val");
      }
    }
  }

  if (Q) {
    for (i=0; i<nidf; i++) {
      printf(">>>: [%d,%d] %-3s %-10s = %s\n", idf[i].row, idf[i].col, idf[i].type, idf[i].key, idf[i].out);
    }
  }
  
  if (hasvalue("out")) {
    istr = stropen(getparam("out"),"w");
    for (i=0; i<nidf; i++) {
      if (i>0 && idf[i].row != idf[i-1].row) 
	fprintf(istr,"\n");
      if (idf[i].col > 1) fprintf(istr," ");
      if (streq(idf[i].type,"qs") && idf[i].out[0] != '\'')
	fprintf(istr,"'%s'",idf[i].out);
      else
	fprintf(istr,"%s",idf[i].out);

    }
    fprintf(istr,"\n");
    if (nopen) {
      csp = idf[nidf-1].cout.nxt;
      while (csp) {
	fprintf(istr,"%s\n",csp->val);	
	csp = csp->nxt;
      }
    }
    strclose(istr);
  }

  if (Qreport) {
    for (i=0; i<nidf; i++) {
      printf("%s=%s\n",idf[i].key,idf[i].out);
    }
  }

}


#endif
