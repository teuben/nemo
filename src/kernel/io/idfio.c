/*
 * IDFIO : Input Directive File  I/O
 *
 *    ...many years in the making...
 *
 *	12-aug-2014	Created - static version 
 */

#include <stdinc.h>

#define IDF_COMMENT "c"
#define IDF_REAL    "r"
#define IDF_INT     "i"
#define IDF_BOOL    "b"
#define IDF_STRING  "s"
#define IDF_QSTRING "qs"

local string idf_types[] = { "c", "r", "i", "b", "s", "qs"};

typedef struct my_idf {
  int itype;     /* one if IDF_xxx types (1,2,3....) */
  string type;   /* same one, but the string value */
  string key;    /* keyword name (part of -- getargv) */
  string val;    /* string value of key */
  int row;       /* 1=1st row    */
  int col;       /* 1=1st column */
  string raw;    /* orginal line (raw entry) */
  string out;    /* new value to be output */
} IDF;

typedef struct my_line {
  int row;       /* 1=1st line */
  string data;   /* point to raw data */
} LINE;


IDF *idf_open_string(string *idfs);
IDF *idf_open_file(string file);
string *line_open_file(string file);
void idf_write(int nidf, IDF *idf, string file);


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
  "idf=\n         Input Directive File",
  "par=\n         Input parameter file",
  "out=\n         Output parameter file",
  "lineno=f\n     Add linenumbers to output?",
  "VERSION=1.0\n  12-aug-2014 PJT", 
  NULL,
};

string testidf[] = {
  "# a comment that is ignored by the parser",
  "c:comment   # a true comment line that will be skipped during parsing",
  "r:a         # the real parameter 'a'",
  "i:n         # an integer value",
  "b:q         # a boolean value",
  "s:s1        # a direct string",
  "qs:s2       # a quoted string",
  "i:n1 i:n2   # two integers",
  "c:comment   # another true parting comment",
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
  char *cp, line[MAX_LINELEN];
  bool Qline = getbparam("lineno");
  int argc, nidf, nidf2, nw;
  IDF *idf;
  string *argv, *av, *w;
  
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

  /* report, pre-parse and count the true IDF's */

  dprintf(1,"report IDF\n");
  nidf = 0;
  for (sp = idf0, lineno=0; *sp; sp++) {
    if (*sp[0] == '#') continue;
    lineno++;
    if (Qline)
      printf("%d: %s\n",lineno,*sp);	
    else
      printf("%s\n",*sp);	
    w = burststring(*sp," \t");
    nw = xstrlen(w,sizeof(string))-1;
    for (i=0; i<nw; i++) {
      if (*w[i] == '#') break;
      cp = strchr(w[i],':');
      if (cp==NULL) error("Missing : on %s",w[i]);
      nidf++;
    }
  } 
  dprintf(0,"Found %d idf in %d lines\n",nidf,lineno);

  /* fully parse IDF again */

  idf = (IDF *) allocate(nidf*sizeof(IDF));
  nidf2 = 0;
  for (sp = idf0; *sp; sp++) {
    if (*sp[0] == '#') continue;
    w = burststring(*sp," \t");
    nw = xstrlen(w,sizeof(string))-1;
    for (i=0; i<nw; i++) {
      if (*w[i] == '#') break;
      cp = strchr(w[i],':');
      if (cp==NULL) error("Missing : on %s",w[i]);
      *cp++ = 0;
      /* now w[i] points to type; cp to keyword */
      idf[nidf2].type = strdup(w[i]);
      idf[nidf2].key  = strdup(cp);
      idf[nidf2].row  = nidf2+1;
      idf[nidf2].col  = i+1;
      nidf2++;
    }
  }
  if (nidf2 != nidf) error("nidf2=%d != nidf=%d\n",nidf2,nidf);


  /* report the full IDF */

  for (i=0; i<nidf; i++) {
    printf("###: [%d,%d] %s %s\n", idf[i].row, idf[i].col, idf[i].type, idf[i].key);
  }

  if (hasvalue("par")) {
    pars = line_open_file(getparam("par"));
    n2 = xstrlen(pars,sizeof(string))-1;
    printf("Found %d lines in par file\n",n2);
    if (n2 != lineno) warning("par file not same as idf");

    for (l=0, i=0; l<lineno; l++) {    /* loop over all lines :   l counts lines, i counts idf's */
      /*  idf[i] is the current IDF */
      /*  pars[l] is the current line */
      /*  idf[i].row should match the current line */
      if (streq(idf[i].type,IDF_COMMENT)) {
	  idf[i].out = strdup(pars[l]);
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
	idf[i].out = strdup(w[iw]);
	i++;
      }
#else
      /* double checking */
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

    for (i=0; i<nidf; i++) {
      printf("###: [%d,%d] %-3s %-10s = %s\n", idf[i].row, idf[i].col, idf[i].type, idf[i].key, idf[i].out);
    }
  }

  /* report the command line tail after the -- */

  dprintf(1,"getargv\n");
  argv = getargv(&argc);
  dprintf(1,"  argc=%d\n",argc);
  if (argc>0) {
    argv++; /* skip the -- */
    for (av = argv; *av; av++) {
      printf("arg: %s\n",*av);
      cp = strchr(*av,'=');
      if (cp) {
	*cp++ = 0;
	for (i=0; i<nidf; i++) {
	  if (streq(*av,idf[i].key)) {
	    dprintf(0,"Patching key=%s with val=%s\n",*av,cp);
	    idf[i].out = strdup(cp);
	    break;
	  }
	}
	if (i==nidf) error("%d: key=%s not registered in IDF",i,*av);
      } else {
	error("cannot handle arugments yet that are not key=val");
      }
    }
  }


  for (i=0; i<nidf; i++) {
    printf(">>>: [%d,%d] %-3s %-10s = %s\n", idf[i].row, idf[i].col, idf[i].type, idf[i].key, idf[i].out);
  }
  


  if (hasvalue("out")) {
    istr = stropen(getparam("out"),"w");
    for (i=0; i<nidf; i++) {
      if (i>0 && idf[i].row != idf[i-1].row) 
	fprintf(istr,"\n");
      fprintf(istr,"%s ",idf[i].out);
    }
    fprintf(istr,"\n");
    strclose(istr);
  }


}


#endif
