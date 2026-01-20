/*
 *  TABTAB:  combine tables, using concat or paste
 *           add line numbers if a single file is input
 */

#include <nemo.h>
#include <extstring.h>
#include <table.h>

string defv[] = {
  "in=???\n         Input tables",
  "out=-\n          Output table",
  "mode=paste\n     cat (by row) or paste (by column)",
  "space=t\n        extra separator?",
  "offset=0\n       line number offset",
  "VERSION=0.4\n    19-jan-2026 PJT",
  NULL,
};

string usage="concatenate or paste conformant tables, or add line numbers";


void nemo_main()
{
  string mode = getparam("mode");
  bool Qspace = getbparam("space");
  int offset = getiparam("offset");
  char sep[8];
  if (Qspace)
    sprintf(sep," ");
  else
    sep[0] = '\0';
  //warning("SEP:'%s'", sep);

  string *infiles = burststring(getparam("in")," ,");
  int nfiles = xstrlen(infiles,sizeof(string))-1;

  dprintf(1,"Processing %d tables\n",nfiles);

  // allocate and open all tables

  tableptr *t = (tableptr *) allocate(nfiles*sizeof(tableptr));
  for (int i=0; i<nfiles; i++)
    t[i] = table_open(stropen(infiles[i],"r"), 0);

  //  make sure tables are conformant, depending on mode
   
  if (*mode == 'c') {          // cat, they all need the same #cols
    for (int i=1; i<nfiles; i++)
      if (table_ncols(t[i]) != table_ncols(t[0]))
	  error("%s not same # cols",infiles[i]);
  } else if (*mode == 'p') {   // paste, they all need the same #rows
    for (int i=1; i<nfiles; i++)
      if (table_nrows(t[i]) != table_nrows(t[0]))
	error("%s not same # rows",infiles[i]);
  } else
    error("unknown mode %s",mode);

  // if the table has no header, we can now do it manually use the space separator
  // but this case can be done with the unix "cat" and "paste" commands as well
  stream o = stropen(getparam("out"),"w");
  
  if (*mode == 'c') {
    for (int i=0; i<nfiles; i++)
      for (int j=0; j<table_nrows(t[i]); j++)
	fprintf(o,"%s\n",table_row(t[i],j)); 
  } else if (*mode == 'p') {
    for (int j=0; j<table_nrows(t[0]); j++) {
      if (nfiles==1) fprintf(o,"%d ",j+offset);
      fprintf(o,"%s ",table_row(t[0],j));      
      for (int i=1; i<nfiles; i++)
	fprintf(o,"%s%s",table_row(t[i],j),sep);
      fprintf(o,"\n");
    }
  }
  
  strclose(o);

  // table_cat
  // output
  
}

