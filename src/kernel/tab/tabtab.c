/*
 *  TABTAB:  combine tables, using concat or paste
 */

#include <nemo.h>
#include <extstring.h>
#include <table.h>

string defv[] = {
  "in=\n            Input tables",
  "out=\n           Output table",
  "mode=paste\n     cat (by row) or paste (by column)",
  "space=t\n        extra separator?",
  "VERSION=0.3\n    2-jun-2022 PJT",
  NULL,
};

string usage="concatenate or paste conformant tables";


void nemo_main()
{
  string mode = getparam("mode");
  bool Qspace = getbparam("space");
  char sep[8];
  if (Qspace)
    sprintf(sep," ");
  else
    sep[0] = '\0';
  warning("SEP:'%s'", sep);

  string *infiles = burststring(getparam("in")," ,");
  int nfiles = xstrlen(infiles,sizeof(string))-1;

  dprintf(0,"Processing %d tables\n",nfiles);
  if (nfiles < 2) error("need > 1 input file");

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

#if 1
  // if the table has no header, we can now do it manually use the space separator
  // but this case can be done with the unix "cat" and "paste" commands as well
  stream o = stropen(getparam("out"),"w");
  
  if (*mode == 'c') {
    for (int i=0; i<nfiles; i++)
      for (int j=0; j<table_nrows(t[i]); j++)
	fprintf(o,"%s\n",table_row(t[i],j));
  } else if (*mode == 'p') {
    for (int j=0; j<table_nrows(t[0]); j++) {
      fprintf(o,"%s ",table_row(t[0],j));      
      for (int i=1; i<nfiles; i++)
	fprintf(o,"%s%s",table_row(t[i],j),sep);
      fprintf(o,"\n");
    }
  }
  
  strclose(o);
#endif

  // table_cat

  
  // output
  
}

