/*
 *  TABTAB:  combine tables, using concat or paste
 */

#include <nemo.h>
#include <extstring.h>
#include <table.h>

string defv[] = {
  "in=\n            Input tables",
  "out=\n           Output table",
  "mode=paste\n     Cat (by row) or Paste (by column)",
  "VERSION=0.1\n    26-apr-2022 PJT",
  NULL,
};

string usage="concatenate or paste conformant tables";


void nemo_main()
{
  string out = getparam("out");
  string mode = getparam("mode");

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

  // table_cat


  // output
  
}

