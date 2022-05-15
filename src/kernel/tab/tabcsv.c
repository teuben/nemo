/*
 *  TABCSV:  convert table to CSV (or any symbol-separated-value table)
 */

#include <nemo.h>
#include <extstring.h>
#include <table.h>

string defv[] = {
  "in=???\n         Input tables",
  "mode=csv\n       separator mode (also: LF, HT,)",
  "addrow=f\n       add a row counter?",
  "out=-\n          Output table",
  "VERSION=0.1\n    9-may-2022 PJT",
  NULL,
};

string usage="convert table to CSV";


void nemo_main()
{
  string mode = getparam("mode");
  string infile = getparam("in");
  string *cols;
  char colsep = mode[0];
  bool Qrow = getbparam("addrow");

  // LF  \n
  // HT  \t
  if (streq(mode,"csv"))
    colsep = ',';
  else if (streq(mode,"LF"))
    colsep = '\n';
  else if (streq(mode,"HT"))
    colsep = '\t';
  dprintf(1,"table mode=%s\n",mode);

  tableptr t = table_open(stropen(infile,"r"), 0);
  int nrows = table_nrows(t);
  int ncols = table_ncols(t);
  
  for (int i=0; i<nrows; i++) {
    cols = table_rowsp(t,i);
    if (Qrow) printf("%d%c",i+1,colsep);
    printf("%s",cols[0]);
    for (int j=1; j<ncols; j++)
      printf("%c%s",colsep,cols[j]);
    printf("\n");
  }
}

