/*
 *  tabhead:    create a header for a table based on columns descriptions
 */

#include <nemo.h>
#include <extstring.h>

string defv[] = {
  "cols=\n        Column names",
  "mode=1\n       Mode:  0=cols names   1=col1,col2,...  2=A,B....AA,....",
  "prefix=#\n     Prefix before column names",
  "separ= \n      Separator between column names",
  "ncols=0\n      Number of columns in case cols= not given",
  "VERSION=0.2\n  3-apr-2023 PJT",
  NULL,
};

string usage="Generate an ascii table header";


local string az="ABCDEFGHIJKLMNOPQRSTUVWXYZ";

void nemo_main()
{
  string cols   = getparam("cols");
  string prefix = getparam("prefix");
  string separ  = getparam("separ");
  int mode  = getiparam("mode");
  int ncols = getiparam("ncols");
  int i, j, k;

  if (hasvalue("cols")) {
    printf("%s", prefix);
    string *colnames = burststring(cols, ", ");
    ncols = xstrlen(colnames,sizeof(string))-1;
    for (i=0; i<ncols; i++) {
      if (i>0) printf("%s ",separ);
      printf("%s", colnames[i]);
    }
    printf("\n");
    free(colnames);
  } else if (mode == 1) {    //   topcat mode (needs prefix="")
    printf("%s", prefix);
    for (i=0; i<ncols; i++) {
      if (i>0) printf("%s ",separ);
      printf("col%d",i+1);
    }
    printf("\n");
  } else if (mode == 2) {
    if (ncols > 702)  // 26*27
      error("spreadheet alphabet mode cannot support more than 702 columns");
    printf("%s", prefix);    
    for (i=0; i<ncols; i++) {
      if (i>0) printf("%s ",separ);
      j = i%26;
      k = i/26 - 1;
      //printf("%d  = %d * 26 + %d ", i,k,j);
      if (k<0)
	printf("%c", az[j]);
      else
	printf("%c%c", az[k],az[j]);
    }
    printf("\n");
  } else
    error("mode=%d not implemented", mode);
}
