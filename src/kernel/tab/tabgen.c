/*
 *  TABGEN:     table with random numbers
 */

#include <nemo.h>

string defv[] = {
    "out=???\n     Output table",
    "nr=10\n       Number of rows",
    "nc=5\n        Number of columns",
    "mode=0\n      Mode",
    "seed=123\n    Random seed",
    "VERSION=0.2\n 13-Jul-2020 XYZ",
    NULL,
};

string usage="Create a table with random numbers";

string cvsid="$Id:$";

void nemo_main()
{
  stream ostr = stropen(getparam("out"),"w");
  int i, nr = getiparam("nr");
  int j, nc = getiparam("nc");
  int mode = getiparam("mode");
  int seed = init_xrandom(getparam("seed"));

  dprintf(1,"seed=%d\n",seed);

  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++)
      fprintf(ostr,"%g ", xrandom(0.0, 1.0));
    fprintf(ostr,"\n");
  }

}
