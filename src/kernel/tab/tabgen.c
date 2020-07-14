/*
 *  TABGEN:     table with random numbers
 */

#include <nemo.h>

string defv[] = {
    "out=???\n     Output table",
    "nr=10\n       Number of rows",
    "nc=5\n        Number of columns",
    "mode=0\n      Mode (0=uniform, 1=normal)",
    "seed=123\n    Random seed",
    "fmt=%g\n      Format statement for output",
    "VERSION=0.3\n 14-Jul-2020 XYZ",
    NULL,
};

string usage="Create a table with random numbers";

string cvsid="$Id:$";

typedef real (*my_real_proc)(real,real);

void nemo_main()
{
  stream ostr = stropen(getparam("out"),"w");
  int i, nr = getiparam("nr");
  int j, nc = getiparam("nc");
  int mode = getiparam("mode");
  int seed = init_xrandom(getparam("seed"));
  char fmt[64];
  real x;
  my_real_proc my_random;

  dprintf(1,"seed=%d\n",seed);

  sprintf(fmt,"%s ",getparam("fmt"));

  if (mode <= 0)
    my_random = xrandom;
  else
    my_random = grandom;

  if (mode < 0) {
    // only produce random numbers, testing production rate
    real *x = (real *) allocate(nc*nr*sizeof(real));
    int k=0;
    for (i=0; i<nr; i++)
      for (j=0; j<nc; j++)
	x[k++] = xrandom(0.0,1.0);
  } else
    // output to file as well
    for (i=0; i<nr; i++) {
      for (j=0; j<nc; j++)
	fprintf(ostr,fmt, my_random(0.0, 1.0));
      fprintf(ostr,"\n");
    }
}
