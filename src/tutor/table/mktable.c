/*   MKTABLE:    create an interesting table with simple C math
 *
 *   compile and link with e.g.
 *
 *               cc -g -o mktable mktable.c -lm
 * 
 *      2-oct-2003     created                    PJT
 */

#include <stdio.h>
#include <math.h>

#define MAXCOL      16
#define MAXROW    1000
#define MAXLINELEN 256

int main(int argc, char *argv[])
{
  int i, j, nxcol, npt;
  int do_sqrt, do_log, do_exp;
  double x[MAXROW], y[MAXROW][MAXCOL];

  npt = 100;
  do_sqrt = 1;
  do_log  = 1;
  do_exp  = 1;

  for (i=0; i<npt; i++) {
    x[i] = (double)i / (double)npt;
    y[i][0] = sqrt(x[i]);
    y[i][1] = log(x[i]);
    y[i][2] = exp(x[i]);
  }



  for (i=0; i<npt; i++)
    printf("%g %g %g %g\n", x[i], y[i][0], y[i][1], y[i][2]);


}


