/*   MKTABLE:    create an 'interesting' table 
 *   
 *   See also: $NEMO/src/tutor/table/mktable.c
 *
 *   compile and link with e.g.
 *
 *               cc -g -o mktable mktable.c -lm
 * 
 *      7-oct-2003     created                  PJT
 */

#include <stdio.h>
#include <math.h>

#define MAXCOL       6
#define MAXROW    1000

/*
 * normalize all elements from an array to the array length
 */

void normalize(double *a, int n)
{
  int i;

  for (i=0; i<n; i++)
    a[i] /= n;            /* same as:    a[i] = a[i]/n   */
}



/*
 *  start of main program
 */

int main(int argc, char *argv[])
{
  int i, j, npt;
  double x[MAXROW], y[MAXCOL][MAXROW], sum0, sum1, sum2;

  npt = 100;         /* set the number of points */

  if (npt >= MAXROW) {      /* make sure program has enough array space */
    printf("npt=%d too large -or- MAXROW=%d too small\n",npt,MAXROW);
    exit(1);
  }

  sum0 = sum1 = sum2 = 0.0;         /* set to 0 */

  for (i=0; i<npt; i++) {               /* compute and accumulate */
    x[i] = (double)(i+1) / (double)npt;
    y[0][i] = sqrt(x[i]);
    y[1][i] = log(x[i]);
    y[2][i] = exp(x[i]);
    sum0 += y[0][i];
    sum1 += y[1][i];
    sum2 += y[2][i];
    y[3][i] = sum0;
    y[4][i] = sum1;
    y[5][i] = sum2;

  }

  normalize(y[3],npt);             /* normalize the cumulatives */
  normalize(y[4],npt);
  normalize(y[5],npt);

  for (i=0; i<npt; i++) {                              /* print out */
    printf("%g", x[i]);                                        /* x */
    printf(" %g %g %g", y[0][i], y[1][i], y[2][i]);    /* functions */
    printf(" %g %g %g", y[3][i], y[4][i], y[5][i]);  /* cumulatives */
    printf("\n");
  }

  return 0;     /* must return an 'int' */

}


