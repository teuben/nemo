/*
 *	Example of a NEMO C program
 *
 *       On *Unix* to be compiled and linked as:
 *       
 *       cc -g -o main main.c $NEMOLIB/libnemo.a -lm
 * or:
 *	 cc -g -o main main.c -lnemo -lm
 *
 * within the NEMO environment.
 *     
 *
 *  22-oct-90	V1.0	Created		    PJT
 *  21-may-92	V1.1	added usage	    PJT
 *  12-feb-95   V1.2    added format=       PJT
 */

#include <stdinc.h>                             /* standard NEMO include */
#include <getparam.h>

string defv[] = {                               /* keywords definitions */
    "nmax=10\n          Number of iterations",
    "format=%20.10f\n	Format to print result",
    "VERSION=1.2\n      12-feb-95 PJT",
    NULL,
};

string usage = "Example C program with nemo_main convention";	/* Usage */

void nemo_main(void)                            /* NEMO's program main entry */
{
    real a;         /* becomes 'float' or 'double' depending on compile flag */
    int   i, nmax;
    char fmt[64];

    nmax = getiparam("nmax");      /* obtain 'nmax' from commmand line */
    if (nmax<1) warning("%d: Unexpected value for nmax",nmax);
    dprintf(1,"Iteration counter = %d\n",nmax);
    
    a = 1.0;
    for (i=0; i<nmax; i++) {              /* loop the loop */
         a = a + a;
    }

    sprintf(fmt,"The sum is %s\n",getparam("format")); /* set format string */
    printf(fmt,a);                                         /* output */
}
