/*
 *	Example of a C program
 *
 *       On *Unix* to be compiled and linked as:
 *       
 *       cc -g -o mainc main.c
 *     
 *       The flags '-g' are defensive programming flags
 *       and are optional (debugging in this case)
 */

#include <stdio.h>		/* standard C include file */

#define NMAX  10

int main(int argc,char *argv[])
{
    float a;
    int   i, nmax;

    if (argc < 2)       /* if no arguments provided, take default */
        nmax = NMAX;
    else                /* else take it from the first argument */
        nmax = atoi(argv[1]);

    a = 1.0;
    for (i=0; i<nmax; i++) {
         a = a + a;
    }

    printf("The sum is %20.10f\n",a);
#ifdef WARNING
    warning("Trying a warning message");
#endif
    return 0;
}
