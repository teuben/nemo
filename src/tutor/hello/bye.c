/*
 *	a very simple (ANSI) C program, no NEMO interface
 *      ANSI C requires main() to be of type int !!!
 */

#include <stdio.h>

int main(int argc, char *argv[])
{
    int i;

    for (i=0; i<argc; i++)      /* report all command line arguments */
        printf("%s ",argv[i]);
    printf("\nBye.\n");         /* bye */
    return 0;                   /* note return is NOT a function */
}
