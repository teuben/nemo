/*
 *	a very simple (ANSI) C program, no NEMO interface
 */

#include <stdio.h>

void main(int argc, char *argv[])
{
    int i;

    for (i=0; i<argc; i++)
        printf("%s ",argv[i]);
    printf("\nBye.\n");
}
