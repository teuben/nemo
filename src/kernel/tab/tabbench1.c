/*
 *     Benchmark some common table I/O operations
 *
 *      3-jul-2020  V0.1    drafted
 */

//1    nbody=1000000      (for 10M there is a bug)
//2    mkplummer - $nbody | snapprint - format=%20.16f    > p1M.tab
//3    /usr/bin/time tabbench1 p1M.tab .
//4    /usr/bin/time tabtranspose p1M.tab p1Mt.tab $nbody
//5    /usr/bin/time tabbench1 p1Mt.tab .


#include <stdlib.h>
#include <ctype.h>
#include <stdinc.h>
#include <getparam.h>
#include "table.h"

string defv[] = {
    "in=???\n	       input file",
    "out=???\n         output file",
    "nmax=10000\n      Default max allocation",
    "mode=1\n          Benchmark mode",
    "VERSION=0.2\n     24-jul-2020 PJT",
    NULL,
};

string usage="table I/O benchmark";
int     nmax;                           /* # lines in file */
int     kmin;                           /* # columns to transpsse */

#ifndef MAX_LINELEN 
#define MAX_LINELEN  2048
#endif

void nemo_main(void)
{
    stream istr, ostr;
    char *line;
    int *select = NULL;
    int nout, next = 0, counter = 0;
    int    i, j;
    string input = getparam("in");
    string output = getparam("out");
    size_t buffer_size = MAX_LINELEN, bufflen;

    line = malloc((MAX_LINELEN) * sizeof(char));
    nmax = nemo_file_lines();

    dprintf(0,"MAX_LINELEN=%d\n",MAX_LINELEN);
    dprintf(0, "Input File: %s\n", input);
    dprintf(0, "Output File: %s\n", output);
    dprintf(0, "Name: %s\n", usage);

    istr = stropen(getparam("in"),"r");
    ostr = stropen(getparam("out"),"w");

    i = 0;
    
    while (getline(&line, &(buffer_size), istr) != -1) {
        counter = 0;
        while(isspace(line[counter]) != 0) {
            ++counter;
        }

        if (line[counter] != '#' && line[counter] != '/' && line[counter] != '!') {
            i++;
	        fputs(line,ostr);
        }
    }

    strclose(istr);
    strclose(ostr);
    free(line);
    line = NULL;
    dprintf(0,"Read %d lines\n",i);
}

