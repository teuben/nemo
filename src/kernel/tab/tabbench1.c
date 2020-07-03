/*
 *     Benchmark some common table operations
 *
 *      3-jul-2020  V0.1    drafted
 */

//    nbody=1000000      (for 10M there is a bug)
//    mkplummer - $nbody | snapprint - > p1M.tab format=%20.16f
//    tabtranspose p1M.tab p1Mt.tab $nbody

#include <stdinc.h>
#include <getparam.h>

string defv[] = {
    "in=???\n	       input file",
    "out=???\n         output file",
    "nmax=10000\n      Default max allocation",
    "mode=1\n          Benchmark mode",
    "VERSION=0.1\n     3-jul-2020 PJT",
    NULL,
};


string usage="table benchmark";

#ifndef MAX_LINELEN 
#define MAX_LINELEN  2048
#endif

void nemo_main(void)
{
    stream istr, ostr;
    char line[MAX_LINELEN];
    int nmax,  *select = NULL;
    int nout, next = 0;
    int    i, j;
    string iname = getparam("in");

    dprintf(0,"MAX_LINELEN=%d\n",MAX_LINELEN);

    istr = stropen(getparam("in"),"r");
    ostr = stropen(getparam("out"),"w");

    i = 0;
    while (fgets(line,MAX_LINELEN,istr) != NULL) {
        i++;
	fputs(line,ostr);
    }
    strclose(istr);
    strclose(ostr);
    dprintf(0,"Read %d lines\n",i);
}

