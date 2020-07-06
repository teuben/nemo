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

string usage="table I/O benchmark";

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

