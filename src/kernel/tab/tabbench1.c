/*
 *     Benchmark some common table I/O operations - see also tabbench2
 *
 *      3-jul-2020  V0.1    drafted
 *     24-jul-2020  V0.2    use getline - Sathvik Ravi
 */

//1    tabgen tab1 1000000   3
//1    tabgen tab2 10000000  3
//1    tabgen tab3 100000000 3
//1    tabgen tab4 3 100000000
//1    /usr/bin/time tabbench1 tab1 .

//2    /usr/bin/time tabtranspose p1M.tab p1Mt.tab $nbody
//2    /usr/bin/time tabbench1 p1Mt.tab .


#include <stdinc.h>
#include <ctype.h>
#include <getparam.h>
#include <table.h>

string defv[] = {
    "in=???\n	       input file",
    "out=???\n         output file",
    "mode=1\n          Benchmark mode (not used yet)",
    "nmax=10000\n      Default max number of lines (in a pipe)",
    "VERSION=0.3\n     25-jul-2020 PJT",
    NULL,
};

string usage="table I/O benchmark";

#ifndef MAX_LINELEN 
#define MAX_LINELEN  2048
#endif

void nemo_main(void)
{
    stream istr, ostr;
    char *line;
    int *select = NULL;
    int nout, next = 0, counter = 0;
    size_t i;
    int nmax = getiparam("nmax");
    string input = getparam("in");
    string output = getparam("out");
    size_t buffer_size = MAX_LINELEN;

    line = malloc((MAX_LINELEN) * sizeof(char));
    nmax = nemo_file_lines(input,nmax);

    dprintf(0,"MAX_LINELEN=%d\n",MAX_LINELEN);
    dprintf(1, "Input File: %s\n", input);
    dprintf(1, "Output File: %s\n", output);

    istr = stropen(getparam("in"),"r");
    ostr = stropen(getparam("out"),"w");

    i = 0;
    while (getline(&line, &(buffer_size), istr) != -1) {
#if 1
        counter = 0;
        while(isspace(line[counter]) != 0) {
            ++counter;
        }
        if (line[counter] != '#' && line[counter] != '/' && line[counter] != '!') {
	    i++;
	    fputs(line,ostr);
        }
#else	
	i++;
#endif
#if 1
	fputs(line,ostr);
#endif
    }

    strclose(istr);
    strclose(ostr);
    
    free(line);
    line = NULL;
    dprintf(0,"Read %ld lines\n",i);
    dprintf(0,"Longest line: %d\n",buffer_size);
}

