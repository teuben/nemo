/*
 * TABTRANSPOSE: transpose a matrix the lazy way
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <table.h>
#include <extstring.h>
#include <ctype.h>

string defv[] = { 
    "in=???\n           input file name(s)",
    "out=???\n          output file name",
    "nmax=10000\n       max space (needed if data in pipe)",
    "VERSION=1.0\n      5-oct-02 PJT",
    NULL
};

string usage = "transpose a table";

string  input, output;			/* file names */
stream  instr, outstr;			/* file streams */
int     nmax;                           /* # lines in file */
int     kmin;                           /* # columns to transpsse */

string *lines;               /* pointer to all lines */
string **words;              /* pointer to all words */

local void setparams(void), do_work(), do_output();

extern  string *burststring(string, string);


nemo_main()
{
    setparams();
    do_work();
    do_output();
}

local void setparams(void)
{
    int i;

    input = getparam("in");
    instr = stropen(input,"r");

    output = getparam("out");
    outstr = stropen(output,"w");

    nmax = nemo_file_lines(input,getiparam("nmax"));
    lines = (string *) allocate (nmax * sizeof(string));
    words = (string **) allocate (nmax * sizeof(string *));
}

local void do_work()
{
  int k, n=0;

  kmin = -1;    /* keep track of min number of columns */

  while ( (n < nmax) && (lines[n] = getaline(instr)) != NULL) {
    if (*lines[n] == '#') continue;
    words[n] = burststring(lines[n]," ,\t");
    k = xstrlen(words[n],sizeof(string))-1;
    kmin = (kmin < 0 ?  k :  MIN(kmin,k));
    dprintf(5,"%d: %s\n",k, lines[n]);
    n++;
  }
  nmax = n;
  dprintf(0,"Read %d lines, %d columns to be transposed\n",n,kmin);
}

local void do_output()
{
  int i, j;

  for (i=0; i<kmin; i++) {
    for (j=0; j<nmax; j++)
      fprintf(outstr,"%s ",words[j][i]);
    fprintf(outstr,"\n");
  }
}
