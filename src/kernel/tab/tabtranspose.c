/*
 * TABTRANSPOSE: transpose a table the lazy way (the whole table in memory)
 *
 *
 *  @todo    burststring() is used, which limits # columns to 2048.  Use strtok() ?
 */

#include <stdinc.h>
#include <getparam.h>
#include <extstring.h>
#include <table.h>

string defv[] = { 
    "in=???\n           input file name(s)",
    "out=???\n          output file name",
    "align=f\n          align the columns?",
    "nmax=0\n           max space (needed if data in pipe)",
    "VERSION=1.1\n      24-jul-2020 PJT",
    NULL
};

string usage = "transpose a table";

string  input, output;			/* file names */
bool alignment;                         /* align the columns ?*/
stream  instr, outstr;			/* file streams */
int     nmax;                           /* # lines in file */
int     kmin;                           /* # columns to transpsse */

string *lines;               /* pointer to all lines */
string **words;              /* pointer to all words */

local void setparams(void);
local void do_work(void);
local void do_output(void);

extern  string *burststring(string, string);


void nemo_main(void)
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

    alignment = getbparam("align");

    nmax = nemo_file_lines(input,getiparam("nmax"));
    
    lines = (string *) allocate (nmax * sizeof(string));
    words = (string **) allocate (nmax * sizeof(string *));
}

local void do_work(void)
{
  int k, n=0;

  kmin = -1;    /* keep track of min number of columns */

  while ( (lines[n] = getaline(instr)) != NULL) {
    if (n == nmax) {
      warning("Too many lines, change nmax=");
      break;
    }
    if (*lines[n] == '#') continue;
    words[n] = burststring(lines[n]," ,\t");
    k = xstrlen(words[n],sizeof(string))-1;
    kmin = (kmin < 0 ?  k :  MIN(kmin,k));
    dprintf(2,"%d: %s\n",k, lines[n]);
    n++;
  }

  nmax = n;
  dprintf(0,"Read %d lines, %d columns to be transposed\n",n,kmin);
}

local void do_output(void)
{
  int i, j, max_spaces = 0, count_spaces;
  int max_space[nmax];

  if (alignment) {
    dprintf(1, "The alignment parameter is true\n"); 
  } else {
    dprintf(1, "The alignment parameter is false\n");
  }

  if (alignment) {
    for (j = 0; j < nmax; ++j) {
      for (i = 0; i < kmin; ++i) {
        if ((int)strlen(words[j][i]) > max_spaces) {
          max_spaces = (int)strlen(words[j][i]);
        }
      }
    /* fprintf(outstr, "%d\n", max_spaces); */
      max_space[j] = max_spaces;
      max_spaces = 0;
    }
  }

  /* for (j = 0; j < nmax; ++j) {
    fprintf(outstr, "%d\n", max_space[j]);
  } */

  for (i=0; i<kmin; i++) {
    for (j=0; j<nmax; j++) {
      fprintf(outstr,"%s ",words[j][i]);
      if (alignment) {
        for (count_spaces = 0; count_spaces < max_space[j] - (int)strlen(words[j][i]); ++count_spaces) {
          fprintf(outstr, " ");
        }
      }
    }
    fprintf(outstr,"\n");
  }
}


