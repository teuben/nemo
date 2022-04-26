/*
 * TABTRANSPOSE: transpose a table the lazy way
 *            (the whole table in memory - is there even another way?)
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <extstring.h>
#include <table.h>

string defv[] = { 
    "in=???\n           input file name",
    "out=???\n          output file name",
    "align=f\n          align the columns?",
    "VERSION=2.0\n      26-apr-2022 PJT",
    NULL
};

string usage = "transpose a table, optional align columns";

string  input, output;			/* file names */
table   *tptr;                          /* table pointer */
bool    alignment;                      /* align the columns ?*/
stream  instr, outstr;			/* file streams */
int     nrow, ncol;                     /* # rows/cols in input file */
string **words;                         /* pointer to all words in all lines  [cells???] */

local void setparams(void);
local void do_work(void);
local void do_output(void);

void nemo_main(void)
{
    setparams();
    do_work();
    do_output();
}

local void setparams(void)
{
    input = getparam("in");
    instr = stropen(input,"r");
    tptr = table_open(instr,0);

    output = getparam("out");
    outstr = stropen(output,"w");

    alignment = getbparam("align");

    nrow = table_nrows(tptr);
    ncol = table_ncols(tptr);
    words = (string **) allocate (nrow * sizeof(string *));
}

local void do_work(void)
{
  int nc;
  
  for (int i=0; i<nrow; i++) {
    words[i] = table_rowsp(tptr,i);
    nc = xstrlen(words[i],sizeof(string))-1;
    if (nc != ncol) error("not a uniform table");
  }
}

local void do_output(void)
{
  int i, j, max_spaces = 0, count_spaces;
#if 0
  int max_space[nrow];   // @todo  stacks are limited (wow, this is allowed now)
  int freemp=0;
  // this segfaults:   tabgen - 10000000 2  | tabtranspose - .
  // stack 8MB according to ulimit -a
#else
  int freemp=1;  
  int *max_space = (int *) allocate(nrow*sizeof(int));
#endif
  
  if (alignment) {
    dprintf(1, "The alignment parameter is true\n"); 
  } else {
    dprintf(1, "The alignment parameter is false\n");
  }

  if (alignment) {
    for (j = 0; j < nrow; ++j) {
      for (i = 0; i < ncol; ++i) {
        if ((int)strlen(words[j][i]) > max_spaces) {
          max_spaces = (int)strlen(words[j][i]);
        }
      }
      max_space[j] = max_spaces;
      max_spaces = 0;
    }
  }

  for (i=0; i<ncol; i++) {
    for (j=0; j<nrow; j++) {
      fprintf(outstr,"%s ",words[j][i]);
      if (alignment) {
        for (count_spaces = 0; count_spaces < max_space[j] - (int)strlen(words[j][i]); ++count_spaces) {
          fprintf(outstr, " ");
        }
      }
    }
    fprintf(outstr,"\n");
  }
  if (freemp) free(max_space);
}


