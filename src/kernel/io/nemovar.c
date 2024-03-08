/*
 *  NEMOVAR:   display a NEMO variable from the $NEMOVAR file (or future shared memory)
 */

#include <nemo.h>
#include <table.h>
#include <extstring.h>


string defv[] = {
  "var=\n             Name of  variable(s)",
  "value=\n           Optional new value",
  "show=f\n           Show all variables",
  "VERSION=0.2\n      7-mar-2024 PJT",
  NULL,
};

string usage="display NEMO variable(s)";


void nemo_main()
{
  string var = getparam("var");
  string value = getparam("value");  
  bool  Qset = hasvalue("value");
  bool Qshow = getbparam("show");
  string nemovar = getenv("NEMOVAR");
  stream tstr = stropen(nemovar,"r");
  table *t = table_open(tstr, -1);      // read table as a dumb set of lines
  int i, nrows = t->nr;
  string s;

  dprintf(1,"nemovar: %s\n", nemovar);
  

  if (Qshow) {
    for (i=0; i<nrows; i++)
	printf("%s\n", table_row(t, i));
    return;
  }



  if (Qset) { 
    // writing code might need flock() to lock file
    warning("setting not allowed yet");
    return;
  }

  // loop over lines and find the keyword, then display it
  
  printf("%s= TBD\n",var);

}
