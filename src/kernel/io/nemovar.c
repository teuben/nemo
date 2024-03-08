/*
 *  NEMOVAR:   display a NEMO variable from the $NEMOVAR file (or future shared memory)
 */

#include <nemo.h>

string defv[] = {
  "var=???\n          Name of  variable(s)",
  "value=\n           Optional new value",
  "VERSION=0.1\n      7-mar-2024 PJT",
  NULL,
};

string usage="display NEMO variable(s)";


void nemo_main()
{
  string var = getparam("var");
  string value = getparam("value");  // could be blank
  bool  Qset = hasvalue("value");
  string nemovar = getenv("NEMOVAR");

  dprintf(0,"nemovar: %s\n", nemovar);

  if (Qset) {
    warning("setting not allowed yet");
    return;
  }
  
  printf("%s= TBD\n",var);

}
