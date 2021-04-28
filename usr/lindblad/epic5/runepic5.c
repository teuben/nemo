/*
 *  RUNEPIC5:    possible frontend for epic5 
 */

#include <nemo.h>

string defv[] = {
    "vfile=vrot\n    some help",
    "bfile=perpot\n    some help",
    "omegap=-18\n    some help",
    "blambd0=8\n    some help",
    "blambdf=7\n    some help",
    "eomeg=5\n    some help",
    "r0=0.01\n    some help",
    "dr=0.5\n    some help",
    "rmax=22\n    some help",
    "ntheta=32\n    some help",
    "VERSION=0.0\n       21-Apr-2021 XYZ",
    NULL,
};

string usage="TEMPLATE program usage -- fill in yourself";

string cvsid="$Id:$";

void nemo_main()
{
  /* your code goes here */
  /* e.g.   string s = getparam("s");  */
  /* e.g.   int    i = getiparam("i"); */
  /* e.g.   real   r = getrparam("r"); */
  /* e.g.   double d = getdparam("d"); */
  /* followed by lots of lovely calculations */
  warning("New template NEMO program");

    string vfile = getparam("vfile");
    string bfile = getparam("bfile");
    string omegap = getparam("omegap");
    string blambd0 = getparam("blambd0");
    string blambdf = getparam("blambdf");
    string eomeg = getparam("eomeg");
    string r0 = getparam("r0");
    string dr = getparam("dr");
    string rmax = getparam("rmax");
    string ntheta = getparam("ntheta");

  /*
   *  other useful things:
   *    if (hasvalue("xxx") && !hasvalue("yyy"))
   *        error("need value for yyy");
   */


  // write the 'epar' file

  // write the 'epic5.stdin' file

  // system('epic5 < epic5.stdin')
}

