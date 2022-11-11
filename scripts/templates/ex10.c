/*
 *  NEMO's getparam in C
 */

#include <nemo.h>

string defv[] = {
    "aaa=1\n             some help for aaa",
    "bbb=10.0\n          some help for bbb",
    "VERSION=1.0\n       11-Nov-2022 XYZ",
    NULL,
};

string usage="example10 NEMO getparam in C";


void nemo_main()
{
  string sa = getparam("aaa");
  string sb = getparam("bbb");
  int    ia = getiparam("aaa");
  real   ra = getrparam("bbb");
  printf("string: aaa=%s bbb=%s   int: aaa=%d   real: bbb=%g\n",
	sa,  sb, ia, ra);
}

