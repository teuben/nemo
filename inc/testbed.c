/*
 *
 *  You can include this in subroutines that have no associated
 *  TESTBED section, but would need one for automatic TESTBED's
 *  as e.g. in the GNU makefiles
 *
 *  #include <testbed.c>
 *
 *  5-nov-93    Created                 peter teuben
 *
 */
 
#ifdef TESTBED

#include <nemo.h>

string defv[] = {
    "message=\n     Hello World!",
    "VERSION=1.0\n  5-nov-93",
    NULL,
};

string usage="blank TESTBED for autocompilers";

nemo_main()
{
    if (hasvalue("message")) printf("%s\n",getparam("message"));
    warning("Only blank TESTBED from $NEMOINC available");
}
#endif
