/***************************************************************/
/* File: loadobj.c                                             */
/*                                                             */
/*  NULL stubbed loadobj() package                             */
/***************************************************************/

#include <stdinc.h>
#include <getparam.h>

/***************************************************************/
/* loadobj(pathname);                                          */
/***************************************************************/

void loadobj(pathname)
string pathname;
{
  warning("loadobjNULL(loadobj): no dynamic object loader");
}
/***************************************************************/
/* fn = findfn(fnname);                                        */
/***************************************************************/

proc findfn(fnname)
string fnname;
{
  error("loadobjNULL(findfn): no dynamic object loader");
}



/***************************************************************/
/* mysymbols(progname);                                        */
/***************************************************************/

void mysymbols(progname)
string progname;
{
  warning("loadobjNULL(findfn): no dynamic object loader");
}

#if defined(TESTBED)
main()
{
  printf("loadobjNULL: dynamic object loader not implemented\n");
}
#endif
