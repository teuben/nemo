/***************************************************************/
/* File: loadobj.c                                             */
/*                                                             */
/*  STUB stubbed loadobj() package                             */
/*  20-jan-94 renamed from loadobjNULL to loadobjSTUB          */
/***************************************************************/

#include <stdinc.h>
#include <getparam.h>

/***************************************************************/
/* loadobj(pathname);                                          */
/***************************************************************/

void loadobj(pathname)
string pathname;
{
  warning("loadobjSTUB(loadobj): no dynamic object loader");
}
/***************************************************************/
/* fn = findfn(fnname);                                        */
/***************************************************************/

proc findfn(fnname)
string fnname;
{
  error("loadobjSTUB(findfn): no dynamic object loader");
  return NULL;
}



/***************************************************************/
/* mysymbols(progname);                                        */
/***************************************************************/

void mysymbols(progname)
string progname;
{
  warning("loadobjSTUB(findfn): no dynamic object loader");
}

#if defined(TESTBED)
main()
{
  printf("loadobjSTUB: dynamic object loader not implemented\n");
}
#endif
