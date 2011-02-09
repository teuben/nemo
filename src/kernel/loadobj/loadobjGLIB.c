/***************************************************************/
/* File: loadobj.c                                             */
/*                                                             */
/*  GLIB based oadobj() package                                */
/***************************************************************/

#include <stdinc.h>
#include <getparam.h>
#include <gmodule.h>

/***************************************************************/
/* loadobj(pathname);                                          */
/***************************************************************/

void loadobj(string pathname)
{
  warning("loadobjGLIB(loadobj): no dynamic object loader yet");
}


/***************************************************************/
/* fn = findfn(fnname);                                        */
/***************************************************************/

proc findfn(string fnname)
  error("loadobjGLIB(findfn): no dynamic object loader");
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
