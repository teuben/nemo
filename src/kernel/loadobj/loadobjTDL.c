/*
  File: loadobjTDL.c                                          

  Interface to the more portable dlopen() etc. package        
  As of this writing, totally untested                        
 
  4-mar-06    trying out....                                  

  Uses env.var:     LTDL_LIBRARY_PATH

*/


#include <stdinc.h>

#include "ltdl.h"

static void *dl_handle = NULL;	/* void pointer to current object file */

/***************************************************************/
/* loadobj(pathname);                                          */
/***************************************************************/

void loadobj(string pathname)
{
    char *err;

    dprintf(1,"loadobj: %s\n",pathname);
    dl_handle = lt_dlopen(pathname); 
    err = lt_dlerror();
    if (err != NULL) error("loadobj: error from dlopen: %s",err);
}
/***************************************************************/
/* fn = findfn(fnname);                                        */
/***************************************************************/

proc findfn(string fnname)
{
    proc ret;

    dprintf(1,"findfn: looking up %s\n",fnname);
    ret = (proc) lt_dlsym(dl_handle,fnname);
    return ret;
}



/***************************************************************/
/* mysymbols(progname);                                        */
/***************************************************************/

void mysymbols(string progname)
{
    dprintf(1,"MySymbols: NULL code in loadobjDL\n");

    if (dl_handle) {
        warning("mysymbols: Closing old dlopen");
	return;
    }
    /* at some point I used:
     * #if defined(sgi) || defined(NEED_MAIN_SYMBOLS)
     * here, but it appears to work on SGI now ...
     */

/* 	
 *  It appears symbols from a 'static' binary are not read 
 *  they must come from a shared object 
 *
 * 	including or excluding the following segment doesn't seem
 *	matter 
 *  yuck: need -rdynamic flag in compiling
 */

}
