/***************************************************************/
/* File: loadobj.c                                             */
/*                                                             */
/* Interface to the so called standard dlopen() etc. package   */
/* SUN OS 4.1 has introduced 				       */
/*                                                             */
/* There were reported problems that dlclose() didn't work     */
/* Indeed, I never call dlclose() now - but it means you can't */
/* cross-call routines,....                                    */
/* Also, when the routine contains initialized data, a .sa     */
/* should be present in addition to the .so file               */
/* This is not implemented yet....			       */
/*                                                             */
/* 10-sep-90   Created         Peter Teuben                    */ 
/* 20-apr-92   Made it working on SUN OS 4.1.2 :               */
/*		cc -pic -c ld-tmp.c                            */
/*              ld -o ld-tmp.so ld-tmp.o                       */
/* 20-jan-94   only works for SUN OS 5.x now 		       */
/*             but also works for ELF linux                    */ 
/*                                                             */
/***************************************************************/

#include <stdinc.h>

#include <dlfcn.h>

static int *dl_handle = NULL;	/* void pointer to current object file */

/***************************************************************/
/* loadobj(pathname);                                          */
/***************************************************************/

void loadobj(string pathname)
{
    char *err;

    dprintf(1,"loadobj: %s\n",pathname);
#if 0
    if (dl_handle) {
        warning("loadobj: Closing old dlopen");
        if (dlclose(dl_handle)) 
            warning("loadobj: Some error closing old handle");
    }
#endif
    dl_handle = dlopen(pathname,1);
    err = dlerror();
    if (err != NULL) error("loadobj: error from dlopen: %s",err);
}
/***************************************************************/
/* fn = findfn(fnname);                                        */
/***************************************************************/

proc findfn(string fnname)
{
    proc ret;

    dprintf(1,"findfn: looking up %s\n",fnname);
    ret = (proc) dlsym(dl_handle,fnname);
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
        if (dlclose(dl_handle)) 
            warning("mysymbols: Some error closing old handle");
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

#if 0
    dl_handle = dlopen(NULL,1);
    err = dlerror();
    if (err != NULL) error("mysymbols: error from dlopen: %s",err);
#endif
}
