/***************************************************************/
/* File: loadobj.c                                             */
/*                                                             */
/* Interface to the GNU dld library                            */
/* See also:  http://www-swiss.ai.mit.edu/~jaffer/DLD.html     */   
/*                                                             */
/* 18-nov-90   Created  - dld V 3.2.1    Peter Teuben          */ 
/* 16-jan-94   Confirmed with dld V 3.2.4 on linux  	PJT    */
/* 13-sep-95   added april95 patches from Stephen Levine       */
/*             to fix multiple loads of same code (dld3.2.3)   */
/*                                                             */
/***************************************************************/

#include <stdinc.h>
#include <dld.h>	/* we keep a copy in  $NEMOINC  */


/***************************************************************/
/* loadobj(pathname);                                          */
/***************************************************************/

void loadobj(string pathname)
{
    int ner;
    
    dprintf(3,"loadobj: %s\n",pathname);
    ner = dld_link(pathname);
    if (ner == DLD_EMULTDEFS) {
        dld_perror("EMULTDEFS Error in loadobj");
    } else if (ner != 0) {
        dld_perror("Error in loadobj");
        error("Could not load %s",pathname);
    }
}
/***************************************************************/
/* fn = findfn(fnname);                                        */
/***************************************************************/

proc findfn(fnname)
string fnname;
{
    proc ret;

    dprintf(3,"findfn: looking up %s\n",fnname);
    ret = (proc) dld_get_func(fnname);
    return(ret);
}



/***************************************************************/
/* mysymbols(progname);                                        */
/***************************************************************/

void mysymbols(progname)
string progname;
{
    dprintf(3,"mysymbols: setting up %s\n",progname);
    if (dld_init(progname) != 0) {
        dld_perror("Error initializing dld_init");
        error("mysymbols(%s)",progname);
    }
}
