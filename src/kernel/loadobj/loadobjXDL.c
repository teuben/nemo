/***************************************************************/
/* File: loadobj.c                                             */
/*                                                             */
/* Interface to the xdl package                                */
/* 	(Gary Flake, peyote@umiacs.umd.edu)                    */
/*                                                             */
/* 28-oct-93   Created         Peter Teuben                    */ 
/*                                                             */
/***************************************************************/

#include <stdinc.h>
#include <getparam.h>
#include <xdl.h>	/* NEMOINC version !!! */

static int loadobj_to_xdl = 0;

/***************************************************************/
/* loadobj(pathname);                                          */
/***************************************************************/

void loadobj(pathname)
string pathname;
{
    string objs[2];

    if (loadobj_to_xdl==0)
        error("loadobj: not properly initialized; call mysymbols()");

    dprintf(1,"loadobj: %s\n",pathname);

    objs[0] = pathname;
    objs[1] = 0;
    if (xdl_link(objs))
        error("loadobj: xdl_link(%s) failed",pathname);    
}
/***************************************************************/
/* fn = findfn(fnname);                                        */
/***************************************************************/

proc findfn(fnname)
string fnname;
{
    proc ret;

    if (loadobj_to_xdl==0)
        error("loadobj: not properly initialized; call mysymbols()");

    dprintf(1,"findfn: looking up %s\n",fnname);
    ret = (proc) xdl_symbol(fnname);
    return ret;
}



/***************************************************************/
/* mysymbols(progname);                                        */
/***************************************************************/

void mysymbols(progname)
string progname;
{
    char *err;

    dprintf(1,"MySymbols: NULL code in loadobjDL\n");


    if (xdl_init(getargv0()))
        error("loadobj: xld_init failed");

    loadobj_to_xdl = 1;
}
