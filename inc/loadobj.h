/***************************************************************/
/* File: loadobj.h                                             */
/* Last modified on Sat Nov 23 08:54:34 1985 by roberts        */
/* 22-feb-94  ansi header PJT                                  */
/* ----------------------------------------------------------- */
/*     This package is used to implement dynamic loading of    */
/* functions from object files.  Any references in the         */
/* object file to previously defined objects are allowed,      */
/* but no additional searching or unresolved references        */
/* are allowed.  The contents of the package are:              */
/*                                                             */
/*     loadobj(pathname)   -- loads object file into memory    */
/*     findfn(fnname)      -- return function pointer          */
/*     mysymbols(progname) -- declare current symbols          */
/***************************************************************/

#ifndef _loadobj_h
#define _loadobj_h

/* from loadobj.c: */
extern void loadobj   (string);
extern proc findfn    (string);
extern void mysymbols (string);

/* from mapsys.c: */
extern void mapsys    (string);

#endif
