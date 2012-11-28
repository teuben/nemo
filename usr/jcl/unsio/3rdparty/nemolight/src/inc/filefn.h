/***************************************************************/
/* File: filefn.h                                              */
/* Last modified on Sat Jan 11 09:26:56 1986 by roberts        */
/* Last modified on Sat Dec 06          1986 by josh           */
/*                  Sun Oct 07 13:23    1990 by Peter          */
/*		ANSI stuff and shorter .h file                 */
/*                  20-sep-93  proper ANSI          pjt        */
/* ----------------------------------------------------------- */
/*     Provides several standard procedures for working with   */
/* files:                                                      */
/*                                                             */
/*         root(filename)                                      */
/*         extension(filename)                                 */
/*         head(filename)                                      */
/*         tail(filename)                                      */
/*         fullname(filename)                                  */
/*         defext(filename, ".xxx")                            */
/*         pathopen(path, filename, mode)                      */
/*         pathfile(path, filename)                            */
/*                                                             */
/* For documentation see filefn.c or manual page filefn.3      */
/***************************************************************/

#ifndef _filefn_h
#define _filefn_h

typedef string (*strfn)(string, string);

extern string root      ( string );
extern string extension ( string );
extern string head      ( string );
extern string tail      ( string );
extern string fullname  ( string );
extern char *defext     ( string, string );
extern stream pathopen  ( string, string, string );
extern string pathfind  ( string, string );
extern string _mappath  ( strfn, string, string, string );

#endif
 



