/***************************************************************/
/* File: strlib.h                                              */
/* Last modified on Wed Dec  4 11:01:25 1985 by roberts        */
/*		20-sep-93  ANSI			pjt            */
/* ----------------------------------------------------------- */
/*     The strlib.h file defines the interface for several     */
/* routines that use dynamically-allocated string storage.     */
/* Since these routines tend to fill up memory, this package   */
/* is not suitable for use in applications which will run for  */
/* extended periods or which require tight memory control.     */
/*                                                             */
/* Contents:                                                   */
/*                                                             */
/*       getmem(nbytes)      malloc with error checking        */
/*       scopy(source)       returns a copy of source          */
/*       sconc(s1,s2)        concatenates its arguments        */
/*       substr(s, p1, p2)   returns substring from p1-p2      */
/*       findstr(text, pat)  finds index of pat in text        */
/***************************************************************/

#ifndef _strlib_h
#define _strlib_h

extern char *getmem    ( int nbytes );
extern string scopy    ( const char *s );
extern string sconc    ( char *s1, char *s2 );
extern string substr   ( char *s, int p1, int p2 );
extern int findstr     ( char *text, char *pat );

#endif
