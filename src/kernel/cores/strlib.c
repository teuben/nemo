
/***************************************************************/
/* File: strlib.c                                              */
/* Last modified on Wed Dec  4 10:52:18 1985 by roberts        */
/* Last modified on Sat Dec  6 15:19:40 1986 by josh           */
/*      22-nov-91   malloc() -> allocate()                     */
/*      25-feb-92   gcc 2.0				       */
/*      20-feb-94   ansi                                       */
/*      20-jun-01   gcc 3                                      */
/* ----------------------------------------------------------- */
/*     The strlib package contains the implementations for     */
/* several routines that use dynamically-allocated string      */
/* storage.  Since these routines tend to fill up memory,      */
/* this package is not suitable for use in applications        */
/* which will run for extended periods or which require        */
/* tight memory control.  Nonetheless, this package does       */
/* provide an easy-to-use interface which is applicable        */
/* to a wide range of programs.                                */
/*                                                             */
/*     The following routines are defined by strlib:           */
/*                                                             */
/*       getmem(nbytes)      malloc with error checking        */
/*       scopy(source)       returns a copy of source          */
/*       sconc(s1,s2)        concatenates its arguments        */
/*       substr(s, p1, p2)   returns substring from p1-p2      */
/*       findstr(text, pat)  finds index of pat in text        */
/***************************************************************/

#include <stdinc.h>
#include <strlib.h>


/***************************************************************/
/* ptr = (type) getmem(nbytes);                                */
/*                                                             */
/*     This routine is exactly like malloc except that (1) it  */
/* checks for no memory errors, and (2) it is defined to take  */
/* an integer rather than an unsigned to keep lint happier.    */
/* See also: allocate()                                        */
/***************************************************************/

char *getmem(int nbytes)
{
    return (char *) allocate(nbytes);
}



/***************************************************************/
/* s = scopy(t);                                               */
/*                                                             */
/*     Copies the string t into dynamically-allocated storage. */
/* see also: strdup()                                          */
/***************************************************************/

string scopy(const_string s)
{
    string result;

    result = (string) getmem((int)strlen(s)+1);
    strcpy(result, s);
    return result;
}



/***************************************************************/
/* s = sconc(s1, s2);                                          */
/*                                                             */
/*     Concatenates two strings and returns the result in      */
/* dynamically-allocated storage.                              */
/***************************************************************/

string sconc(string s1, string s2)
{
    int l;
    string result;

    result = (string) getmem((l = (int)strlen(s1)) + (int)strlen(s2) + 1);
    strcpy(result, s1);
    strcpy(result+l, s2);
    return result;
}



/***************************************************************/
/* s = substr(s, p1, p2);                                      */
/*                                                             */
/*     Returns the substring of s extending from the integer   */
/* indices p1 and p2 (inclusive).  The following edge cases    */
/* apply:                                                      */
/*                                                             */
/*      if p1 < 0 then p1 <- 0;                                */
/*      if p2 > strlen(s) then p2 <- strlen(s);                */
/*      if p1 > p2 then return "";                             */
/***************************************************************/

string substr(string s, int p1, int p2)
{
    int l, i;
    string result;

    l = strlen(s);
    if (p1 < 0) p1 = 0;
    if (p2 >= l) p2 = l - 1;
    if ((l = p2 - p1 + 1) <= 0) return ("");
    result = (string) getmem(l + 1);
    for (i = 0; i < l; i++)
	result[i] = s[p1+i];
    result[l] = 0;
    return result;
}



/***************************************************************/
/* p = findstr(text, pat);                                     */
/*                                                             */
/*     Searches for the string pat in text and returns the     */
/* first index at which it appears, or -1 if no match is       */
/* found.  This function executes a simple compare and         */
/* advance algorithm, and is inappropriate if text contains    */
/* a very long string.                                         */
/***************************************************************/

int findstr(string text, string pat)
{
    register string s;
    int nch;

    nch = strlen(pat);
    for (s = text; *s; s++)
	if (strncmp(s, pat, (unsigned)nch) == 0)
	    return (s - text);
    return -1;
}
