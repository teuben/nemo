/****************************************************************************/
/* GETPARAM.H: include file for command-line processing.                    */
/* Copyright (c) 1999 by Joshua E. Barnes, Tokyo, JAPAN.                    */
/****************************************************************************/

#ifndef _getparam_h
#define _getparam_h

void initparam(string *, string *);

string getparam(string);

int getiparam(string);
bool getbparam(string);
double getdparam(string);

/*
 * Macros to obtain name and version of program.
 */

#define getargv0()    (getparam("argv0"))
#define getversion()  (getparam("VERSION"))

/*
 * Function to inquire about parameter status.
 */

int getparamstat(string);

/* Bits for parameter status flags. */

#define DEFPARAM        001                     /* param has default value  */
#define REQPARAM        002                     /* must be given a value    */
#define ARGPARAM        004                     /* value reset by argument  */
#define INPARAM         010                     /* param used for input     */
#define OUTPARAM        020                     /* param used for output    */
#define HIDPARAM        040                     /* parameter is hidden      */

#endif  /* ! _getparam_h */
