/*
 * DPRINTF.C: print debug messages, if requested debug level high enough
 *
 * xx-Mar-88  Peter Teuben: created
 * 14-Jun-88  Josh Barnes: moved to seperate source file.
 * 10-sep-90  PJT: routines to set_ and get_ the debug level
 *  8-jul-91  pjt: using VFPRINTF now - gcc should work again
 * 12-oct-91  pjt: added __BORLANDC__ to vfprintf users
 * 25-feb-92  pjt: happy gcc2.0
 *  6-apr-93  pjt: took debug_level back into this routine
 * 21-sep-93  pjt: ANSI
 *  4-mar-96  pjt: format string now a 'const'
 *
 *	GCC complains about va_start() being not OK
 *	GCC v1.39 stopped complaining.... but actually doesn't work
 */

#include <stdinc.h>

int debug_level=0;	/* needs to be global; see also getparam.c */

#if defined(__STDC__) || defined(__cplusplus)
/* Posix compliant */
#include <stdarg.h>


/*
 * DPRINTF: printf-style debugging messages, controlled by debug_level
 *	    as set by the user interface (debug=)
 */

void dprintf(int debug, const string fmt, ...)
{
    va_list ap;


    if (debug <= debug_level) {		/* print this debug message? */
        va_start(ap, fmt);	
#if defined(NEED_DOPRNT)
	_doprnt(fmt, ap, stderr);	/*   use low-level interface */
#else
	vfprintf(stderr, fmt, ap);	/* this should be ok on sun for now */
#endif
					/*   may become vprintf in ANSI C */
	fflush(stderr);			/*   drain std error buffer */
        va_end(ap);
    }
}
#else

/* Old K&R version */

#include <varargs.h>

/*VARARGS1*/
void dprintf(debug, va_alist)
int debug;
va_dcl
{
    va_list l;
    string fmt;

    va_start(l);			/* l starts with string */
    fmt = va_arg(l, string);		/* point to string */ /*PPAP*/
					/* l now has remainder */	
    if (debug <= debug_level) {		/* print this debug message? */
#if defined(NEED_DOPRNT)
	_doprnt(fmt, l, stderr);	/*   use low-level interface */
#else
	vfprintf(stderr, fmt, l);	/* this should be ok on sun for now */
#endif
					/*   may become vprintf in ANSI C */
	fflush(stderr);			/*   drain std error buffer */
    }
    va_end(l);
}

#endif
