/*
 * KPRINTF.C: store output "parameters" for optional output
 *
 *	real eps = 0.025;
 *	kprintf("eps","%6.3f", eps);
 *
 */

#include <stdinc.h>

#if defined(__STDC__) || defined(__cplusplus)
/* Posix compliant */
#include <stdarg.h>


/*
 * KPRINTF: printf-style messages, parameter lookup done via getparam.c
 *		
 */

void kprintf(string key, const string fmt, ...)
{
    va_list ap;
    string magic;

    if (!isakey(key)) return;

    magic = secret_location(key);

    va_start(ap, fmt);	
    vsprintf(magic, fmt, ap);	/* this should be ok on sun for now */
    va_end(ap);

}
#else

/* Old K&R version ------------------------ needs work */

#include <varargs.h>

/*VARARGS1*/
void kprintf(key, va_alist)
string key;
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
