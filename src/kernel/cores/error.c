/*
 * ERROR: scream (on stderr) and die (quickly?).
 * WARNING: scream  (on stderr)
 * RECOVER: set a recover() function when error() called
 * STOP:    exit to system with proper shell status (does not return)
 *
 *      See: dprintf() for a similar teachnique that issues a
 *      warning message and then continues.
 *
 *      xx-xxx-86  - origal program for Nemo                    JEB
 *      xx-xxx-88  - proper <varargs.h>
 *      12-sep-90  - first merged NEMO and Starlab program      PJT
 *	16-nov-90  - no #ifdef MICRO anymore			PJT
 *	22-mar-91  - ERROR env.var. introduced - via stop()	PJT
 *	 8-jul-91 V1.2  using vfprintf()                        pjt
 *	12-oct-91  - added __BORLANDC__ to the vfprintf() users pjt
 *	25-feb-92    made all of them void'd			pjt
 *	 6-apr-93 V1.3 moved stop() in here, plus local data    pjt
 *	20-sep-93 V1.4 ansi requires > 0 named arguments 	pjt
 *	26-feb-94 V1.5 keep old <varags> code around too	pjt
 *      21-mar-01 V1.6 removed old <varargs> code, only allow ANSI compilers 	PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <stdarg.h>


extern int debug_level;    /* see also user interface getparam.c for this */

int error_level = 0;	   /* yuck, but this needs to be globally visible */
local proc cleanup=NULL;   /* if not NULL this points to error cleanup proc */
local int error_count=0;

void error(string fmt, ...)
{
    va_list ap;

    fprintf(stderr,"### Fatal error [%s]: ",getargv0());  /* report name */
    
    va_start(ap, fmt);              /* ap starts with string 'fmt' */

#if defined(NEED_DOPRNT)
    _doprnt(fmt, ap, stderr);       /* print out on stderr */
#else
    vfprintf(stderr, fmt, ap);      /* print out on stderr */
#endif

    if (fmt[strlen(fmt)-1] != '\n') /* be nice if no newline supplied */
        fprintf(stderr,"\n");       /* and print it anyhow */
    fflush(stderr);                 /* flush it NOW */
    va_end(ap);                      /* end varargs */
    
    if (cleanup==NULL){		    /* if no cleanup set */
      if (debug_level>5) abort();   /* produce coredump if requested */
      stop(-1);                     /* close the shop and pass this to parent */
    } else {                        /* else say so : */
	fprintf(stderr,"### Recoverable error ....\n");
	(*cleanup)();               /* clean up */
    }                               /* and proceed as if nothing happened */
}

void fatal(string fmt, ...)
{
    va_list ap;

    fprintf(stderr,"### Fatal error [%s]: ",getargv0());  /* report name */
    
    va_start(ap, fmt);              /* ap starts with string 'fmt' */

#if defined(NEED_DOPRNT)
    _doprnt(fmt, ap, stderr);       /* print out on stderr */
#else
    vfprintf(stderr, fmt, ap);      /* print out on stderr */
#endif

    if (fmt[strlen(fmt)-1] != '\n') /* be nice if no newline supplied */
        fprintf(stderr,"\n");       /* and print it anyhow */
    fflush(stderr);                 /* flush it NOW */
    va_end(ap);                     /* end varargs */
    abort();                        /* nasty, but writes a core dump and dies */
}

void warning(string fmt, ...)
{
    va_list ap;

    fprintf(stderr,"### Warning [%s]: ",getargv0());  /* report name */

    va_start(ap, fmt);               /* ap starts with string 'fmt' */

#if defined(NEED_DOPRNT)
    _doprnt(fmt, ap, stderr);        /* print out on stderr */
#else
    vfprintf(stderr, fmt, ap);       /* print out on stderr */
#endif

    if (fmt[strlen(fmt)-1] != '\n') /* be nice if no newline supplied */
	fprintf(stderr,"\n");       /* and print it anyhow */
    fflush(stderr);                 /* flush it NOW */
    va_end(ap);                      /* end varargs */
}



void recover(cl)
proc cl;
{
    if (cl)
	dprintf(1,"Setting recoverable error\n");
    else
	dprintf(1,"Resetting recoverable error\n");
    cleanup = cl;
}


void stop(lev)
int lev;
{
    if (lev<0)
        if (error_count++ < error_level) {
            warning("[%d/%d] error ignored",error_count,error_level);
            return;
        }
    finiparam();
    exit(lev);
    /*NOTREACHED*/
}


#ifdef TESTBED
string defv[]={
    "recover=f\n    Try a recoverable error",
    "count=1\n	    Loopcount calling error",
    "VERSION=1.5\n  26-feb-94 PJT",
    NULL,
};

nemo_main()
{
    int fun(), n;

    if (getbparam("recover"))
        recover( (proc) fun );
    n = getiparam("count");
    while(n-- > 0)
      error("error: foo=%f  bar=%d  fum=\"%s\"", 3.1415, 32768, "waldo");
}

fun()
{
    printf("Fun Fun Fun - a recoverable error occurred\n");
}
#endif

