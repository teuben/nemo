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
 *      28-nov-01 V1.7 allow to set an exit level               pjt
 *       8-dec-01 V1.8 added errno reporting	 		pjt
 *      16-jan-01 V1.8a  calling abort() will be announced      pjt
 *      13-feb-03 V1.8b  revoking errno reporting		pjt
 *      06-jun-08 V1.8c  nemo_exit                              wd
 *      12-jun-08 V1.8d  report MPI proc                        wd
 *      16-sep-08 V1.8e  removed nemo_exit (see stdinc.h)       wd
 */

#include <stdinc.h>
#include <getparam.h>
#include <stdarg.h>
#include <errno.h>

extern int debug_level;    /* see also user interface getparam.c for this */

int error_level = 0;	   /* yuck, but this needs to be globally visible */
local proc cleanup=NULL;   /* if not NULL this points to error cleanup proc */
local int error_count=0;
local int exit_level=0;
local int last_errno=0;

static void report_errno(void)
{
#if 0
	/* something wrong here, this thing is lying half the time  */
    if (errno)
        fprintf(stderr,"### Fatal errno %d errmsg=%s\n",
            errno, strerror(errno));
#endif
}

/* 
 * commented out WD 10-09-2008, also in stdinc.h
void errorn(string fmt, ...)
{
    error("errorn is sadly not implemented; use debug>0 to get errno messages");
}
*/

/* Start changes WD 12/06/2008 */
extern bool mpi_proc;   /* dprintf.c */
extern int  mpi_rank;   /* dprintf.c */
/* End changes WD 12/06/2008 */

void error(string fmt, ...)
{
    va_list ap;

    report_errno();
    fprintf(stderr,"### Fatal error [%s]: ",getargv0());  /* report name */

    /* Start changes WD 12/06/2008 */
    if(mpi_proc)
      fprintf(stderr,"@%d: ",mpi_rank);
    /* End changes WD 12/06/2008 */
    
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
      if (debug_level>5) {          /* produce coredump if requested */
        fprintf(stderr,"Now aborting....\n");
	fflush(stderr);                
	abort();
      }
      if (exit_level)
        stop(exit_level);
      else
        stop(-1);                   /* close the shop and pass this to parent */
    } else {                        /* else say so : */
	fprintf(stderr,"### Recoverable error ....\n");
	(*cleanup)();               /* clean up */
    }                               /* and proceed as if nothing happened */
}

void fatal(string fmt, ...)
{
    va_list ap;

    report_errno();
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
    fprintf(stderr,"Now aborting....\n");
    fflush(stderr);                
    abort();                        /* nasty, but writes a core dump and dies */
}

void warning(string fmt, ...)
{
    va_list ap;

    fprintf(stderr,"### Warning [%s]: ",getargv0());  /* report name */

    /* Start changes WD 12/06/2008 */
    if(mpi_proc)
      fprintf(stderr,"@%d: ",mpi_rank);
    /* End changes WD 12/06/2008 */

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



void recover(proc cl)
{
    if (cl)
	dprintf(1,"Setting recoverable error\n");
    else
	dprintf(1,"Resetting recoverable error\n");
    cleanup = cl;
}


void stop(int lev)
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

void set_exit_level(int lev)
{
    exit_level = lev;
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

