/*
 * (NEMO_)DPRINTF.C: print debug messages, if requested debug level enough
 *
 *	Note dprintf(int,char *,...) is now used by the snprintfv package
 *	see http://www.gnu.org/software/autogen/
 *	MacOS is also known to have this function, with exactly the same
 *	purpose as NEMO's original dprintf()
 *	The official name for NEMO's function is now nemo_dprintf()
 *	although a macro will accept dprintf() since their prototype
 *	*happens* to be the same as ours. You just don't have acccess
 *	to their code.
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
 * 20-jun-01  pjt: gcc0, no more non-ansi code
 * 26-sep-01  name now nemo_dprintf(), to avoid conflict with new dprintf(3?)
 *  8-may-04  return -1 if debug level not high enough
 *  4-may-05  nemo_set_debug, why was this not needed before?
 * 22-feb-06  WD added blurb ("### nemo Debug info: ") before print out
 *  5-aug-06  WD print_head: "### ... " printed not more than once in a line
 * 12-jun-08  WD nemo_dprintf is macro in stdinc.h, report [file:line]
 */

#include <stdinc.h>
#include <stdarg.h>

int debug_level=0;	/* needs to be global; see also getparam.c */

static char *nemo_file = "dprintf.c: debugging stuff";

bool nemo_debug(int debug) 
{
  return debug <= debug_level;
}

void nemo_set_debug(int debug) 
{
  debug_level = debug;
}

/* START changes WD 12th June 2008 */
/* code to support reporting of library name, file name, and line number */
static const_string debug_file = 0;
static int          debug_line = 0;
#ifndef   DEBUG_LEVEL_FOR_REPORTING_FILE
#  define DEBUG_LEVEL_FOR_REPORTING_FILE 4
#endif
/* code to support reporting of MPI process number */
bool mpi_proc = 0;
int  mpi_rank = 0;
void set_mpi_rank(int rank)
{
  mpi_proc = 1;
  mpi_rank = rank;
}
/* END   changes WD 12th June 2008 */

/*
 * DPRINTF: printf-style debugging messages, controlled by debug_level
 *	    as set by the user interface (debug=)
 */

int __nemo_dprintf(int debug, const_string fmt, ...)
{
    va_list ap;
    int nret = -1;
    static bool print_head = TRUE;

    if (debug <= debug_level) {		/* print this debug message? */

        /* START changes WD 12th June 2008 */
        if(print_head) {
	    if(mpi_proc)
	        fprintf(stderr,"### nemo Debug Info @%d: ",mpi_rank);
	    else
	        fprintf(stderr,"### nemo Debug Info: ");
	    if(debug_file && debug_level >= DEBUG_LEVEL_FOR_REPORTING_FILE)
	        fprintf(stderr,"[%s:%d]: ",debug_file,debug_line);
	}
	/* END   changes WD 12th June 2008 */

        va_start(ap, fmt);	
#if defined(NEED_DOPRNT)
	_doprnt(fmt, ap, stderr);	/*   use low-level interface */
#else
	nret= vfprintf(stderr, fmt, ap);/*   this should be ok on sun for now */
#endif
					/*   may become vprintf in ANSI C */
	fflush(stderr);			/*   drain std error buffer */
        va_end(ap);
	print_head = fmt && fmt[strlen(fmt)-1] == '\n';
    }
    return nret;
}

/* START changes WD 12th June 2008 */
dprintf_pter get_dprintf(const_string file, int line)
{
    debug_file = file;
    debug_line = line;
    return &__nemo_dprintf;
}
/* END   changes WD 12th June 2008 */

