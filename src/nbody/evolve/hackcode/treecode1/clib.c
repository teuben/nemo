/****************************************************************************/
/* CLIB.C: assorted C routines with prototypes in stdinc.h.                 */
/* Copyright (c) 1999 by Joshua E. Barnes, Tokyo, JAPAN.                    */
/****************************************************************************/

#include "stdinc.h"
#include "getparam.h"
#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <stdarg.h>

/*
 * ALLOCATE: memory allocation, with error checking.
 */

void *allocate(int nb)
{
    void *mem;

    mem = calloc(nb, 1);                /* allocate, also clearing memory   */
    if (mem == NULL)
        error("allocate in %s: not enuf memory (%d bytes)\n",
              getargv0(), nb);
    return (mem);
}

/*
 * CPUTIME: compute total process CPU time in minutes.
 */

double cputime(void)
{
    struct tms buffer;

    if (times(&buffer) == -1)
        error("cputime in %s: times() call failed\n", getargv0());
    return ((buffer.tms_utime + buffer.tms_stime) / (60.0 * HZ));
}

/*
 * ERROR: scream and die quickly.
 */

void error(string fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);          /* invoke interface to printf       */
    fflush(stderr);                     /* drain std error buffer           */
    va_end(ap);
    exit(1);                            /* quit with error status           */
}

/*
 * EPRINTF: scream, but don't die yet.
 */

void eprintf(string fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);          /* invoke interface to printf       */
    fflush(stderr);                     /* drain std error buffer           */
    va_end(ap);
}

/*
 * SCANOPT: scan string of the form "word1,word2,..." for a match.
 * Words must be separated by commas only -- no spaces allowed!
 */

bool scanopt(string opt, string key)
{
    char *op, *kp;

    op = (char *) opt;                  /* start scan of option strings     */
    while (*op != NULL) {               /* loop while words left to check   */
        kp = key;                       /* (re)start scan of key word       */
        while ((*op != ',' ? *op : (char) NULL) == *kp) {
                                        /* char by char, compare word, key  */
            if (*kp++ == NULL)          /* reached end of key word, so...   */
                return (TRUE);          /* indicate success                 */
            op++;                       /* else go on to next char          */
        }
        while (*op != NULL && *op++ != ',')
                                        /* scan for start of next word      */
            continue;
    }
    return (FALSE);                     /* indicate failure                 */
}

/*
 * STROPEN: open a STDIO stream like fopen, with these extensions: (1)
 * existing files cannot be opend for writing unless mode == "w!"  or
 * mode == "a", (2) names of form "-" map to stdin/stdout, depending on
 * mode, and (3) names of the form "-num" up a stream to read/write file
 * descriptor num.
 */

stream stropen(string name, string mode)
{
    bool inflag;
    int fds;
    stream res;
    struct stat buf;

    inflag = streq(mode, "r");
    if (name[0] == '-') {
        if (streq(name, "-")) {
            fds = dup(fileno(inflag ? stdin : stdout));
            if (fds == -1)
                error("stropen in %s: cannot dup %s\n",
                      getargv0(), inflag ? "stdin" : "stdout");
        } else
            fds = atoi(&name[1]);
        res = fdopen(fds, streq(mode, "w!") ? "w" : mode);
        if (res == NULL)
            error("stropen in %s: cannot open f.d. %d for %s\n",
                  getargv0(), fds, inflag ? "input" : "output");
    } else {
        if (streq(mode, "w") && stat(name, &buf) == 0)
            error("stropen in %s: file \"%s\" already exists\n",
                  getargv0(), name);
        res = fopen(name, streq(mode, "w!") ? "w" : mode);
        if (res == NULL)
            error("stropen in %s: cannot open file \"%s\" for %s\n",
                  getargv0(), name, inflag ? "input" : "output");
    }
    return (res);
}
