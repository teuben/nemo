/*
 * STDINC.H: standard include file for C programs.
 * Copyright (c) 1986  Joshua Edward Barnes  Princeton, NJ.
 *
 * xx-jul-89     stdinc_h added, 'stream' defined as *FILE   PJT
 * 10-sep-90     void  not always needed to be #defined	     PJT
 * 16-sep-90     merging header files Nemo with Starlab	     PJT/PH
 * 19-oct-90	 bcopy/memcpy
 *  1-mar-91     added strname()			     PJT
 * 25-may-91     deleted accidental trigraph                 PJT
 * 18-nov-91	 bool never got defined on AIX               PJT
 * 22-nov-91     special #undef/#def for NULL on AIX	     PJT
 * 24-feb-92     bool is defined in <curses.h> as char ...   PJT
 *               and this causes severe headaches.....
 */

                    /* protection against re-entry */
#ifndef nemo_stdinc_h
#define nemo_stdinc_h

      /* various operating system dependant thingo's can go in options.h */
#include <options.h>

/*
 * If not already loaded, include <stdio.h> 
 *     (is use of FILE really POSIX OK?)
 */

#if !defined(FILE)
#  include <stdio.h>
#endif

#if !defined(MACHINE)
#   define MACHINE "UNDEFINED"
#endif

#if !defined(NOMATH)
#   include <math.h>
#endif

/*
 * NULL: a long name for nothing.
 */

#if defined(aix)
# undef NULL
#endif

#if !defined(NULL) 
#  define NULL 0
#endif

/*
 * BOOL, TRUE and FALSE: standard names for logical values.
 *      18-nov-91:  disabled the #ifdef to always get 'bool' (AIX)
 *	24-feb-92:  bool as char ? - sendoff mesg to Josh
 */

/*  #if !defined(TRUE)  */

#if 0
   typedef short int bool;      /* do not use char for logicals */
#else
#if !defined(CURSES_H)
   typedef char bool;              /* see also curses.h .... */
#endif
#endif
#if !defined(FALSE)
#  define FALSE (0)
#  define TRUE  (1)
#endif

/*  #endif  */

/*
 * VOID: type specifier used to declare procedures called for
 * side-effect only.  Note: this slightly kinky substitution
 * is used to so that one need not declare something to be
 * void before using it. Since ANSI is around, we should not need
 * to define void anymore... some old machine may still need it
 */

#ifdef mc68k
#define void int
#endif

/*
 * BYTE: a short name for a handy chunk of bits.
 */

typedef unsigned char byte;

/*
 * STRING: for null-terminated strings which are not taken apart.
 */

typedef char *string;

/*
 * STREAM: a replacement for FILE *.
 */
#if 0
typedef struct _iobuf *stream;        /* assume stdio.h included 1st */
#endif
typedef FILE *stream;		        /* more general *? */

stream stropen(/* string name, string mode */);
string strname(/* stream str */); 
void   strclose(/* stream str */);
int    strdelete(/* stream str, bool scratch */);
bool   strseek(/* stream str */);

/*
 * REAL: default type is double; if single precision calculation is
 * supported and desired, compile with -DSINGLEPREC. But really the
 * whole of the package should be recompiled to be sure.
 */

#if !defined(SINGLEPREC)
  typedef  double  real, *realptr;
#else
  typedef  float   real, *realptr;
#endif

/*
 * PROC: pointer to a "function" object which returns no value.
 * BPROC: pointer to a boolean-valued function
 * IPROC: pointer to an integer-valued function.
 * RPROC: pointer to a real-valued function.
 */

typedef void (*proc)();
/*typedef bool (*bproc)();*/	/* problems on sun3 ??? aix too */
typedef int (*iproc)();
typedef real (*rproc)();

/*
 * LOCAL: declare something to be local to a file.
 * PERMANENT: declare something to be permanent data within a function.
 */

#define local     static
#define permanent static

#if !defined(__STDC__)
#   define const
#   define volatile
#endif

/*
 * STREQ: handy string-equality macro.
 */

#define streq(x,y) (strcmp((x), (y)) == 0)

/*
 *  PI, etc.  --  mathematical constants
 */

#define   PI         3.14159265358979323846

#define   TWO_PI     6.28318530717958647693
#define   FOUR_PI   12.56637061435917295385
#define   HALF_PI    1.57079632679489661923
#define   FRTHRD_PI  4.18879020478639098462

/*
 *  POS_ANGLE, SYM_ANGLE: bring angular variables into standard form.
 *    To map an angular variable 'phi' into the range [0, 2pi),
 *    keeping 'phi' invariant modulo 2pi, use the statement
 *            phi = pos_angle(phi);
 *    To map an angular variable 'phi' into the range [-pi, pi),
 *    keeping 'phi' invariant modulo 2pi, use the statement
 *            phi = sym_angle(phi);
 *
 *  Note: be aware of side effects here
 */

#define  pos_angle(phi)    ((phi) - TWO_PI * floor((phi)/TWO_PI ))
#define  sym_angle(phi)    ((phi) - TWO_PI * floor(((phi)+PI)/TWO_PI ))

/*
 *  ABS: returns the absolute value of its argument
 *  SGN: returns the sign of its argument
 *  MAX: returns the argument with the highest value
 *  MIN: returns the argument with the lowest value
 *
 *  Note: be aware of side effects here
 */

#define   ABS(x)       (((x) < 0) ? -(x) : (x))
#define   SGN(x)       (((x) < 0) ? (-1) : ((x) > 0) ? 1 : 0)
#define   MAX(x,y)     (((x) > (y)) ? (x) : (y))
#define   MIN(x,y)     (((x) < (y)) ? (x) : (y))

/*
 *      On some older BSD (4.2?) implementations the strchr() and strrchr()
 *      may not be implemented, and are known under the name 
 *      index() and rindex()
 *	Also on these machines <strings.h> may not exist, and Nemo should
 *	define "strings.h" to #include <string.h>
 *	Here we adhere to the ANSI standard and leave out this stuff
 */

#if 0
#   define index strchr
#   define rindex strrchr
#endif

/* unicos has no random, srandom ?? */
#if defined(unicos)
#   define random rand
#   define srandom srand
#endif

/* We have to decide which memory copy/compare function we're supporting */
/* The BSTRING(3) or MEMORY(3) functions */
/* Posix says:  MEMORY ? */
#ifdef mc68k
# define bcopy(a,b,c) memcpy(b,a,c)
#endif


/* Declare some often used functions, that gcc would else complain
 * about - at some point we need to ANSI-fy things, and this will go
 * into more formal header files.. PJT 25-feb-92
 */

void error(/* string ... */);
void warning(/* string ... */);
void dprintf(/* int, string ... */);

#if 0
int printf(/* string ... */);
int fprintf(/* FILE*, string ... */);
#endif

#endif /* _nemo_stdinc */
