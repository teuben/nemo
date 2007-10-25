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
 *               and this conflicts with our typedef
 *               also g++ seems to define it (2.6.3 at least)
 * 21-oct-93     defined ARGS to aid ANSI prototype def's    PJT
 *  5-nov-93	 added #ifndef/MAXLINES for ascii tables     PJT
 * 19-feb-94     added stdlib.h and C++/ANSIfied this header PJT
 *		note the PROTO kludge at the end for some bad acc bugs
 *  3-aug-94    while linux installation:
 *              - sym_angle, pos_angle macro's now in UPPER case -
 * 12-apr-95    cleanup for release 
 *              - got rid of all ARGS and PROTO; force prototypes
 *  5-jul-95    static version ID
 *  4-mar-96    test const argument to dprintf() for UNC's GCC 2.7.2
 * 13-feb-97    sizeof redefinition for allocate()
 *  2-apr-97    optional HUGE ???
 *  8-jan-98    temp. added common.c until laptop version merged
 *  6-apr-99    merged in AutoConf for 2.5 release
 * 22-jun-01    added some new ZENO macros for compatibility (-DNEMO)
 *  7-sep-01    added maxsizes.h
 *  8-dec-01    added errno.h
 *  7-may-02    removed some old starlab crap, added strneq
 * 13-mar-03    macros MIN,MAX deleted before redefinition
 * 12-jul-03    WD's macro name changed backfitted
 * 24-nov-03    include mathfns.h, which includes forced math.h
 * 22-jun-04    added some more useful astr constants
 * 25-jan-05    added powd/powi
 *  1-aug-05    nemo_string, nemo_stream for starlab interfaces
 * 29-nov-06    within() now using double's
 * 23-oct-07    added some more PI's and LN's from ZENO/Koda's mathfns.h
 */

#ifndef _stdinc_h      /* protect against re-entry */
#define _stdinc_h

#if !defined(__STDC__)
ERROR!  Sorry, NEMO now requires an ANSI C compiler
#endif

#define NEMO  1

/* If we're not using GNU C, elide __attribute__ */
#ifndef __GNUC__
#  define  __attribute__(x)  /*NOTHING*/
#endif

#include <maxsizes.h>   /* should be in $NEMOLIB - made during install */
#include <version.h>	/* our static version id - made during install */
#include <config.h>     /* should be in $NEMOLIB - made during install */
#include <options.h>    /* our private options   - manually edited still */

/* 	some backwards compat things */
#if defined(HAVE_LIBGSL)
#ifndef  HAVE_GSL
#define HAVE_GSL  1
#endif
#endif


/*
 * Always include stdio.h and stdlib.h and ....
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#if !defined(NOMATH)
#   include <math.h>
#endif

/*
 * our own math functions
 */

/*
 * string stuff, hopefully handled via the configure script (config.h)
 */

#if STDC_HEADERS
#include <string.h>
#else
# ifndef HAVE_STRCHR
#   define strchr  index
#   define strrchr rindex
# endif

  char *strchr(), *strrchr();

# ifndef HAVE_MEMCPY
#  define memcpy(d,s,n)  bcopy((s),(d),(n))
#  define memmove(d,s,n) bcopy((s),(d),(n))
# endif

#endif

/*
 * NULL: a long name for nothing....
 *	Since there has been a debate if NULL should be 0 or (void *)0,
 *	there is potential trouble. AIX is an example of trouble.
 *	According to ANSI either is OK. We recommend 0, though AIX
 *	defines another one which we can't use.
 */

#if defined(aix)
# undef NULL
#endif

#if !defined(NULL) 
#  define NULL 0
#endif

/*
 * BOOL, TRUE and FALSE: standard names for logical values.
 * The standard guarantees that the result of logical operations is
 * 1 if TRUE and 0 if FALSE.
 * (see  6.3.9, 6.3.10, 6.3.13, 6.3.14)
 *      18-nov-91:  disabled the #ifdef to always get 'bool' (AIX)
 *	24-feb-92:  bool as char ? - sendoff mesg to Josh
 *
 * WARNING: bool is now an official type in C++.....
 */

/*  probably should not use in C++ since bool is a new type in C++ */
#if !defined(__cplusplus)
#if !defined(CURSES_H) && !defined(_CURSES_H_)
    typedef unsigned char bool;
#endif
#endif

#if !defined(FALSE)
#  define FALSE (0)      
#  define TRUE  (1)
#endif


/*
 * BYTE: a short name for a handy chunk of bits. 
 *       (void *) could be fine too and is to be preferred in ANSI
 */

typedef unsigned char byte;

/*
 * STRING: for null-terminated strings which are not taken apart.
 * CONST_STRING:  and read-only strings (useful for C++)
 */

typedef char       *string;
typedef char       *nemo_string;
typedef const char *const_string;


/*
 * STREAM: a replacement for 'FILE *'.
 * ==> typedef struct _iobuf *stream        <== older version
 */
 
typedef FILE *stream;
typedef FILE *nemo_stream;

typedef struct _mstr {     /* see mstropen(3) */
  string tmplate;
  string filename;
  int status;
  int count;
  int mode;
  stream ostr;
} mstr;



/*
 * REAL, REALPTR: default type is double; if single precision calculation is
 * supported and desired, compile with -DSINGLEPREC. 
 *      DOUBLEPREC:     everything (variables & functions) is double.
 *      MIXEDPREC:      user values are float, -lm functions are double.
 *      SINGLEPREC:     everything (variables & functions) is float.
 * Note: The whole package must be compiled in one floating point
 *       type, although some routines are compiled in both single
 *       and double. Using this conmpile directive we can keep one
 *       version of the source code.
 * The "D" language defines a "real" as a floating point with the largest
 * number of bits the hardware supports (e.g. 80 bit on intel). Don't confuse
 * that with the real here, which is either float *or* double.
 */

/*   first off, by default, NEMO will be in DOUBLEPREC mode */

#if !defined(MIXEDPREC) && !defined(SINGLEPREC) && !defined(DOUBLEPREC)
#define DOUBLEPREC
#endif


#if defined(DOUBLEPREC)
#undef SINGLEPREC
#undef MIXEDPREC
typedef double real, *realptr;
#define Precision "DOUBLEPREC"
#endif

#if defined(MIXEDPREC)
#undef DOUBLEPREC
#undef SINGLEPREC
typedef float *realptr, real;
#define Precision "MIXEDPREC"
#endif

#if defined(SINGLEPREC)
#undef DOUBLEPREC
#undef MIXEDPREC
typedef float real, *realptr;
#define Precision "SINGLEPREC"
#endif

/* NOTE: to be deprecated
 * The following conveniences cannot be used in full ANSI C or C++:
 * one needs a new type for each kind of function pointer, i.e.
 * not only dependant on the return type, but also on the arguments.
 * Note that C++ does not allow a varargs-only argument list, C does.
 * 
 * PROC: pointer to a "function" object which returns no value.
 * BPROC: pointer to a boolean-valued function
 * IPROC: pointer to an integer-valued function.
 * RPROC: pointer to a real-valued function.
 */

typedef void (*proc)();
typedef bool (*bproc)();	/* problems on sun3 ??? aix too */
typedef int (*iproc)();
typedef real (*rproc)();

typedef int  (*qsort_proc)(const void *, const void *);
typedef real (*real_proc)(real);
typedef real (*double_proc)(double);
typedef real (*float_proc)(float);

/*
 * LOCAL: declare something to be local to a file.
 * PERMANENT: declare something to be permanent data within a function.
 */

#define local     static
#define permanent static

#if 0
#if !defined(__STDC__)
#   define const
#   define volatile
#endif
#endif

/*
 * STREQ, STRNEQ, STRNULL: handy string-equality macros
 */

#define streq(x,y)    (strcmp((x), (y)) == 0)
#define strneq(x,y,z) (strncmp((x), (y), (z)) == 0)
#define strnull(x)    (strcmp((x), "") == 0)

/*
 *  PI, etc.  --  mathematical constants
 *  Note: some <math.h> already define PI, the official name is M_PI and M_PIl
 *  See also: http://3.141592653589793238462643383279502884197169399375105820974944592.com
 */

#ifndef PI
#define   PI         3.141592653589793238462643
#endif
#define   TWO_PI     6.283185307179586476925287
#define   FOUR_PI   12.56637061435917295385
#define   HALF_PI    1.57079632679489661923
#define   FRTHRD_PI  4.18879020478639098462
#define   INV_PI     0.31830988618379067154
#define   TENSVN_IPI 0.45472840883398667363
#define   ONETRD     0.33333333333333333333
#define   TWOTRD     0.66666666666666666667
#define   FORTRD     1.33333333333333333333
#define   ONESVN     0.14285714285714285714
#define   ONESIX     0.16666666666666666667

/*
 * angular conversion factors (multiplicative)
 *   DR2H   radians to hours
 *   DH2R   hours   to radians
 *   DR2D   radians to degrees
 *   DD2R   degrees to radians
 *   DR2AS  radians to arcsec
 *   DAS2R  arcsec  to radians
 *
 */
#define DR2H    3.819718634205488058453210
#define DH2R    0.2617993877991494365385536
#define DR2D    57.29577951308232087679815
#define DD2R    1.745329251994329576923691e-2
#define DR2AS   206264.8062470963551564733
#define DAS2R   4.848136811095359935899141e-6

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
 *  3-aug-94: changed lower case to upper case for consistency
 */

#define  POS_ANGLE(phi)    ((phi) - TWO_PI * floor((phi)/TWO_PI ))
#define  SYM_ANGLE(phi)    ((phi) - TWO_PI * floor(((phi)+PI)/TWO_PI ))

/*
 * Some "useful" macro's:    this is where C just falls down,
 *		             everybody pretty much comes up with these names....
 *
 *  ABS: returns the absolute value of its argument
 *  SGN: returns the sign of its argument
 *  MIN: returns the argument with the lowest value
 *  MAX: returns the argument with the highest value
 *  RND: roundup a number (now: ROUNDUP)
 *  RNG: range a number between a low and high value (now: RANGE)
 *
 *  Note: be aware of side effects here
 */

#ifdef ABS
#  undef ABS
#endif
#define   ABS(x)       (((x) < 0) ? -(x) : (x))

#ifdef MIN	/* WD 13/03/03: added to prevent clash with other definitions */
#  undef MIN
#endif
#define   MIN(x,y)     (((x) < (y)) ? (x) : (y))

#ifdef MAX      /* WD 13/03/03: added to prevent clash with other definitions */
#  undef MAX
#endif
#define   MAX(x,y)     (((x) > (y)) ? (x) : (y))

#define   SGN(x)       (((x) < 0) ? (-1) : ((x) > 0) ? 1 : 0)
#if 0
/* nrutil.h defines SIGN(a,b) */
#define   SIGN(x)      (((x) < 0) ? (-1) : ((x) > 0) ? 1 : 0)
#endif
#define   RND(x,y)     ((y)*(((x)+(y)-1)/(y)))
#define   ROUNDUP(x,y) ((y)*(((x)+(y)-1)/(y)))
#define   RNG(x,l,h)   (((x) > (h)) ? (h) : (((x) < (l)) ? (l) : (x)))
#define   RANGE(x,l,h) (((x) > (h)) ? (h) : (((x) < (l)) ? (l) : (x)))


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

/* SYSV unicos has no random, srandom -- see xrandom.c for better solutions */
#if defined(unicos)
#   define random rand
#   define srandom srand
#endif

/* We have to decide which memory copy/compare function we're supporting */
/* The BSTRING(3) or MEMORY(3) functions */
/* Posix says:  MEMORY : see string.h */
/* let's keep it out now, see above with string handling */
#if 0
#if defined(mc68k) || defined(SYSV)
# define bcopy(from,to,n) memcpy(to,from,n)
#endif
#endif

/*
 * Prototypes for some of NEMO's frequently used kernel functions
 * should really use some more descriptive .h files
 */

#if defined(__cplusplus)
extern "C" {
#endif

/* various math */

#include <mathfns.h>

/* io/stropen.c */

extern stream stropen(const_string, string);
extern int    strdelete(stream, bool);
extern string strname(stream);
extern bool   strseek(stream);

/* io/mstropen.c */
extern mstr  *mstr_init(string tmplate);
extern stream mstr_open(mstr *mp, string mode);
extern void   mstr_close(mstr *mp);
extern int    mstr_count(mstr *mp);
extern int    mstr_multi(mstr *mp);

/* io/endian.c */
extern bool BigEndian(void);
extern bool LittleEndian(void);

/* io/filesecret.c */
extern void   strclose(stream);


/* error.c dprintf.c */
/* C99 stdargs example of macro usage:   #define HELLO(a,...)  error(a,__VA_ARGS__)   */
/* GNU stdargs (deprecated now)          #define HELLO(a,args...)  error(a,##args)    */
void error(string, ...);
void errorn(string, ...);
void warning(string, ...);
int progress(double dtime, string, ...);
int nemo_dprintf(int, const_string, ...);    /* NEMO has same name as libc */
bool nemo_debug(int);

#ifndef HAVE_DPRINTF
#define dprintf		nemo_dprintf
#else
#define dprintf		nemo_dprintf
#endif

/* eprintf is ZENO's "warning" */
#define eprintf warning

/* core/allocate.c */
extern void *allocate(size_t);
extern void *reallocate(void *, size_t);
/* this is to shut up e.g. gcc when allocate(n*sizeof(T)) is used */
#define  sizeof  (size_t)sizeof

/* core/common.c */
void set_common(int id, int byte_size);
int get_common(int id, int elt_size, int bucket_size);
byte *open_common(int id);
void close_common(int id);

extern bool scanopt(string, string);


/* core/cputime.c */
extern double cputime(void);

  /* misc/sqr.c */
extern double sqr(double);
extern double qbe(double);
extern double dex(double);

  /* misc/pow.c */
extern double powi(double,int);
extern double powd(double,double);

  /* misc/log2.c */
  /* cygwin defines log2 as a macro */
#if defined(log2)
#define HAVE_LOG2
#endif
#if !defined(HAVE_LOG2)
extern double log2(double);
#endif

  /* cores/bswap.c */
extern void bswap(void *vdat, int len, int cnt);
extern void bswap_litend(void *vdat, int len, int cnt);
extern void bswap_bigend(void *vdat, int len, int cnt);
#define bswapr(p,cnt)  bswap(p,sizeof(real),cnt)
#define bswapd(p,cnt)  bswap(p,sizeof(double),cnt)
#define bswapf(p,cnt)  bswap(p,sizeof(float),cnt)
#define bswapl(p,cnt)  bswap(p,sizeof(long),cnt)
#define bswapi(p,cnt)  bswap(p,sizeof(int),cnt)
#define bswaps(p,cnt)  bswap(p,sizeof(short),cnt)

  /* misc/within.c */
extern bool within(double val, string range, double fuzz);


#if defined(__cplusplus)
}
#endif

/* for tables: (ought to go into maxsizes.h) */

#ifndef MAXLINES
#define MAXLINES 10000
#endif

/* for solaris compiler SC4.2, doesn't appear to know HUGE .... */
/* should fit double and single precision floating point        */
#ifndef HUGE
#define HUGE 1e30
#endif

#endif /* _stdinc_h */
