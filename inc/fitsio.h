/*
 *  fitsio.h: header for fitsio.c
 *      sep-91      added DOUBLE                          PJT
 *      dec-91      added BYTE                            PJT
 *	mar-92      added some forgotten  declaration	  PJT
 *		    deleted all IMCALC references
 *   25-jun-92	repeaired nemo_stdinc_h -> _stdinc_h check  PJT
 *   21-feb-94  ansi + nemo -> now needs <stdinc.h>
 *    5-jan-98		added fitwrhd
 *   21-mar-00	offset/skip confusion cleared		    pjt
 *    7-aug-01  clarified ncards
 */

#ifndef FLOAT
    typedef float FLOAT;
#endif
#ifndef DOUBLE
    typedef double DOUBLE;
#endif

/* relic due to the fact that ANSI fortran cannot handle > 7 dims */
#ifndef MAXNAX
#  define MAXNAX 7
#endif

#ifndef private
#  define private static
#endif

#ifndef TRUE
#  define TRUE 1
#endif

#ifndef FALSE
#  define FALSE 0
#endif



#ifndef BSD

	/* 	Make sure if BSD is not set      */
#if defined(unicos) || defined(mc68k)
	/*	It`s really not needed 		 */
#else
#define BSD
#endif

#endif

#if !defined(_stdinc_h)		/* NEMO already defined this */
typedef unsigned char byte;
#define ABS(x)  ( (x)>0 ? (x) : -(x) )
#endif
/*	Some remaining fixed definitions and structure definitions */

#define TYPE_8INT   1
#define TYPE_16INT  2
#define TYPE_32INT  3
#define TYPE_64INT  4
#define TYPE_FLOAT  5
#define TYPE_DOUBLE 6

#define STATUS_OLD 1		/* reading a fits file */
#define STATUS_NEW 2		/* writing, but still in header */
#define STATUS_NEW_WRITE 3	/* writing, but now in data */

typedef struct { 
    int ncards;		/* is now 1-based !!! */
    int naxis;
    int axes[MAXNAX];
    int offset;		/* will/can change during I/O */
    int skip;		/* fixed */
    int type;
    int bytepix;
    int status;		/* STATUS  _OLD, _NEW, _NEW_WRITE   */
    stream fd;
    FLOAT bscale,bzero; 
} FITS;

/* public functions declared in fitsio.c */

FITS *fitopen (string, string, int, int *);
void fitclose (FITS *),
     fitsetpl (FITS *, int, int *),
     fitread  (FITS *, int, FLOAT *),
     fitwrite (FITS *, int, FLOAT *),
     fitrdhdr (FITS *, string, FLOAT *, FLOAT),
     fitrdhdi (FITS *, string, int *, int),
     fitrdhda (FITS *, string, string, string),
     fitwrhdr (FITS *, string, FLOAT),
     fitwrhdi (FITS *, string, int),
     fitwrhdl (FITS *, string, int),
     fitwrhda (FITS *, string, string),
     fitwrhd  (FITS *, string, string),
     fitwra   (FITS *, string, string);

void fit_setbitpix (int),
     fit_setscale  (FLOAT, FLOAT),
     fit_setblocksize (int);
int  fitexhd (FITS *, string);
