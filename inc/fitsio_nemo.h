/*
 *  fitsio_nemo.h: header for NEMO's fitsio.c, with optional wrapper for CFITSIO
 *                 (this file used to be called fitsio.h, but CFITSIO is using it also)
 *    	
 *      sep-91      added DOUBLE                          PJT
 *      dec-91      added BYTE                            PJT
 *	mar-92      added some forgotten  declaration	  PJT
 *		    deleted all IMCALC references
 *   25-jun-92	repeaired nemo_stdinc_h -> _stdinc_h check  PJT
 *   21-feb-94  ansi + nemo -> now needs <stdinc.h>
 *    5-jan-98		added fitwrhd
 *   21-mar-00	offset/skip confusion cleared		    pjt
 *    7-aug-01  clarified ncards
 *   18-dec-01  renamed this file from fitsio.h to fitsio_nemo.h
 *              and added optional CFITSIO wrapper stuff
 *   23-jul-02  add fitresize
 */

#ifndef _fitsio_nemo_h
#define _fitsio_nemo_h

/*  NOTE:   this also assumed cfitsio was installed in its own directory/namespace ... */
#ifdef HAVE_LIBCFITSIO
#include <cfitsio/fitsio.h>
#endif

#ifndef FLOAT
    typedef float FLOAT;
#endif
#ifndef DOUBLE
    typedef double DOUBLE;
#endif

/* relic due to the fact that ANSI fortran cannot handle > 7 dims */
/* CFITSIO often uses 9                                           */
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

#ifdef SUPPORT_64BIT_INTEGERS
/* this is from cfitsio */
#endif


#define STATUS_OLD 1		/* reading a fits file */
#define STATUS_NEW 2		/* writing, but still in header */
#define STATUS_NEW_WRITE 3	/* writing, but now in data */

#ifdef HAVE_LIBCFITSIO

typedef fitsfile FITS;

#else

typedef struct { 
    int ncards;		/* is now 1-based !!! */
    int naxis;
    int axes[MAXNAX];
    size_t offset;		/* will/can change during I/O */
    size_t skip;		/* fixed */
    int type;
    int bytepix;
    int status;		/* STATUS  _OLD, _NEW, _NEW_WRITE   */
    stream fd;
    FLOAT bscale,bzero; 
} FITS;

#endif

/* public functions declared in fitsio.c */

FITS *fitopen (string, string, int, int *);
void fitclose (FITS *),
     fitresize(FITS *, int, int *),
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
#endif
