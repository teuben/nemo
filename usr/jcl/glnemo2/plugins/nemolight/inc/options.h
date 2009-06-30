/*
 * Global options for NEMO programs
 *	This include file is automatically included by <stdinc.h>
 *	introduced	September 1989 - but not well advertised.....
 *
 *	27-mar-90	make it understand NEMO
 *	11-sep-90       DLOPEN
 *	16-jan-91	error, warning with line numbers..
 *	25-may-91	deleted some accidental trigraphs - added -Dsun3/4
 *	 6-apr-99       NEWIO
 *	 1-apr-01	LOADOBJ3 is the new .so loaded, instead of .o files
 *			this is for NEMO V3
 */

/* Possible CC machine options which can/have been used in #if defined(): */

/* TRIGGERS that can be used are:
    sun && sparc  sun4
    sun           sun3 (if sun present and sparc absent, a sun4 is assumed)
    mc68k         Unix PC (3b1 - the define 'u3b' doesn't seem to work)
    mips          various MIPS machines: IRIS workstation
    unicos        The Cray, with UNICOS of course
    vms           A VAX running VMS
    _trace_       The MultiFlow
    alliant       The Alliant 
    bsdvax        ultrix?
    NeXT	  NeXT boxes, note the 3 upper cases, and 1 lower case 'e'
 */
    
  
/*
 *   If in bodytrans() you want the functions to be saved, set SAVE_OBJ
 */
#define SAVE_OBJ

/*
 *   This is the new way of loading .so files instead of old .o files
 */
#define LOADOBJ3

/*
 *   If to create standalone subsystems...  (what that this mean PJT?)
 *   I think this is to signal programs that NEMO has been used, i.e.
 *   more fancy dprintf() can be used, nemoinp(), etc. ?
 *	== deprecated == ??
 */
#define NEMO_INC

/*
 *   If you have so-called standard dlopen() (loadobjDL.c)
 */
/*#define DLOPEN    */

/* 
 *  If you have NUMREC, you may also have to provide -I$NUMREC_DIR for cc
 */
#define NUMREC

/*
 *  If debugging of malloc() etc. is needed via the MNEM routines
 *  in $NEMO/src/kernel/core/mnem
 *
 *  DONT USE THIS - NOT DEBUGGED
 */
/* #define MNEMOSYNE */

#if defined(MNEMOSYNE)
# include <mnemosyne.h>
#endif

#if defined(sun)
#if defined(sparc)
#define sun4 1
#else
#define sun3 1
#endif
#endif


/* new experimental snapshot I/O - faster style */
/* careful, there appear to be some problems    */

/* #define NEWIO */
