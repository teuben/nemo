/*
** tr2latex - troff to LaTeX converter
** $Id$
** COPYRIGHT (C) 1987 Kamal Al-Yahya, 1991,1992 Christian Engel
**
** Module: setups.h
**
** setup file
*/

#ifdef __TURBOC__
#  ifndef __COMPACT__
#    error "The COMPACT model is needed"
#  endif
#endif

#include	<stdio.h>
#include	<ctype.h>
#include	<errno.h>
#include	<stdlib.h>

#include	<time.h>


#if defined(__TURBOC__) | defined(MSC)
#  include	<io.h>          /* for type declarations */
#endif

#ifdef UCB
#  include	<strings.h>
#else
#  include	<string.h>
#endif

#ifdef VMS
#  define GOOD 1
#else
#  define GOOD 0
#endif

#ifdef __TURBOC__
/*--- ONE MEGABYTE for EACH BUFFER? Really? - not for a brain-dead cpu
      like the 8086 or 80286. Would my account on the VAX let me have
      this much?
      Maybe, one day when we all have an 80386+, or a SPARC, or a real
      cpu....                                                           ---*/
#define	MAXLEN	(0xfff0)		/* maximum length of document */
								/* A segment is 64k bytes! */
#else
/* ... but no problem on virtual operating systems */
#define	MAXLEN	(1024*1024)		/* maximum length of document */
#endif

#define	MAXWORD	250				/* maximum word length */
#define	MAXLINE	500				/* maximum line length */
#define	MAXDEF	200				/* maximum number of defines */
#define MAXARGS 128				/* maximal number of arguments */

#define EOS		'\0'			/* end of string */

#if defined(__TURBOC__) | defined (VAXC)
#  define index strchr
#endif

#ifdef MAIN				/* can only declare globals once */
#  define GLOBAL
#else
#  define GLOBAL extern
#endif

GLOBAL int math_mode,	/* math mode status */
	       de_arg,		/* .de argument */
	       IP_stat,		/* IP status */
	       QP_stat,		/* QP status */
           TP_stat;		/* TP status */

#ifdef DEBUG
GLOBAL int debug_o;
GLOBAL int debug_v;
#endif

GLOBAL struct defs {
	char *def_macro;
	char *replace;
	int illegal;
} def[MAXDEF];

GLOBAL struct mydefs {
	char *def_macro;
	char *replace;
	int illegal;
	int arg_no;
	int par;		/* if it impiles (or contains) a par break */
} mydef[MAXDEF];

GLOBAL struct measure {
	char old_units[MAXWORD];	float old_value;
	char units[MAXWORD];		float value;
	char def_units[MAXWORD];	/* default units */
	int def_value;			/* default value: 0 means take last one */
} linespacing, indent, tmpind, space, vspace;

typedef unsigned char bool;

extern int errno;
