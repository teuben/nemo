/*	machdep.h	1.11	90/05/22	*/

/*
 * Copyright (c) 1988, 1990 by Sam Leffler.
 * All rights reserved.
 *
 * This file is provided for unrestricted use provided that this
 * legend is included on all tape media and as a part of the
 * software program in whole or part.  Users may copy, modify or
 * distribute this file at will.
 */

#ifndef _MACHDEP_
#define	_MACHDEP_
/*
 * Machine dependent definitions:
 *   o floating point formats
 *   o byte ordering
 *   o signed/unsigned chars
 *
 * NB, there are lots of assumptions here:
 *   - 32-bit natural integers	(sign extension code)
 *   - native float is 4 bytes	(floating point conversion)
 *   - native double is 8 bytes	(floating point conversion)
 */

#if defined(sun) || defined(sparc) || defined(stellar) || defined(MIPSEB) || defined(hpux) || defined(apollo) || defined(NeXT)
#define	BIGENDIAN	1
#define	IEEEFP
#endif /* sun || sparc || stellar || MIPSEB || hpux || apollo || NeXT */

/* MIPSEL = MIPS w/ Little Endian byte ordering (e.g. DEC 3100) */
#if defined(MIPSEL)
#define	BIGENDIAN	0
#define	IEEEFP
#endif /* MIPSEL */

#ifdef transputer
#define	BIGENDIAN	0
#define	UCHARS			/* yech! unsigned characters */
#define	IEEEFP			/* IEEE floating point supported */
#endif /* transputer */

#if defined(m386) || defined(M_I86) || defined(i386)
#define BIGENDIAN	0
#define IEEEFP			/* IEEE floating point supported */
#endif /* m386 || M_I86 || i386 */

#ifdef mips
#define	UCHARS			/* unsigned characters */
#endif /* mips */

#ifdef IEEEFP
typedef	struct ieeedouble nativedouble;
typedef	struct ieeefloat nativefloat;
#define	ISFRACTION(e)	(1022 - 4 <= (e) && (e) <= 1022 + 15)
#define	EXTRACTFRACTION(dp, fract) \
   ((fract) = ((unsigned long)(1<<31)|((dp)->native.mant<<11)|\
      ((dp)->native.mant2>>21)) >> (1022+16-(dp)->native.exp))
#define	EXTRACTEXPONENT(dp, exponent) ((exponent) = (dp)->native.exp)
#define	NATIVE2IEEEFLOAT(fp)
#define	IEEEFLOAT2NATIVE(fp)
#define	IEEEDOUBLE2NATIVE(dp)

#define	TIFFSwabArrayOfFloat(fp,n)  TIFFSwabArrayOfLong((unsigned long *)fp,n)
#define	TIFFSwabArrayOfDouble(dp,n) TIFFSwabArrayOfLong((unsigned long *)dp,2*n)
#endif /* IEEEFP */

#ifdef tahoe
#define	BIGENDIAN	1

typedef	struct {
	unsigned	sign:1;
	unsigned	exp:8;
	unsigned	mant:23;
	unsigned	mant2;
} nativedouble;
typedef	struct {
	unsigned	sign:1;
	unsigned	exp:8;
	unsigned	mant:23;
} nativefloat;
#define	ISFRACTION(e)	(128 - 4 <= (e) && (e) <= 128 + 15)
#define	EXTRACTFRACTION(dp, fract) \
    ((fract) = ((1<<31)|((dp)->native.mant<<8)|((dp)->native.mant2>>15)) >> \
	(128+16-(dp)->native.exp))
#define	EXTRACTEXPONENT(dp, exponent) ((exponent) = (dp)->native.exp - 2)
/*
 * Beware, over/under-flow in conversions will
 * result in garbage values -- handling it would
 * require a subroutine call or lots more code.
 */
#define	NATIVE2IEEEFLOAT(fp) { \
    if ((fp)->native.exp) \
        (fp)->ieee.exp = (fp)->native.exp - 129 + 127;	/* alter bias */\
}
#define	IEEEFLOAT2NATIVE(fp) { \
    if ((fp)->ieee.exp) \
        (fp)->native.exp = (fp)->ieee.exp - 127 + 129; 	/* alter bias */\
}
#define	IEEEDOUBLE2NATIVE(dp) { \
    if ((dp)->native.exp = (dp)->ieee.exp) \
	(dp)->native.exp += -1023 + 129; \
    (dp)->native.mant = ((dp)->ieee.mant<<3)|((dp)->native.mant2>>29); \
    (dp)->native.mant2 <<= 3; \
}
/* the following is to work around a compiler bug... */
#define	SIGNEXTEND(a,b)	{ char ch; ch = (a); (b) = ch; }

#define	TIFFSwabArrayOfFloat(fp,n)  TIFFSwabArrayOfLong((unsigned long *)fp,n)
#define	TIFFSwabArrayOfDouble(dp,n) TIFFSwabArrayOfLong((unsigned long *)dp,2*n)
#endif /* tahoe */

#ifdef vax
#define	BIGENDIAN	0

typedef	struct {
	unsigned	mant1:7;
	unsigned	exp:8;
	unsigned	sign:1;
	unsigned	mant2:16;
	unsigned	mant3;
} nativedouble;
typedef	struct {
	unsigned	mant1:7;
	unsigned	exp:8;
	unsigned	sign:1;
	unsigned	mant2:16;
} nativefloat;
#define	ISFRACTION(e)	(128 - 4 <= (e) && (e) <= 128 + 15)
#define	EXTRACTFRACTION(dp, fract) \
    ((fract) = ((1<<31)|((dp)->native.mant1<<16)|(dp)->native.mant2)) >> \
	(128+16-(dp)->native.exp))
#define	EXTRACTEXPONENT(dp, exponent) ((exponent) = (dp)->native.exp - 2)
/*
 * Beware, these do not handle over/under-flow
 * during conversion from ieee to native format.
 */
#define	NATIVE2IEEEFLOAT(fp) { \
    float_t t; \
    if (t.ieee.exp = (fp)->native.exp) \
	t.ieee.exp += -129 + 127; \
    t.ieee.sign = (fp)->native.sign; \
    t.ieee.mant = ((fp)->native.mant1<<16)|(fp)->native.mant2; \
    *(fp) = t; \
}
#define	IEEEFLOAT2NATIVE(fp) { \
    float_t t; int v = (fp)->ieee.exp; \
    if (v) v += -127 + 129;		/* alter bias of exponent */\
    t.native.exp = v;			/* implicit truncation of exponent */\
    t.native.sign = (fp)->ieee.sign; \
    v = (fp)->ieee.mant; \
    t.native.mant1 = v >> 16; \
    t.native.mant2 = v;\
    *(fp) = v; \
}
#define	IEEEDOUBLE2NATIVE(dp) { \
    double_t t; int v = (dp)->ieee.exp; \
    if (v) v += -1023 + 129; 		/* if can alter bias of exponent */\
    t.native.exp = v;			/* implicit truncation of exponent */\
    v = (dp)->ieee.mant; \
    t.native.sign = (dp)->ieee.sign; \
    t.native.mant1 = v >> 16; \
    t.native.mant2 = v;\
    t.native.mant3 = (dp)->mant2; \
    *(dp) = t; \
}

#define	TIFFSwabArrayOfFloat(fp,n)  TIFFSwabArrayOfLong((unsigned long *)fp,n)
#define	TIFFSwabArrayOfDouble(dp,n) TIFFSwabArrayOfLong((unsigned long *)dp,2*n)
#endif /* vax */

/*
 * These unions are used during floating point
 * conversions.  The macros given above define
 * the conversion operations.
 */

typedef	struct ieeedouble {
#if BIGENDIAN == 1
	unsigned	sign:1;
	unsigned	exp:11;
	unsigned long	mant:20;
	unsigned	mant2;
#else
#if !defined(vax)
#ifdef INT_16_BIT	/* MSDOS C compilers */
	unsigned long	mant2;
	unsigned long	mant:20;
	unsigned long	exp:11;
	unsigned long	sign:1;
#else /* 32 bit ints */
	unsigned	mant2;
	unsigned long	mant:20;
	unsigned	exp:11;
	unsigned	sign:1;
#endif /* 32 bit ints */
#else
	unsigned long	mant:20;
	unsigned	exp:11;
	unsigned	sign:1;
	unsigned	mant2;
#endif /* !vax */
#endif
} ieeedouble;
typedef	struct ieeefloat {
#if BIGENDIAN == 1
	unsigned	sign:1;
	unsigned	exp:8;
	unsigned long	mant:23;
#else
#ifdef INT_16_BIT	/* MSDOS C compilers */
	unsigned long	mant:23;
	unsigned long	exp:8;
	unsigned long	sign:1;
#else /* 32 bit ints */
	unsigned long	mant:23;
	unsigned	exp:8;
	unsigned	sign:1;
#endif /* 32 bit ints */
#endif
} ieeefloat;

typedef	union {
	ieeedouble	ieee;
	nativedouble	native;
	char		b[8];
	double		d;
} double_t;

typedef	union {
	ieeefloat	ieee;
	nativefloat	native;
	char		b[4];
	float		f;
} float_t;

/*
 * Sign extend the low order byte of one
 * integer to form a 32-bit value that's
 * placed in the second.
 */
#if defined(UCHARS) && !defined(SIGNEXTEND)
/* there may be a better way on some machines... */
#define	SIGNEXTEND(a,b)	{ if (((b) = (a)) & 0x80) (b) |= 0xffffff00; }
#else
#if !defined(SIGNEXTEND)
#define	SIGNEXTEND(a,b)	((b) = (char)(a))
#endif
#endif /* UCHARS */
#endif /* _MACHDEP_ */
