/* $Header */
/*
 * This file defines the machine/compiler C environment. It defines
 * pre-processor macros that tell what C features are supported.
 *
 * #define HAS_UNSIGNED_SHORT	Implies unsigned shorts are supported
 * #define CHAR_IS_SIGNED	Implies chars are signed
 * #define HAS_UNSIGNED_CHAR	Implies unsigned chars are supported
 * #define HAS_UNSIGNED_LONG	Implies unsigned longs are supported
 * #define BITS_PER_CHAR n	Number of bits in a char
 * #define BITS_PER_INT n	Number of bits in an int
 * #define BITS_PER_LONG n	Number of bits in a long
 * #define BITS_PER_POINTER n	Number of bits in a pointer
 * #define BITS_PER_SHORT n	Number of bits in a short
 * #define HAS_VOID		Implies void function type is supported
 */

#ifndef	STANDARD_ENVIRON	/* prevent multiple inclusions	*/

#if defined(vax) || defined(pyr) || defined(sun) || defined(__STDC__)
#   define HAS_UNSIGNED_SHORT
#   define CHAR_IS_SIGNED
#   define HAS_UNSIGNED_CHAR
#   define HAS_UNSIGNED_LONG
#   define HAS_VOID

#   define BITS_PER_CHAR 8
#   define BITS_PER_SHORT 16
#   define BITS_PER_INT 32
#   define BITS_PER_LONG 32
#   define BITS_PER_POINTER 32
#   define STANDARD_ENVIRON
#endif

#endif STANDARD_ENVIRON

/* make sure a known processor type was	specified */
#ifndef	STANDARD_ENVIRON
#   include "Processor type unknown or unspecified"
#endif STANDARD_ENVIRON
