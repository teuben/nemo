/* $Header */
/* Standard machine independent type definitions. */

#ifndef	STANDARD_TYPE	/* prevent multiple inclusions	*/
#define	STANDARD_TYPE

/*
 * Integers
 *     Tiny/UnsignedTiny	8+ bit integers
 *     Short/UnsignedShort	16+ bit integers
 *     Integer/UnsignedInteger	natural machine integer size
 *     Long/UnsignedLong	32+ bit integers
 *
 * Bits
 *     TinyBits		8+ bits
 *     Bits		16+ bits
 *     LongBits		32+ bits
 *
 * Booleans
 *     TinyBoolean
 *     Boolean
 *
 * Void
 *
 * Storage Classes
 *     Export		Seen in other compilation units
 *     Import		Supplied by another compilation unit
 *     Local		Unseen outside compilation unit
 */

/*
 * Each of the following sections for the integer types defines both
 * a base type and an extraction macro for the value.
 */

typedef char	Tiny;
/* Not all machines have signed characters so we may have to simulate them. */
#ifdef CHAR_IS_SIGNED
#   define TINY(x) (x)
#else
#   define TINY(x) (((x) & MIN_TINY) ? (~MAX_TINY | (x)) : (x))
#endif CHAR_IS_SIGNED

/* Not all compilers support unsigned chars so we may have to simulate them. */
#ifdef HAS_UNSIGNED_CHAR
    typedef unsigned char UnsignedTiny;
#else
    typedef char UnsignedTiny;
#endif HAS_UNSIGNED_CHAR
#ifdef CHAR_IS_SIGNED
#   define UNSIGNED_TINY(x) ((x) & MAX_UNSIGNED_TINY)
#else
#   define UNSIGNED_TINY(x) (x)
#endif

/*
 * All compilers have signed short integers.  This type is included
 * for lexical consistency.
 */
typedef short Short;

/* Not all compilers support unsigned shorts so we may have to simulate them. */
#ifdef HAS_UNSIGNED_SHORT
    typedef unsigned short UnsignedShort;
#else
    typedef short UnsignedShort;
#endif
#   define UNSIGNED_SHORT(x) ((unsigned)(x) & MAX_UNSIGNED_SHORT)

/* These types are solely for lexical consistency. */
typedef int Integer;
typedef	unsigned int UnsignedInteger;

typedef long Long;

/* Not all compilers support unsigned longs so we may have to simulate them. */
#ifdef HAS_UNSIGNED_LONG
    typedef unsigned long UnsignedLong;
#   define UNSIGNED_LONG(s) ((UnsignedLong)(x))
#else
    typedef long UnsignedLong;
#   define UNSIGNED_LONG(x) ((long)(x) & MAX_LONG)
#endif HAS_UNSIGNED_LONG

/* Boolean types take on only the values TRUE or FALSE. */
typedef	char TinyBoolean;
typedef	short Boolean;

/* This type is included for lexical consistency. */
typedef char Character;

/* Bit types are used only for bit set, clear and test operations. */
typedef	char TinyBits;
typedef	short Bits;
typedef	long LongBits;

/* Not all compilers support void functions so we may have to simulate it. */
#ifdef HAS_VOID
#   define Void void
#else
    typedef int Void;
#endif

/* Storage classes. */
#define	Export
#define Import extern
#define	Local static

#endif	STANDARD_TYPE
