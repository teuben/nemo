#define BIGENDIAN   1               /* make this 1 for left to right dump */
#define ASCII       1               /* make this 1 for ASCII machines */
#define MSDOS       0               /* make this 1 for MSDOS machines */
#define IBMPC_CHAR  0               /* make this 1 for IBM PC character set */
#define CBSIZE      16              /* bytes per line to dump */
#define FBSIZE (CBSIZE*3+1)         /* output array size */

#if IBMPC_CHAR
typedef unsigned char CHAR;         /* chars are unsigned for IBMPC */
#else
typedef char CHAR;                  /* they are signed on others */
#endif

/*
    I make my own version of isprint() here so I don't have to #include
    all the other stuff in <ctype.h>, like the tables for isascii(),
    isalpha(), islower(), etc.  Since I don't know of an easy way to make
    this macro work on non-ASCII machines, I just use <ctype.h> for them.

    I special-case IBMPC_CHAR since all characters above ' ' are printable
    on an IBM PC screen and printer.  If your machine or printer does not
    support this feature, set IBMPC_CHAR to 0.
*/

#if ASCII & ( IBMPC_CHAR == 0 )
#define isprint(c)  ( c >= ' ' && c <= '~' )
#endif

#if ASCII & IBMPC_CHAR
#define isprint(c)  ( c >= ' ' )
#endif

#if ASCII == 0
#include <ctype.h>
#endif
