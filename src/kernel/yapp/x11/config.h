#include <stdinc.h>
#include <ctype.h>
#include <math.h>
#include <yapp.h>

#ifdef POSIX
#include <unistd.h>
#endif

#if 0
	/* don't use these anymore - stdinc.h and/or ANSI headers
	 * are supposed to pick them up
         */
#if defined(USG) || defined(STDC_HEADERS)
#include <string.h>
#define index strchr
#define rindex strrchr
#define bcopy(from, to, len) memcpy ((to), (from), (len))
#define bzero(s, n) memset ((s), 0, (n))
#else /* USG or STDC_HEADERS */
#include <strings.h>
extern char *malloc();
extern char *realloc();
extern double atof();
#endif /* USG or STDC_HEADERS */

#endif
