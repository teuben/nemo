static char *trwsccs = "@(#)shift.c	1.1 (TRW) 1/14/86";
#include <ctype.h>

/*
 * Downshifts a string in place.
 */
char *string_downshift(s)
char *s;
{

	register char *x = s;
	for (; *x; x++)
		if (isupper(*x))
			*x = tolower(*x);
	return(s);
}

/*
 * Upshifts a string in place.
 */
char *string_upshift(s)
char *s;
{

	register char *x = s;
	for (; *x; x++)
		if (islower(*x))
			*x = toupper(*x);
	return(s);
}
