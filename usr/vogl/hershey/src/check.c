#include <stdio.h>

extern int	hLoaded;

/*
 * check_loaded
 *
 * 	Checks and prints out a message if the font isn't loaded.
 */
void
check_loaded(who)
	char	*who;
{
	if (!hLoaded) {
		fprintf(stderr, "%s: no hershey font loaded.\n", who);
		exit(1);
	}
}
