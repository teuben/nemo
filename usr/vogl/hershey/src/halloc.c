
#include <stdio.h>

extern	char	*malloc();

/*
 * hallocate
 *
 *	Allocate some memory, barfing if malloc returns NULL.
 */
char *
hallocate(size)
	unsigned	size;
{
	char	*p;

	if ((p = (char *)malloc(size)) == (char *)NULL) {
		fprintf(stderr,"hallocate: request for %d bytes returned NULL", size);
		exit(1);
	}

	return (p);
}
