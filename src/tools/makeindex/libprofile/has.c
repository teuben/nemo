static char *trwsccs = "@(#)has.c	1.1 (TRW) 1/14/86";
#include <stdio.h>
#include "profile.h"

PROFILE_MARKER *profile_has_marker (s, m)
PROFILE_STANZA *s;
char *m;
{
	PROFILE_MARKER *x;
	int result;

	for (x = s->marker; x; x = x->next)
		if(glob_match(x->text, m) > 0)
			return(x);
	return((PROFILE_MARKER *)0);
}

/*
 * read down a linked list of stanzas looking
 * for a stanza that has the requested markers
 */
PROFILE_STANZA *profile_has_stanza(s, marker)
PROFILE_STANZA *s;
char *marker[];		/* terminated by a null pointer */
{
	int i;
	PROFILE_STANZA *x;

	if (s == NULL)
		return(s);
	x = s;
	do {
		for (i = 0; marker[i] != NULL; i++)
			if (profile_has_marker(x, marker[i]) == NULL)
				break;
		if (marker[i] == NULL)
			return(x);
		x = x->next;
	} while (x != s && x != NULL);

	return((PROFILE_STANZA *)NULL);
}

PROFILE_BINDING *profile_has_binding (s, b)
PROFILE_STANZA *s;
char *b;
{
	PROFILE_BINDING *x;

	for (x = s->binding; x; x = x->next)
		if (glob_match(x->name, b) > 0)
			return(x);
	return((PROFILE_BINDING *)0);
}
