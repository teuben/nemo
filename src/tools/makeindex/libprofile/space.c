static char *trwsccs = "@(#)space.c	1.1 (TRW) 1/14/86";
#include "profile.h"

extern char *calloc();

PROFILE_STANZA *profile_stanza_space ()
{
	return((PROFILE_STANZA *)calloc(1, sizeof(PROFILE_STANZA)));
}

PROFILE_MARKER *profile_marker_space (n)
int n;
{
	char *space;
	PROFILE_MARKER *m = (PROFILE_MARKER *)0;

	if (space = calloc(1, sizeof(PROFILE_MARKER) + n + 1)) {
		m = (PROFILE_MARKER *)space;
		m->text = space + sizeof(PROFILE_MARKER);
	}
	return(m);
}

PROFILE_BINDING *profile_binding_space (n)
int n;		/* length of binding name in characters */
{
	char *space;
	PROFILE_BINDING *b = (PROFILE_BINDING *)0;

	if (space = calloc(1, sizeof(PROFILE_BINDING) + n + 1)) {
		b = (PROFILE_BINDING *)space;
		b->name = space + sizeof(PROFILE_BINDING);
	}
	return(b);
}

PROFILE_VALUE *profile_value_space (n)
int n;
{
	char *space;
	PROFILE_VALUE *v = (PROFILE_VALUE *)0;

	if (n > 0) {
		if (space = calloc(1, sizeof(PROFILE_VALUE) + n + 1)) {
			v = (PROFILE_VALUE *)space;
			v->value.s = space + sizeof(PROFILE_VALUE);
		}
	} else
		v = (PROFILE_VALUE *)calloc(1, sizeof(PROFILE_VALUE));
	return(v);
}
