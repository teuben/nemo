static char *trwsccs = "@(#)free.c	1.1 (TRW) 1/14/86";
#include "profile.h"

profile_free_profile (s)
PROFILE_STANZA *s;
{
	PROFILE_STANZA *x;

	for (x = s; x != (PROFILE_STANZA *)0 && x != s; x = x->next)
		profile_free_stanza(x);
}

profile_free_stanza (s)
PROFILE_STANZA *s;
{
	free_markers(s->marker);
	free_bindings(s->binding);
	free(s);
}

static free_markers (m)
PROFILE_MARKER *m;
{
	PROFILE_MARKER *x;

	for (; m; m = x) {
		x = m->next;
		free(m);
	}
}

static free_bindings (b)
PROFILE_BINDING *b;
{
	PROFILE_BINDING *x;

	for (; b; b = x) {
		x = b->next;
		free_values(b->value);
		free(b);
	}
}

static free_values (v)
PROFILE_VALUE *v;
{
	PROFILE_VALUE *x;

	for (; v; v = x) {
		x = v->next;
		free(v);
	}
}
