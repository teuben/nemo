static char *trwsccs = "@(#)boolean.c	1.1 (TRW) 1/14/86";
#include <ctype.h>
#include "profile.h"

static char *Yes[] = {
	"yes",
	"on",
	"true",
	"enable",
	"available",
	"present",
	0
};

static char *No[] = {
	"no",
	"off",
	"false",
	"disable",
	"unavailable",
	"absent",
	0
};

int profile_boolean (v)
PROFILE_VALUE *v;
{
	char x[16];
	int i;

	if (v == 0)
		return(0);

	switch (v->class) {
	case PROFILE_OTHER:
	case PROFILE_STRING:
		strncpy(x, v->value.s, sizeof(x)-1);
		x[sizeof(x)-1] = 0;
		downshift(x);
		for (i = 0; Yes[i]; i++)
			if (strcmp(x, Yes[i]) == 0)
				return(1);
			else if (strcmp(x, No[i]) == 0)
				return(0);
		return(-1);	/* unknown string */

	case PROFILE_HEX:
	case PROFILE_INTEGER:
	case PROFILE_OCTAL:
		return(v->value.i != 0);

	case PROFILE_CHARACTER:
		return(v->value.c != 0);

	case PROFILE_FLOAT:
		return(v->value.f != 0.0);

	default:
		return(-1);	/* unknown class */
	}
}

/* downshift a string in place */
static downshift (s)
char *s;
{
	for (; *s; s++)
		if (isupper(*s))
			*s = tolower(*s);
}
