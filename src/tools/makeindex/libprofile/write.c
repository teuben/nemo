static char *trwsccs = "@(#)write.c	1.1 (TRW) 1/14/86";
#include <stdio.h>
#include <ctype.h>
#include "profile.h"

profile_write_stanza (f, s)
FILE *f;
PROFILE_STANZA *s;
{
	write_markers(f, s->marker);
	fprintf(f, "{\n");
	write_bindings(f, s->binding);
	fprintf(f, "}\n");
}

static write_markers (f, m)
FILE *f;
PROFILE_MARKER *m;
{
	for (; m; m = m->next)
		fprintf(f, "%s\n", m->text);
}

static write_bindings (f, b)
FILE *f;
PROFILE_BINDING *b;
{
	while (b) {
		fprintf(f, "\t%s", b->name);
		write_values(f, b->value);
		fputc('\n', f);
		b = b->next;
	}
}

static write_values (f, v)
FILE *f;
PROFILE_VALUE *v;
{
	char scratch[PROFILE_MAX_TEXT+1];

	for (; v; v = v->next)
		switch (v->class) {
		case PROFILE_INTEGER:
			fprintf(f, " %D", v->value.i);
			continue;
		case PROFILE_FLOAT:
			fprintf(f, " %G", v->value.f);
			continue;
		case PROFILE_STRING:
			unparse_string(v->value.s, scratch);
			fprintf(f, " \"%s\"", scratch);
			continue;
		case PROFILE_CHARACTER:
			unparse_character(v->value.c, scratch);
			fprintf(f, " '%s'", scratch);
			continue;
		case PROFILE_OCTAL:
			fprintf(f, " 0o%O", v->value.i);
			continue;
		case PROFILE_HEX:
			fprintf(f, " 0x%X", v->value.i);
			continue;
		case PROFILE_OTHER:
			fprintf(f, " %s", v->value.s);
			continue;
		}
}

static int unparse_string (from, to)
char *from;
char *to;
{
	char *x = to;

	for (; *from; from++)
		switch (*from) {
		case '\b':		/* backspace */
			*x++ = '\\';
			*x++ = 'b';
			continue;
		case '\f':		/* formfeed */
			*x++ = '\\';
			*x++ = 'f';
			continue;
		case '\n':		/* newline */
			*x++ = '\\';
			*x++ = 'n';
			continue;
		case '\r':
			*x++ = '\\';
			*x++ = 'r';
			continue;
		case '\t':		/* horizontal tab */
			*x++ = '\\';
			*x++ = 't';
			continue;
		case '\\':		/* backslash */
			*x++ = '\\';
			*x++ = '\\';
			continue;
		case '"':		/* double quote */
			*x++ = '\\';
			*x++ = '"';
			continue;
		case '^':
			*x++ = '\\';
			*x++ = '^';
			continue;
		default:
			if (isascii(*from))
				if (iscntrl(*from)) {
					sprintf(x, "^%c", *from == '\177' ? '?' : *from + '@');
					x += 2;
				} else
					*x++ = *from;
			else {
				sprintf(x, "\\%03o", *from);
				x += 4;
			}
			continue;
		}
	*x = '\0';
	return(x - to);
}

static int unparse_character (from, to)
char from;
char *to;
{
	char *x = to;

	switch (from) {
	case '\b':		/* backspace */
		*x++ = '\\';
		*x++ = 'b';
		break;
	case '\f':		/* formfeed */
		*x++ = '\\';
		*x++ = 'f';
		break;
	case '\n':		/* newline */
		*x++ = '\\';
		*x++ = 'n';
		break;
	case '\r':
		*x++ = '\\';
		*x++ = 'r';
		break;
	case '\t':		/* horizontal tab */
		*x++ = '\\';
		*x++ = 't';
		break;
	case '\\':		/* backslash */
		*x++ = '\\';
		*x++ = '\\';
		break;
	case '\'':		/* single quote */
		*x++ = '\\';
		*x++ = '\'';
		break;
	case '^':
		*x++ = '\\';
		*x++ = '^';
		break;
	default:
		if (isascii(from))
			if (iscntrl(from)) {
				sprintf(x, "^%c", from == '\177' ? '?' : from + '@');
				x += 2;
			} else
				*x++ = from;
		else {
			sprintf(x, "\\%03o", from);
			x += 4;
		}
		break;
	}
	*x = '\0';
	return(x - to);
}

/*
 * write out a linked list of stanzas
 */
profile_write_profile(f, s)
FILE *f;
PROFILE_STANZA *s;
{
	PROFILE_STANZA *x;

	for (x = s; x != NULL; x = x->next) {
		profile_write_stanza(f, x);
		if (x->next == s)
			break;
	}
}
