static char *trwrcs = "$Header$";
static char *trwsccs = "@(#)read.c	1.1 (TRW) 1/14/86";
#include <stdio.h>
#include <ctype.h>
#include "profile.h"

#define isoctal(d) ('0' <= d && d <= '7')
#define ishex(x) (isdigit(x) || ('a' <= x && x <= 'f') || ('A' <= x && x <= 'F'))
#define isprime(c) (c == '\'')
#define isbackslash(c) (c == '\\')
#define iscaret(c) (c == '^')

extern char *strcpy();
extern PROFILE_STANZA *profile_stanza_space();
extern PROFILE_MARKER *profile_marker_space();
extern PROFILE_BINDING *profile_binding_space();
extern PROFILE_VALUE *profile_value_space();

static PROFILE_BINDING *get_binding();
static PROFILE_BINDING *get_bindings();
static PROFILE_MARKER *get_marker();
static PROFILE_MARKER *get_markers();
static PROFILE_VALUE *get_value();
static PROFILE_VALUE *get_values();
static char parse_character();
static PROFILE_VALUE *parse_value();

PROFILE_STANZA *profile_read_stanza (f)
FILE *f;
{
	PROFILE_STANZA *stanza;

	stanza = profile_stanza_space();
	if (stanza == NULL)
		return(NULL);
	stanza->marker = get_markers(f);
	if (get_open_bindings(f))
		stanza->binding = get_bindings(f);
	else {
		profile_free_stanza(stanza);
		return(NULL);
	}
	if (get_close_bindings(f))
		return(stanza);
	else {
		profile_free_stanza(stanza);
		return(NULL);
	}
}

/* Returns the list of markers at the head of a stanza. */
static PROFILE_MARKER *get_markers (f)
FILE *f;
{
	PROFILE_MARKER *head, *tail, *m;

	head = tail = NULL;
	while (m = get_marker(f))
		if (tail) {
			tail->next = m;
			m->previous = tail;
			tail = m;
		} else
			head = tail = m;
	return(head);
}

/* Returns the next marker from the head of the stanza. */
static PROFILE_MARKER *get_marker (f)
FILE *f;
{
	int n;
	PROFILE_MARKER *m;
	char scratch[PROFILE_MAX_TEXT+1];

	for (;;)
		if (n = get_name_text(f, scratch)) {
			if ((m = profile_marker_space(n)) == NULL)
				return(NULL);
			strcpy(m->text, scratch);
			return(m);
		} else if (get_end_of_line(f))
			continue;
		else
			return(NULL);
}

/* Returns the list of bindings in the body of the stanza. */
static PROFILE_BINDING *get_bindings (f)
FILE *f;
{
	PROFILE_BINDING *head, *tail, *b;

	head = tail = NULL;
	while (b = get_binding(f))
		if (tail) {
			tail->next = b;
			b->previous = tail;
			tail = b;
		} else
			head = tail = b;
	return(head);
}

/* Returns the next binding in the body of the stanza. */
static PROFILE_BINDING *get_binding (f)
FILE *f;
{
	int n;
	PROFILE_BINDING *b;
	char scratch[PROFILE_MAX_TEXT+1];

	for (;;)
		if (n = get_name_text(f, scratch)) {
			if ((b = profile_binding_space(n)) == NULL)
				return(NULL);
			strcpy(b->name, scratch);
			break;
		} else if (get_end_of_line(f))
			continue;
		else
			return(NULL);
	b->value = get_values(f);
	return(b);
}

/* Returns the list of values following the name of the binding. */
static PROFILE_VALUE *get_values (f)
FILE *f;
{
	PROFILE_VALUE *head, *tail, *v;

	head = tail = NULL;
	while (v = get_value(f))
		if (tail) {
			tail->next = v;
			v->previous = tail;
			tail = v;
		} else
			head = tail = v;
	return(head);
}

/* Returns the next value in the binding. */
static PROFILE_VALUE *get_value (f)
FILE *f;
{
	char text[PROFILE_MAX_TEXT+1];
	int n;

	for (;;)
		if (n = get_value_text(f, text))
			return(parse_value(text, n));
		else if (get_end_of_line(f))
			return(NULL);
		else
			return(NULL);
}

/*
 * Reads the text of the next value (if any) in the binding.  Returns
 * the length of the literal text in characters.
 */
static int get_value_text (f, text)
FILE *f;
char *text;
{
	register int c;
	char *s = text;

	while ((c = getc(f)) != EOF)
		switch (c) {
		case '\b': case '\f':
		case '\r': case '\t': case ' ':
			/* white space terminates any text gathered so far */
			if (s > text) {
				*s = '\0';
				return(s - text);
			}
			continue;

		case '\n':
			/* newline terminates a binding */
			ungetc(c, f);
			*s = '\0';
			return(s - text);

		case '#':
			/* gobble up the comment */
			while ((c = getc(f)) != EOF && c != '\n')
				continue;
			if (c == '\n')
				ungetc(c, f);
			*s = '\0';
			return(s - text);

		case '{': case '}':
			ungetc(c, f);
			*s = '\0';
			return(s - text);

		case '"':		/* string quotes */
			ungetc(c, f);
			if (s > text) {
				*s = '\0';
				return(s - text);
			} else
				return(get_string(f, s));

		case '\'':		/* character quotes */
			ungetc(c, f);
			if (s > text) {
				*s = '\0';
				return(s - text);
			} else
				return(get_character(f, s));

		case '\\':		/* newline escape */
			c = getc(f);
			if (c == '\n') {
				if (s > text) {
					*s = '\0';
					return(s - text);
				}
				continue;	/* just like a blank */
			}
			ungetc(c, f);
			if (s < &text[PROFILE_MAX_TEXT])
				*s++ = '\\';
			continue;

		default:
			if (s < &text[PROFILE_MAX_TEXT])
				*s++ = c;
			continue;
		}
	*s = '\0';
	return(s - text);
}

/* Digests the raw value text returning a new value structure. */
static PROFILE_VALUE *parse_value (s, length)
char *s;		/* literal text */
int length;		/* in characters */
{
	PROFILE_VALUE *v;

	if (is_integer(s)) {
		if ((v = profile_value_space(0)) == NULL)
			return(NULL);
		v->class = PROFILE_INTEGER;
		sscanf(s, "%D", &v->value.i);
	} else if (is_octal(s)) {
		if ((v = profile_value_space(0)) == NULL)
			return(NULL);
		v->class = PROFILE_OCTAL;
		/* skip the `0o' prefix */
		sscanf(s+2, "%O", &v->value.i);
	} else if (is_hex(s)) {
		if ((v = profile_value_space(0)) == NULL)
			return(NULL);
		v->class = PROFILE_HEX;
		/* skip the `0x' prefix */
		sscanf(s+2, "%X", &v->value.i);
	} else if (is_string(s)) {
		/* be careful when dealing with the empty string "" */
		if ((v = profile_value_space(length > 2 ? length - 2 : 1)) == NULL)
			return(NULL);
		v->class = PROFILE_STRING;
		/* erase the terminating double quote */
		s[length - 1] = '\0';
		/* skip past the initial quote */
		parse_string(s + 1, v->value.s);
	} else if (is_character(s)) {
		if ((v = profile_value_space(0)) == NULL)
			return(NULL);
		v->class = PROFILE_CHARACTER;
		/* erase the end single quote */
		s[length - 1] = '\0';
		v->value.c = parse_character(s + 1);
	} else if (is_float(s)) {
		if ((v = profile_value_space(0)) == NULL)
			return(NULL);
		v->class = PROFILE_FLOAT;
		sscanf(s, "%E", &v->value.f);
	} else {
		if ((v = profile_value_space(length)) == NULL)
			return(NULL);
		v->class = PROFILE_OTHER;
		strcpy(v->value.s, s);
	}
	return(v);
}

/* Converts a string literal to the internal representation. */
static parse_string (source, result)
char *source;
char *result;
{
	for (; *source; source++)
		if (*source == '\\')
			switch (*++source) {
			case 'b':			/* backspace */
				*result++ = '\b';
				continue;
			case 'f':			/* formfeed */
				*result++ = '\f';
				continue;
			case 'n':			/* newline */
				*result++ = '\n';
				continue;
			case 'r':			/* carriage return */
				*result++ = '\r';
				continue;
			case 't':			/* horizontal tab */
				*result++ = '\t';
				continue;
			case '\'':			/* single quote */
				*result++ = '\'';
				continue;
			case '"':			/* double quote */
				*result++ = '"';
				continue;
			case '\\':			/* backslash */
				*result++ = '\\';
				continue;			
			case '^':			/* caret */
				*result++ = '^';
				continue;
			case '0': case '1':		/* octal constant */
			case '2': case '3':
			case '4': case '5':
			case '6': case '7':
				source += parse_octal(source, result) - 2;
				result++;
				continue;
			default:
				*result++ = *source;	/* ignore backslash */
			}
		else if (*source == '^') {	/* control escape */
			char c = *++source;
			*result++ = ('@' <= c && c <= '_') ? c - '@' :
				    (c == '?') ? '\177' : c;
			continue;
		} else
			*result++ = *source;
	*result = '\0';
}

/* Converts a character literal to the internal representation. */
static char parse_character (source)
char *source;
{
	char c;

	if (*source == '\\')
		switch (*++source) {
		case 'b':			/* backspace */
			return('\b');
		case 'f':			/* formfeed */
			return('\f');
		case 'n':			/* newline */
			return('\n');
		case 'r':			/* carriage return */
			return('\r');
		case 't':			/* horizontal tab */
			return('\t');
		case '\'':			/* single quote */
			return('\'');
		case '\\':			/* backslash */
			return('\\');
		case '^':
			return('^');
		case '0': case '1':		/* octal constant */
		case '2': case '3':
		case '4': case '5':
		case '6': case '7':
			parse_octal(source, &c);
			return(c);
		default:
			return(*source);	/* ignore backslash */
		}
	else if (*source == '^') {	/* control escape */
		c = *++source;
		return(('@' <= c && c <= '_') ? c - '@' : (c == '?') ? '\177' : c);
	} else
		return(*source);
}

/* Converts an octal escape `\ddd' to its byte representation. */
static int parse_octal (source, result)
char *source;
char *result;
{
	int count;
	char byte = '\0';
	char digit;

	for (count = 1; count <= 3; count++) {
		digit = *source++;
		if ('0' <= digit && digit <= '7')
			byte = (byte * 8) + (digit - '0');
		else
			break;
	}
	*result = byte;
	return(count);
}

/*
 * Reads the literal text for markers and binding names.  Returns the
 * length in characters of the literal text.
 */
static int get_name_text (f, text)
FILE *f;
char *text;
{
	register int c;
	char *s = text;

	while ((c = getc(f)) != EOF)
		switch (c) {
		case '\b': case '\f':
		case '\r': case '\t': case ' ':
			/* white space terminates text gathered so far */
			if (s > text) {
				*s = '\0';
				return(s - text);
			}
			continue;

		case '\n':
			ungetc(c, f);
			*s = '\0';
			return(s - text);

		case '#':
			/* gobble up the comment */
			while ((c = getc(f)) != EOF && c != '\n')
				continue;
			if (c == '\n')
				ungetc(c, f);
			*s = '\0';
			return(s - text);

		case '{': case '}':
			ungetc(c, f);
			*s = '\0';
			return(s - text);

		case '[':
			/* sets may contain embedded white space */
			if (s + 1 < &text[PROFILE_MAX_TEXT]) {
				*s++ = '[';
				if ((c = getc(f)) != EOF)
				*s++ = c;
			}
			while ((c = getc(f)) != EOF) {
				if (s < &text[PROFILE_MAX_TEXT])
					*s++ = c;
				if (c == ']')
					break;
			}
			continue;

		case '\\':
			c = getc(f);
			if (c == '\n') {
				if (s > text) {
					*s = '\0';
					return(s - text);
				}
				continue;	/* just like a blank */
			}
			ungetc(c, f);
			if (s < &text[PROFILE_MAX_TEXT])
				*s++ = '\\';
			continue;

		default:
			if (s < &text[PROFILE_MAX_TEXT])
				*s++ = c;
			continue;
		}
	*s = '\0';
	return(s - text);
}

/* Returns non-zero on end of line and zero otherwise. */
static int get_end_of_line (f)
FILE *f;
{
	int c;

	if ((c = getc(f)) == '\n')
		return(1);
	ungetc(c, f);
	return(0);
}

/* Returns non-zero on seeing `{' and zero otherwise. */
static int get_open_bindings (f)
FILE *f;
{
	int c;

	if ((c = getc(f)) == '{')
		return(1);
	ungetc(c, f);
	return(0);
}

/* Returns non-zero on seeing `}' and zero otherwise. */ 
static int get_close_bindings (f)
FILE *f;
{
	int c;

	if ((c = getc(f)) == '}')
		return(1);
	ungetc(c, f);
	return(0);
}

/* Reads a string literal returning the length of the literal text in characters */
static int get_string (f, text)
FILE *f;
char *text;
{
	register int c;
	char *s = text;

	/* the first double quote is guaranteed */
	*s++ = getc(f);
	while ((c = getc(f)) != EOF)
		switch (c) {
		case '\\':
			if (s < &text[PROFILE_MAX_TEXT])
				*s++ = c;
			c = getc(f);
			if (c == EOF)
				return(s - text);
			else if (c == '\n') {
				ungetc(c, f);
				return(s - text);
			} else if (s < &text[PROFILE_MAX_TEXT])
				*s++ = c;
			continue;

		case '"':
			if (s < &text[PROFILE_MAX_TEXT])
				*s++ = c;
			*s = '\0';
			return(s - text);

		case '\n':
			ungetc(c, f);
			*s = '\0';
			return(s - text);

		default:
			if (s < &text[PROFILE_MAX_TEXT])
				*s++ = c;
			continue;
		}
	*s = '\0';
	return(s - text);
}

/* Reads a character literal returning the length of the literal text in characters. */
static int get_character (f, text)
FILE *f;
char *text;
{
	register int c;
	char *s = text;

	/* the first single quote is guaranteed */
	*s++ = getc(f);
	while ((c = getc(f)) != EOF)
		switch (c) {
		case '\\':
			if (s < &text[PROFILE_MAX_TEXT])
				*s++ = c;
			c = getc(f);
			if (c == EOF)
				return(s - text);
			else if (c == '\n') {
				ungetc(c, f);
				return(s - text);
			} else if (s < &text[PROFILE_MAX_TEXT])
				*s++ = c;
			continue;
		case '\'':
			if (s < &text[PROFILE_MAX_TEXT])
				*s++ = c;
			*s = '\0';
			return(s - text);
		case '\n':
			ungetc(c, f);
			*s = '\0';
			return(s - text);
		default:
			if (s < &text[PROFILE_MAX_TEXT])
				*s++ = c;
			continue;
		}
	*s = '\0';
	return(s - text);
}

/* all regular expressions below are in lex notation */

/* returns non-zero iff -?[0-9]+ matches */
static int is_integer (s)
char *s;
{
	char *x;

	/* -? */
	if (*s == '-')
		s++;
	/* [0-9]+ */
	for (x = s; isdigit(*s); s++)
		continue;
	return(s > x && !*s);
}

/* returns non-zero iff 0[oO][0-7]+ matches */
static int is_octal (s)
char *s;
{
	char *x;

	/* 0 */
	if (*s == '0')
		s++;
	else
		return(0);
	/* [oO] */
	if (*s == 'o' || *s == 'O')
		s++;
	else
		return(0);
	/* [0-7]+ */
	for (x = s; isoctal(*s); s++)
		continue;
	return(s > x && !*s);
}

/* returns non-zero iff 0[xX][0-9a-fA-F]+ matches */
static int is_hex (s)
char *s;
{
	char *x;

	/* 0 */
	if (*s == '0')
		s++;
	else
		return(0);
	/* [xX] */
	if (*s == 'x' || *s == 'X')
		s++;
	else
		return(0);
	/* [0-9a-fA-F]+ */
	for (x = s; ishex(*s); s++)
		continue;
	return(s > x && !*s);
}

/* returns non-zero iff [eE][-+]?[0-9]+ matches */
static int is_exponent (s)
char *s;
{
	char *x;

	/* [eE] */
	if (*s == 'e' || *s == 'E')
		s++;
	else
		return(0);
	/* [-+]? */
	if (*s == '-' || *s == '+')
		s++;
	/* [0-9]+ */
	for (x = s; isdigit(*s); s++)
		continue;
	return(s > x && !*s);
}

static int is_float (s)
char *s;
{
	return(is_integer_part_float(s) ||
	       is_fractional_part_float(s) ||
	       is_power_float(s));
}

/* returns non-zero iff -?[0-9]+"."[0-9]*({exponent})? matches */
static int is_integer_part_float (s)
char *s;
{
	char *x;

	/* -? */
	if (*s == '-')
		s++;
	/* [0-9]+"." */
	for (x = s; isdigit(*s); s++)
		continue;
	if (x == s || *s != '.')
		return(0);
	/* [0-9]* */
	for (s++; isdigit(*s); s++)
		continue;
	/* ({exponent})? */
	return(*s ? is_exponent(s) : 1);
}

/* returns non-zero iff -?"."[0-9]+({exponent})? matches */
static int is_fractional_part_float (s)
char *s;
{
	char *x;

	/* -? */
	if (*s == '-')
		s++;
	/* "." */
	if (*s == '.')
		s++;
	else
		return(0);
	/* [0-9]+({exponent})? */
	for (x = s; isdigit(*s); s++)
		continue;
	return(s > x ? !*s || is_exponent(s) : 0);
}

/* returns non-zero iff -?[0-9]+{exponent} matches */
static int is_power_float (s)
char *s;
{
	char *x;

	/* -? */
	if (*s == '-')
		s++;
	/* [0-9]+{exponent} */
	for (x = s; isdigit(*s); s++)
		continue;
	return(s > x ? is_exponent(s) : 0);
}

/* returns non-zero iff '[^^\]' | '\^.' | '\\\\' | '\\'' | '\[0-7]{1-3}' matches */
static int is_character (s)
char *s;
{
	char *x;

	if (isprime(*s))
		s++;
	else
		return(0);
	if (isbackslash(*s)) {
		s++;
		if ((isbackslash(s[0]) || isprime(s[0]) || !isdigit(s[0])) &&
		    isprime(s[1]) && !s[2])
			return(1);
		for (x = s; isoctal(*s); s++)
			continue;
		return(x < s && s < (x+4) && isprime(s[0]) && !s[1]);
	} else if (iscaret(*s))
		s++;
	return(isprint(s[0]) && isprime(s[1]) && !s[2]);
}

/* returns non-zero iff s is a string constant */
static int is_string (s)
char *s;
{
	char *x;

	if (*s != '"')
		return(0);
	for (s++; *s; s++) {
		if (*s == '"')
			return(!*++s);	/* quote must be followed by null */
		if (isbackslash(*s) || iscaret(*s)) {
			if (*++s)
				continue;	/* legal escape */
			return(0);	/* null follows \ or ^ */
		}
	}
	return(0);
}

/*
 * read an entire profile, making a bidirectional
 * circularly linked list
 * returns pointer to the first stanza or NULL on error
 */
PROFILE_STANZA *profile_read_profile(f)
FILE *f;
{
	PROFILE_STANZA *head = NULL;
	PROFILE_STANZA *tail = NULL;
	PROFILE_STANZA *x = NULL;

	while ((x = profile_read_stanza(f)) != NULL) {
		if (head == NULL)
			head = tail = x;
		else {
			tail->next = x; 
			x->previous = tail;
			tail = x;
		}
	}
	if (head != NULL) { 
		tail->next = head;
		head->previous = tail;
	}
	return(head);
}
