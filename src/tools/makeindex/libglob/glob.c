static char *trwsccs= "@(#)glob.c	1.1 (TRW) 1/14/86";
#include "glob.h"

#define SLOP 5
#define MAX_SET 0177

/* control codes for regular expression evaluation */
#define PATTERN_ANY '?'
#define PATTERN_CHARACTER 'X'
#define PATTERN_END '$'
#define PATTERN_SET '['
#define PATTERN_SET_MEMBER 'M'
#define PATTERN_SET_RANGE '-'
#define PATTERN_STAR '*'

/*
 * Examples (=> denotes `compiles into')
 *
 *	a	=>	Xa
 *	?	=>	?
 *	[x0-9]	=>	[^EMx-09	(^E is control-E)
 *	*	=>	*
 *	END	=>	$
 *
 *	a?[x0-9]* => Xa?[^EMx-09*$
 */

glob_compile (pattern, buffer)
char *pattern;
char *buffer;	/* compiled pattern */
{
	char *x;	/* pointer into compiled pattern */
	int c;
	int result;

	if (pattern == 0 || pattern[0] == 0)
		return(GLOB_PATTERN_EMPTY);

	x = buffer;
	while (x < &buffer[GLOB_MAX_PATTERN - SLOP]) {
		c = *pattern++;
		if (c == 0) {
			*x++ = PATTERN_END;
			return(GLOB_OK);
		}

		switch (c) {
		case '?':
			*x++ = PATTERN_ANY;
			continue;

		case '[':
			if ((result = compile_set(pattern, x, &buffer[GLOB_MAX_PATTERN - SLOP])) < 0)
				return(result);
			pattern += result + 1;
			x += x[1] + 2;
			continue;

		case '*':
			*x++ = PATTERN_STAR;
			continue;

		default:
			*x++ = PATTERN_CHARACTER;
			*x++ = c;
			continue;
		}
	}
	return(GLOB_PATTERN_TOO_BIG);
}

int glob_execute (pattern, s)
char *pattern;	/* compiled pattern */
char *s;	/* string to be matched against */
{
	char *current;
	int result;

	for (;;)
		switch (*pattern++) {
		case PATTERN_ANY:
			if (*s++)
				continue;
			return(0);

		case PATTERN_CHARACTER:
			if (*pattern++ == *s++)
				continue;
			return(0);

		case PATTERN_END:
			return(*s == 0);

		case PATTERN_SET:
			if ((result = in_set(pattern, *s++)) == 1) {
				pattern += *pattern + 1;
				continue;
			}
			return(result);

		case PATTERN_STAR:
			current = s;
			while (*s++)
				continue;
			do {
				s--;
				if (result = glob_execute(pattern, s))
					return(result);
			} while (s > current);
			return(0);

		default:
			return(GLOB_EXECUTION_ERROR);
		}
}

int glob_match (pattern, s)
char *pattern;
char *s;
{
	int result;
	char buffer[GLOB_MAX_PATTERN];

	if ((result = glob_compile(pattern, buffer)) < 0)
		return(result);
	else
		return(glob_execute(buffer, s));
}

/* returns 1 if character c is member of set and 0 otherwise */
static int in_set (set, c)
char *set;	/* compiled set pattern */
char c;
{
	int n;

	if (c == 0)
		return(0);
	n = *set++;
	while (n > 0)
		switch (*set++) {
		case PATTERN_SET_MEMBER:
			if (*set++ == c)
				return(1);
			n -= 2;
			continue;

		case PATTERN_SET_RANGE:
			if (*set++ <= c && c <= *set++)
				return(1);
			n -= 3;
			continue;

		default:
			return(GLOB_EXECUTION_ERROR);
		}
	return(0);
}

#define IS_RANGE(s) (s[1] && s[2] && s[1] == '-' && s[2] != ']')

/* compiles a set returning the number of pattern characters consumed */
static int compile_set (pattern, x, limit)
char *pattern;
char *x;
char *limit;
{
	char *slot;	/* size of set goes here */
	int size;	/* number of bytes in compiled set */
	char *start = pattern;

	if (*pattern == 0)
		return(GLOB_BRACKET_MISSING);

	*x++ = PATTERN_SET;
	slot = x++;
	size = 0;

	if (IS_RANGE(pattern)) {
		if (pattern[0] > pattern[2])	/* pattern[1] == '-' */
			return(GLOB_RANGE_INVERTED);
		*x++ = PATTERN_SET_RANGE;
		*x++ = pattern[0];
		*x++ = pattern[2];
		pattern += 3;
		size += 3;
	} else {
		*x++ = PATTERN_SET_MEMBER;
		*x++ = *pattern++;
		size += 2;
	}

	while (*pattern != ']' && x < limit) {
		if (*pattern == 0)
			return(GLOB_BRACKET_MISSING);
		if (IS_RANGE(pattern)) {
			if (pattern[0] > pattern[2])	/* pattern[1] == '-' */
				return(GLOB_RANGE_INVERTED);
			*x++ = PATTERN_SET_RANGE;
			*x++ = pattern[0];
			*x++ = pattern[2];
			pattern += 3;
			size += 3;
		} else {
			*x++ = PATTERN_SET_MEMBER;
			*x++ = *pattern++;
			size += 2;
		}
	}
	if (size > MAX_SET)
		return(GLOB_SET_TOO_BIG);
	*slot = size;
	return(pattern - start);
}
