/* @(#)profile.h	1.1 (TRW) 1/14/86 */
typedef struct PROFILE_VALUE {
	char class;
	union {
		long int i;
		double f;
		char c;
		char *s;
	} value;
	struct PROFILE_VALUE *previous;
	struct PROFILE_VALUE *next;
} PROFILE_VALUE;

typedef struct PROFILE_BINDING {
	char *name;
	PROFILE_VALUE *value;
	struct PROFILE_BINDING *previous;
	struct PROFILE_BINDING *next;
} PROFILE_BINDING;

typedef struct PROFILE_MARKER {
	char *text;
	struct PROFILE_MARKER *previous;
	struct PROFILE_MARKER *next;
} PROFILE_MARKER;

typedef struct PROFILE_STANZA {
	PROFILE_MARKER *marker;
	PROFILE_BINDING *binding;
	struct PROFILE_STANZA *previous;
	struct PROFILE_STANZA *next;
} PROFILE_STANZA;

/* classes */
#define PROFILE_INTEGER 01
#define PROFILE_FLOAT 02
#define PROFILE_STRING 03
#define PROFILE_CHARACTER 04
#define PROFILE_OTHER 05
#define PROFILE_OCTAL 06
#define PROFILE_HEX 07

/* no single lexical element may exceed this size in characters */
#define PROFILE_MAX_TEXT 255

PROFILE_STANZA *profile_read_stanza();
PROFILE_STANZA *profile_read_profile();
PROFILE_MARKER *profile_has_marker();
PROFILE_STANZA *profile_has_stanza();
PROFILE_BINDING *profile_has_binding();
PROFILE_STANZA *profile_stanza_space();
PROFILE_MARKER *profile_marker_space();
PROFILE_BINDING *profile_binding_space();
PROFILE_VALUE *profile_value_space();
