/* $Header */
/* Yacc parser for LaTeX index processor */
/* Roman numeral code written by John Renner (adobe!renner@decwrl.dec.com) */
%{
#include <stdio.h>
#include "standard.h"
#include "profile.h"
#include <ctype.h>
#include <pwd.h>
#define TABLEINCREMENT	50		/* Number of additional entries added when expanding a table */
#define eq(s,t)		(!strcmp((s),(t)))
#define odd(i)		(((i) % 2) == 1)
#define ITEMDEPTH	3		/* Number of nestings of \item's, \subitem's, \subsubitem's, etc. */
char *ItemRep[] = { "\\item", "\\subitem", "\\subsubitem", NULL }; /* and their representation */
char	*calloc(), *realloc();
enum TokenType	{controlword, controlsymbol, string, integer, roman, comma, obrace, cbrace, whitespace};
struct IndexEntry {
    char		literal[81];	/* Literal representation of index entry */
    char		alphabetic[81];	/* Alphabetic representation for sorting of index entry */
    struct Token	*tokenlist;	/* Doubly linked token list */
    struct IndexEntry	*subitem;	/* Pointer to subitem table, subsubitem table, etc */
    int			subitemcount;	/* Number of items in subitem table */
    int			subitemtabsize;	/* Subitem table size currently allocated */
    struct PageNoTable	*pagenos;	/* List of page numbers */
    int			pagetablecount;	/* Number of items in page number table */
    int			pagetablesize;	/* Size of page number table currently allocated */
};
struct Token {
    enum TokenType	type;		/* Token type */
    char		lexeme[81];	/* Representation of all the token types */
    struct Token	*prev, *next;
};
struct PageNoTable {
    int			number;		/* Page number */
    Boolean		range;		/* True if this is the beginning of a range */
    Boolean		isroman;	/* True if this was a roman numeral */
};
struct IndexEntry	*IndexTable = NULL;	/* Table of primary index entries */
int			IndexTableCount = 0;	/* Count of number of elements used in index table */
int			IndexTableSize = 0;	/* Current allocated size of index table */
int			ExitStatus = SUCCEED;	/* Guess */
int			LineNo = 1;		/* Line number at start of token */
int			EndLineNo = 1;		/* Line number at end of token */
Boolean			Label = FALSE;		/* True if -l option given */
Boolean			Range;			/* True if this \indexentry is a range */
PROFILE_STANZA		*SortStanza = NULL;	/* Alphabetize stanza */
extern int		optind;			/* From getopt(3) */
extern char		*optarg;
char			*Whoami;		/* argv[0] */
char			*Usage = "Usage: %s [-l] [-f alphabetizefile] [file...]\n";
char			*Marker[] = { "alphabetize", NULL };	/* Markers for alphabetize stanza */
char			IdxFileName[81];	/* .idx file name */
char			Literal[81];		/* Literal string of key */
char			Alphabetic[81];		/* Alphabetic string of key */
FILE			*InputFile;		/* Current input file */
FILE			*OutputFile;		/* Current output file */
struct Token		*CurKey;		/* Current key we are constructing */
struct IndexEntry	**CurSearchTable;	/* Current table to search for match */
struct IndexEntry	*CurEntry;		/* Current table entry */
struct IndexEntry	*PrevEntry;		/* Previous Entry */
%}
%union {
    char		value[81];
    struct Token	*t;
}
%token <value> ROMAN CONTROLSEQUENCE INTEGER WHITESPACE STRING INDEXENTRY
%type <t> noncommaelement anyelement anyelements
%%
indexfile :
    optwhitespace
    indexentries
	{
	    sort(IndexTable, IndexTableCount);
	    fprintf(OutputFile, "\\begin{theindex}\n");
	    if (Label) {
		fprintf(OutputFile,
 "\\newcommand{\\largeletter}[1]{{\\pagebreak[2]\\Large\\hspace{-.5in}\\parbox[t]{.5in}{\\makebox[.35in][r]");
		fprintf(OutputFile, "{\\uppercase{#1}}}\\nopagebreak[4]\\vspace{-1.5ex}}}\n");
	    }
	    printindexentries(IndexTable, IndexTableCount, 1);
	    fprintf(OutputFile, "\\end{theindex}\n");
	}
    ;

indexentries :
    indexentries
    indexentry
    |
    indexentry
    ;

indexentry :
    INDEXENTRY
	{
	    CurSearchTable = &IndexTable, PrevEntry = NULL;
	    CurKey = NULL;
	    Range = FALSE;
	}
    optwhitespace
    '{'
    keys
    '}'
    optwhitespace
    '{'
    optwhitespace
    anumber
    optwhitespace
    '}'
    optwhitespace
    ;

anumber : 
    INTEGER
	{
	    struct PageNoTable	*p;

	    if (!(p = findpage(CurEntry->pagenos, CurEntry->pagetablecount, atoi($1), FALSE))) {
		if (CurEntry->pagetablecount >= CurEntry->pagetablesize) {
		    if (!(CurEntry->pagenos = (struct PageNoTable *)reallocate(CurEntry->pagenos, CurEntry->pagetablesize,
		    TABLEINCREMENT, sizeof(struct PageNoTable)))) {
			yyerror("memory allocation failure");
			exit(FAIL);
		    }
		    CurEntry->pagetablesize += TABLEINCREMENT;
		}
		CurEntry->pagenos[CurEntry->pagetablecount].number = atoi($1);
		CurEntry->pagenos[CurEntry->pagetablecount].isroman = FALSE;
		CurEntry->pagenos[CurEntry->pagetablecount].range = Range;
		CurEntry->pagetablecount++;
	    } else
		p->range = Range;
	}
      |
      ROMAN
	{
	    struct PageNoTable	*p;

	    if (!(p = findpage(CurEntry->pagenos, CurEntry->pagetablecount, rmtoi($1), TRUE))) {
		if (CurEntry->pagetablecount >= CurEntry->pagetablesize) {
		    if (!(CurEntry->pagenos = (struct PageNoTable *)reallocate(CurEntry->pagenos, CurEntry->pagetablesize,
		    TABLEINCREMENT, sizeof(struct PageNoTable)))) {
			yyerror("memory allocation failure");
			exit(FAIL);
		    }
		    CurEntry->pagetablesize += TABLEINCREMENT;
		}
		CurEntry->pagenos[CurEntry->pagetablecount].number = rmtoi($1);
		CurEntry->pagenos[CurEntry->pagetablecount].isroman = TRUE;
		CurEntry->pagenos[CurEntry->pagetablecount].range = Range;
		CurEntry->pagetablecount++;
	    } else
		p->range = Range;
	}
	
	
keys :
    multiplekeys
    key
    {
	struct Token	*t;

	for (t = CurKey; t->next; t = t->next)
	    ;
	if (t->type == string)
	    if (t->lexeme[strlen(t->lexeme) - 1] == '-') {
		t->lexeme[strlen(t->lexeme) - 1] = '\0';
		Range = TRUE;
	    }
	goto installkey;
    }
    ;

multiplekeys :
    multiplekeys
    key
    ','
	{
	    struct Token	*t;

installkey: strcpy(Literal, literalstring(CurKey));
	    strcpy(Alphabetic, alphabetizestring(CurKey, SortStanza));
	    if (!*CurSearchTable) {
		if (!(*CurSearchTable = (struct IndexEntry *)reallocate(*CurSearchTable, 0, TABLEINCREMENT,
		sizeof(struct IndexEntry)))) {
		    yyerror("memory allocation failure");
		    exit(FAIL);
		}
		if (!PrevEntry)
		    IndexTableSize = TABLEINCREMENT;
		else
		    PrevEntry->subitemtabsize = TABLEINCREMENT;
	    }
	    if (!(CurEntry = findentry(*CurSearchTable, PrevEntry ? PrevEntry->subitemcount : IndexTableCount, Literal))) {
		if (!PrevEntry) {
		    if (IndexTableCount >= IndexTableSize) {
			if (!(*CurSearchTable = (struct IndexEntry *)reallocate(*CurSearchTable, IndexTableSize, TABLEINCREMENT,
			sizeof(struct IndexEntry)))) {
			    yyerror("memory allocation failure");
			    exit(FAIL);
			}
			IndexTableSize += TABLEINCREMENT;
		    }
		    CurEntry = (*CurSearchTable + IndexTableCount);
		    IndexTableCount++;
		} else {
		    if (PrevEntry->subitemcount >= PrevEntry->subitemtabsize) {
			if (!(*CurSearchTable = (struct IndexEntry *)reallocate(*CurSearchTable, PrevEntry->subitemtabsize,
			TABLEINCREMENT, sizeof(struct IndexEntry)))) {
			    yyerror("memory allocation failure");
			    exit(FAIL);
			}
			PrevEntry->subitemtabsize += TABLEINCREMENT;
		    }
		    CurEntry = (*CurSearchTable + PrevEntry->subitemcount);
		    PrevEntry->subitemcount++;
		}
		strcpy(CurEntry->literal, Literal);
		strcpy(CurEntry->alphabetic, Alphabetic);
		CurKey->prev = CurEntry->tokenlist, CurEntry->tokenlist = CurKey;
		CurEntry->subitem = NULL, CurEntry->subitemcount = CurEntry->subitemtabsize = 0;
		CurEntry->pagenos = NULL, CurEntry->pagetablecount = CurEntry->pagetablesize = 0;
	    }
	    CurSearchTable = &CurEntry->subitem;
	    PrevEntry = CurEntry;
	    CurKey = NULL;
	}
    |
    /* epsilon */
    ;

key :
    key
    noncommaelement
    |
    noncommaelement
    ;

noncommaelement :
    CONTROLSEQUENCE
	{
	    if (!($$ = (struct Token *)calloc(1, sizeof(struct Token)))) {
		yyerror("memory allocation failure");
		exit(FAIL);
	    }
	    $$->type = isalpha($1[1]) ? controlword : controlsymbol;
	    strcpy($$->lexeme, $1);
	    $$->next = NULL;
	    if (!CurKey)
		$$->prev = CurKey, CurKey = $$;
	    else {
		struct Token *p;
		
		for (p = CurKey; p->next; p = p->next)
		    ;
		p->next = $$, $$->prev = p;
	    }
	    $$ = CurKey;
	}
    |
    ROMAN
	{
	    if (!($$ = (struct Token *)calloc(1, sizeof(struct Token)))) {
		yyerror("memory allocation failure");
		exit(FAIL);
	    }
	    $$->type = roman;
	    strcpy($$->lexeme, $1);
	    $$->next = NULL;
	    if (!CurKey)
		$$->prev = CurKey, CurKey = $$;
	    else {
		struct Token *p;
		
		for (p = CurKey; p->next; p = p->next)
		    ;
		p->next = $$, $$->prev = p;
	    }
	    $$ = CurKey;
	}
    |
    INTEGER
	{
	    if (!($$ = (struct Token *)calloc(1, sizeof(struct Token)))) {
		yyerror("memory allocation failure");
		exit(FAIL);
	    }
	    $$->type = integer;
	    strcpy($$->lexeme, $1);
	    $$->next = NULL;
	    if (!CurKey)
		$$->prev = CurKey, CurKey = $$;
	    else {
		struct Token *p;
		
		for (p = CurKey; p->next; p = p->next)
		    ;
		p->next = $$, $$->prev = p;
	    }
	    $$ = CurKey;
	}
    |
    WHITESPACE
	{
	    if (!($$ = (struct Token *)calloc(1, sizeof(struct Token)))) {
		yyerror("memory allocation failure");
		exit(FAIL);
	    }
	    $$->type = whitespace;
	    strcpy($$->lexeme, $1);
	    $$->next = NULL;
	    if (!CurKey)
		$$->prev = CurKey, CurKey = $$;
	    else {
		struct Token *p;
		
		for (p = CurKey; p->next; p = p->next)
		    ;
		p->next = $$, $$->prev = p;
	    }
	    $$ = CurKey;
	}
    |
    STRING
	{
	    if (!($$ = (struct Token *)calloc(1, sizeof(struct Token)))) {
		yyerror("memory allocation failure");
		exit(FAIL);
	    }
	    $$->type = string;
	    strcpy($$->lexeme, $1);
	    $$->next = NULL;
	    if (!CurKey)
		$$->prev = CurKey, CurKey = $$;
	    else {
		struct Token *p;
		
		for (p = CurKey; p->next; p = p->next)
		    ;
		p->next = $$, $$->prev = p;
	    }
	    $$ = CurKey;
	}
    |
    '{'
	{
	    if (!($$ = (struct Token *)calloc(1, sizeof(struct Token)))) {
		yyerror("memory allocation failure");
		exit(FAIL);
	    }
	    $$->type = obrace;
	    strcpy($$->lexeme, "{");
	    $$->next = NULL;
	    if (!CurKey)
		$$->prev = CurKey, CurKey = $$;
	    else {
		struct Token *p;

		for (p = CurKey; p->next; p = p->next)
		    ;
		p->next = $$, $$->prev = p;
	    }
	}
    anyelements
    '}'
	{
	    if (!($$ = (struct Token *)calloc(1, sizeof(struct Token)))) {
		yyerror("memory allocation failure");
		exit(FAIL);
	    }
	    $$->type = cbrace;
	    strcpy($$->lexeme, "}");
	    $$->next = NULL;
	    if (!CurKey)
		$$->prev = CurKey, CurKey = $$;
	    else {
		struct Token *p;
		
		for (p = CurKey; p->next; p = p->next)
		    ;
		p->next = $$, $$->prev = p;
	    }
	    $$ = CurKey;
	}
    ;

anyelements :
    anyelements
    anyelement
	{
	    $$ = $2;
	}
    |
    anyelement		/* Default action is $$ = $1 */
    ;

anyelement :
    noncommaelement 	/* Default action is $$ = $1 */
    |
    ','
	{
	    if (!($$ = (struct Token *)calloc(1, sizeof(struct Token)))) {
		yyerror("memory allocation failure");
		exit(FAIL);
	    }
	    $$->type = comma;
	    strcpy($$->lexeme, ",");
	    $$->next = NULL;
	    if (!CurKey)
		$$->prev = CurKey, CurKey = $$;
	    else {
		struct Token *p;
		
		for (p = CurKey; p->next; p = p->next)
		    ;
		p->next = $$, $$->prev = p;
	    }
	    $$ = CurKey;
	}
    ;

optwhitespace :
    WHITESPACE
    |
    ;
    
%%
#include "indexlex.c"

main(argc, argv)
int	argc;
char	*argv[];
{
    int			c;
    Boolean		sortfilegiven = FALSE;
    char		sortfilename[81];
    char		indfilename[81];
    struct passwd	*pwentry;
    FILE		*stanzafileptr;

    Whoami = argv[0];
    pwentry = getpwuid(geteuid());
    sprintf(sortfilename, "%s/.alphabetize", pwentry->pw_dir);
    while ((c = getopt(argc, argv, "f:l")) != EOF)
	switch (c) {
	case 'l':
	    Label = TRUE;
	    break;
	case 'f':
	    strcpy(sortfilename, optarg);
	    sortfilegiven = TRUE;
	    break;
	case '?':
	    fprintf(stderr, Usage, Whoami);
	    exit(FAIL);
	}
    stanzafileptr = fopen(sortfilename, "r");
    if (sortfilegiven && !stanzafileptr) {
	fprintf(stderr, "%s: cannot open alphabetization file %s\n", Whoami, sortfilename);
	exit(FAIL);
    }
    if (stanzafileptr) {
	if (!(SortStanza = profile_read_profile(stanzafileptr))) {
	    fprintf(stderr, "%s: file %s is not in stanza format\n", Whoami, sortfilename);
	    fclose(stanzafileptr);
	    exit(FAIL);
	}
	if (!(SortStanza = profile_has_stanza(SortStanza, Marker))) {
	    fprintf(stderr, "%s: file %s does not contain a stanza with marker %s\n", Whoami, sortfilename, Marker[0]);
	    fclose(stanzafileptr);
	    exit(FAIL);
	}
	fclose(stanzafileptr);
    }
    checkstanza(SortStanza);
    if (optind == argc) {
	InputFile = stdin;
	OutputFile = stdout;
	strcpy(IdxFileName, "stdin");
    }
    do {
	if (InputFile != stdin) {
	    strcpy(IdxFileName, argv[optind]);
	    if (!(InputFile = fopen(argv[optind], "r"))) {
		strcpy(IdxFileName, argv[optind]);
		strcat(IdxFileName, ".idx");
		if (!(InputFile = fopen(IdxFileName, "r"))) {
		    fprintf(stderr, "%s: cannot open %s\n", Whoami, IdxFileName);
		    ExitStatus = FAIL;
		    continue;
		}
	    }
	    if (strlen(IdxFileName) >= 4 && eq(&IdxFileName[strlen(IdxFileName)-4], ".idx"))
		sprintf(indfilename, "%.*s.ind", strlen(IdxFileName)-4, IdxFileName);
	    else
		sprintf(indfilename, "%s.ind", IdxFileName);
	    if (!(OutputFile = fopen(indfilename, "w"))) {
		fprintf(stderr, "%s: cannot open output file %s\n", Whoami, indfilename);
		fclose(InputFile);
		ExitStatus = FAIL;
		continue;
	    }
	} else
	    strcpy(IdxFileName, "stdin");
	if (yyparse() != 0)
	    ExitStatus = FAIL;
	fclose(InputFile);
	fclose(OutputFile);
	freetables(IndexTable, IndexTableCount);
	IndexTable = NULL, IndexTableCount = IndexTableSize = 0, LineNo = EndLineNo = 1;
	yysptr = yysbuf;	/* Resets Lex lookahead buffer */
    } while (++optind < argc);
    exit(ExitStatus);
}

yyerror(s)
char *s;
{
    fprintf(stderr, "\"%s\", line %d: %s\n", IdxFileName, LineNo, s);
}

/* Allocates additional space for tables. Returns NULL if memory allocation failure or inconsistent parameters */
char *reallocate(table, current, increment, elementsize)
char	*table;			/* pointer to current table */
int	current;		/* current size of table */
int	increment;		/* additional entries to add */
int	elementsize;		/* size of an element in the table */
{
    char	*calloc(), *realloc();
    char	*p;

    if ((!table && current > 0) || current < 0 || increment < 0 || elementsize < 0)
	return NULL;
    if (increment == 0 || elementsize == 0)
	return table;
    if (current == 0)
	if (!(p = calloc(increment, elementsize)))
	    return NULL;
	else
	    return p;
    else
	if (!(p = realloc(table, (current + increment) * elementsize)))
	    return NULL;
	else
	    return p;
}

/* Frees the space allocated for all the tables */
freetables(index, noentries)
struct IndexEntry	*index;		/* index table */
int			noentries;	/* number of entries in table */
{
    struct Token	*t, *ttemp;
    int			i;

    if (!index || noentries == 0)
	return;
    for (i = 0; i < noentries; i++)
    if (index[i].subitem)
	freetables(index[i].subitem, index[i].subitemcount);	/* recursion! */
    for (t = index[i].tokenlist; t; t = ttemp)
	ttemp = t->next, free(t);
    if (index[i].pagenos)
	free(index[i].pagenos);
    free(index);
}

/* Checks alphabetize stanza for validity */
checkstanza(ps)
PROFILE_STANZA	*ps;
{
    PROFILE_BINDING	*pb;
    PROFILE_VALUE	*pv;
    int			count;

    if (!ps)
	return;
    if (pb = profile_has_binding(ps, "skipchars"))
	for (pv = pb->value; pv; pv = pv->next)
	    if (pv->class != PROFILE_CHARACTER)
		switch (pv->class) {
		case PROFILE_INTEGER:
		    fprintf(stderr, "%s: illegal integer constant %d in skipchars binding\n", Whoami, pv->value.i);
		    break;
		case PROFILE_HEX:
		    fprintf(stderr, "%s: illegal hex constant 0x%x in skipchars binding\n", Whoami, pv->value.i);
		    break;
		case PROFILE_OCTAL:
		    fprintf(stderr, "%s: illegal octal constant 0%o in skipchars binding\n", Whoami, pv->value.i);
		    break;
		case PROFILE_FLOAT:
		    fprintf(stderr, "%s: illegal float constant %f in skipchars binding\n", Whoami, pv->value.f);
		    break;
		case PROFILE_STRING:
		case PROFILE_OTHER:
		    fprintf(stderr, "%s: illegal string constant %s in skipchars binding\n", Whoami, pv->value.s);
		    break;
		}
    if (pb = profile_has_binding(ps, "mapctrlsequence")) {
	for (count = 0, pv = pb->value; pv; pv = pv->next, count++)
	    if (pv->class != PROFILE_OTHER || pv->class != PROFILE_STRING)
		switch (pv->class) {
		case PROFILE_INTEGER:
		    fprintf(stderr, "%s: illegal integer constant %d in mapctrlsequence binding\n", Whoami, pv->value.i);
		    break;
		case PROFILE_HEX:
		    fprintf(stderr, "%s: illegal hex constant 0x%x in mapctrlsequence binding\n", Whoami, pv->value.i);
		    break;
		case PROFILE_OCTAL:
		    fprintf(stderr, "%s: illegal octal constant 0%o in mapctrlsequence binding\n", Whoami, pv->value.i);
		    break;
		case PROFILE_FLOAT:
		    fprintf(stderr, "%s: illegal float constant %f in mapctrlsequence binding\n", Whoami, pv->value.f);
		    break;
		case PROFILE_CHARACTER:
		    fprintf(stderr, "%s: illegal character constant %c in mapctrlsequence binding\n", Whoami, pv->value.c);
		    break;
		}
	if (odd(count))
	    fprintf(stderr, "%s: must have an even number of string values for mapctrlsequence binding\n", Whoami);
    }
    if (pb = profile_has_binding(ps, "mapindexentry")) {
	for (count = 0, pv = pb->value; pv; pv = pv->next, count++)
	    if (pv->class != PROFILE_OTHER || pv->class != PROFILE_STRING)
		switch (pv->class) {
		case PROFILE_INTEGER:
		    fprintf(stderr, "%s: illegal integer constant %d in mapindexentry binding\n", Whoami, pv->value.i);
		    break;
		case PROFILE_HEX:
		    fprintf(stderr, "%s: illegal hex constant 0x%x in mapindexentry binding\n", Whoami, pv->value.i);
		    break;
		case PROFILE_OCTAL:
		    fprintf(stderr, "%s: illegal octal constant 0%o in mapindexentry binding\n", Whoami, pv->value.i);
		    break;
		case PROFILE_FLOAT:
		    fprintf(stderr, "%s: illegal float constant %f in mapindexentry binding\n", Whoami, pv->value.f);
		    break;
		case PROFILE_CHARACTER:
		    fprintf(stderr, "%s: illegal character constant %c in mapindexentry binding\n", Whoami, pv->value.c);
		    break;
		}
	if (odd(count))
	    fprintf(stderr, "%s: must have an even number of string values for mapindexentry binding\n", Whoami);
    }
}

/* Returns the literal string of a token list */
char	*literalstring(t)
struct Token	*t;
{
    static char	literal[81];

    strcpy(literal, "");
    for (t = CurKey; t; t = t->next)
	    strcat(literal, t->lexeme);
    return literal;
}

/* Returns alphabetization string for a token list and a stanza */
char	*alphabetizestring(tokenlist, stanza)
struct Token	*tokenlist;
PROFILE_STANZA	*stanza;
{
    char		litstring[81];
    char		ctrlstring[21];
    char		c[2];
    static char		alphastring[81];
    int			i;
    Boolean		add;
    struct Token	*t;
    PROFILE_BINDING	*pb, *pbchars, *pbctrlsequence;
    PROFILE_VALUE	*pv;
    
    if (!tokenlist)
	return NULL;
    strcpy(alphastring, "");
    if (!stanza) {
	for (t = tokenlist; t; t = t->next)
	    switch (t->type) {
	    case string:
	    case integer:
	    case roman:
	    case comma:
	    case obrace:
	    case cbrace:
		strcat(alphastring, t->lexeme);
		break;
	    }
	return alphastring;
    } else {
	if (pb = profile_has_binding(stanza, "mapindexentry")) {
	    strcpy(litstring, literalstring(tokenlist));
	    for (pv = pb->value; pv && pv->next; pv = pv->next, pv = pv->next)
		if ((pv->class == PROFILE_STRING || pv->class == PROFILE_OTHER) && (pv->next->class == PROFILE_STRING ||
		pv->next->class == PROFILE_OTHER))
		    if (eq(litstring, pv->value.s)) {
			strcpy(alphastring, pv->next->value.s);
			return alphastring;
		    }
	} /* end if there is a mapindexentry binding */
	pbchars = profile_has_binding(stanza, "skipchars");
	pbctrlsequence = profile_has_binding(stanza, "mapctrlsequence");
	c[1] = '\0';
	for (t = tokenlist; t; t = t->next)
	    switch (t->type) {
	    case controlword:
	    case controlsymbol:
		if (pbctrlsequence)
		    for (pv = pbctrlsequence->value; pv && pv->next; pv = pv->next, pv = pv->next)
			if ((pv->class == PROFILE_STRING || pv->class == PROFILE_OTHER) && (pv->next->class == PROFILE_STRING ||
			pv->next->class == PROFILE_OTHER))
			    if (strlen(pv->value.s) > 0) {
				if (pv->value.s[0] != '\\')
				    sprintf(ctrlstring, "\\%s", pv->value.s);
				else
				    strcpy(ctrlstring, pv->value.s);
				if (eq(ctrlstring, t->lexeme))
				    strcat(alphastring, pv->next->value.s);
			    }
		break;
	    case string:
	    case integer:
	    case roman:
		for (i = 0; t->lexeme[i]; i++)
		    if (pbchars) {
			for (add = TRUE, pv = pbchars->value; pv && add; pv = pv->next)
			    if (pv->class == PROFILE_CHARACTER)
				if (pv->value.c == t->lexeme[i])
				    add = FALSE;
			if (add) {
			    c[0] = t->lexeme[i];
			    strcat(alphastring, c);
			}
		    } else {
			c[0] = t->lexeme[i];
			strcat(alphastring, c);
		    }
		break;
	    case comma:
		c[0] = ',';
		goto insert;
	    case obrace:
		c[0] = '{';
		goto insert;
	    case cbrace:
		c[0] = '}';
insert:		if (pbchars) {
		    for (add = TRUE, pv = pbchars->value; pv && add; pv = pv->next)
			if (pv->class == PROFILE_CHARACTER)
			    if (pv->value.c == c[0])
				add = FALSE;
		    if (add)
			strcat(alphastring, c);
		} else
		    strcat(alphastring, c);
		break;
	    }
	return alphastring;
    }
}

/* Finds an entry in a table. Returns NULL if not found. */
struct IndexEntry	*findentry(table, noentries, string)
struct IndexEntry 	*table;
int			noentries;
char			*string;
{
    int	i;

    if (noentries <= 0)
	return NULL;
    for (i = 0; i < noentries; i++)
	if (eq(string, table[i].literal))
	    return &table[i];
    return NULL;
}

/* Returns pointer to page number if found, NULL otherwise */
struct PageNoTable	*findpage(pagearray, elements, pageno, aromannum)
struct PageNoTable	*pagearray;
int			elements;
int			pageno;
Boolean			aromannum;
{
    int	i;

    if (!pagearray)
	return NULL;
    for (i = 0; i < elements; i++)
	if ((pagearray[i].number == pageno) &&
	    (pagearray[i].isroman == aromannum))
	    return &pagearray[i];
    return NULL;
}

/* Sorts the entries in the structures */
sort(base, numberelements)
struct IndexEntry	*base;
int			numberelements;
{
    int	i;
    int numericcompare();
    int alphacompare();

    for (i = 0; i < numberelements; i++) {
	if (base[i].pagenos)
	    qsort(base[i].pagenos, base[i].pagetablecount, sizeof(struct PageNoTable), numericcompare);
	if (base[i].subitem)
	    sort(base[i].subitem, base[i].subitemcount);	/* recursion! */
    }
    qsort(base, numberelements, sizeof(struct IndexEntry), alphacompare);
}

/* Prints out the index entries */
printindexentries(base, noelements, level)
struct IndexEntry	*base;
int			noelements;
int			level;
{
	int             i, j;
	Boolean         prevoutput = FALSE;
	Boolean         prevrange = FALSE;
	char            c;
	char            letter = '\0';

	if (level > ITEMDEPTH)
		return;
	for (i = 0; i < noelements; i++) {
		if (level == 1)
			if (strlen(base[i].alphabetic) > 0)
				if (isalpha(base[i].alphabetic[0])) {
					if (isupper(c = base[i].alphabetic[0]))
						c = tolower(c);
					if (!letter) {
						if (Label) {
							fprintf(OutputFile, "\\indexspace\n");
							fprintf(OutputFile, "\\largeletter{%c}\n", c);
						} else if (prevoutput)
							fprintf(OutputFile, "\\indexspace\n");
					} else if (letter != c) {
						fprintf(OutputFile, "\\indexspace\n");
						if (Label)
							fprintf(OutputFile, "\\largeletter{%c}\n", c);
					}
					letter = c;
				}
		prevoutput = TRUE;
		for (j = 1; j < level; j++)
			fprintf(OutputFile, "  ");
		fprintf(OutputFile, "%s %s ", ItemRep[level - 1], base[i].literal);
		if (base[i].pagenos) {
			for (j = 0; j < base[i].pagetablecount; j++) {
				if (j == base[i].pagetablecount - 1) {
					if (base[i].pagenos[j].isroman == FALSE)
						fprintf(OutputFile, "%d\n", base[i].pagenos[j].number);
					else {
						fprintf(OutputFile, "{\\romannumeral %d}\n", base[i].pagenos[j].number);
					}
				} else if (base[i].pagenos[j].range) {
					if (!prevrange) {
						if (base[i].pagenos[j].isroman == FALSE)
							fprintf(OutputFile, "%d--", base[i].pagenos[j].number);
						else {
							fprintf(OutputFile, "{\\romannumeral %d}--", base[i].pagenos[j].number);
						}
					}
				} else {
					if (base[i].pagenos[j].isroman == FALSE)
						fprintf(OutputFile, "%d, ", base[i].pagenos[j].number);
					else {
						fprintf(OutputFile, "{\\romannumeral %d}, ", base[i].pagenos[j].number);
					}
				}
				prevrange = base[i].pagenos[j].range;
			}
			if (prevrange)
				fprintf(stderr, "%s: file %s, %s %s ends with a range\n", Whoami, IdxFileName, ItemRep[level - 1],
					base[i].literal);
		} else
			fprintf(OutputFile, "\n");
		if (base[i].subitem)
			printindexentries(base[i].subitem, base[i].subitemcount, level + 1);	/* recursion! */
	}
}

int	numericcompare(e1, e2)
struct PageNoTable	*e1, *e2;
{
    if ((e1->isroman == TRUE) && (e2->isroman == FALSE))
        return -1;
    if ((e1->isroman == FALSE) && (e2->isroman == TRUE))
        return 1;
   		/* else either both roman or both integers */
    if (e1->number == e2->number)
	return 0;
    else if (e1->number < e2->number)
	return -1;
    else
	return 1;
}

int alphacompare(e1, e2)
struct IndexEntry	*e1, *e2;
{
    char	s1[81], s2[81];

    strcpy(s1, e1->alphabetic), strcpy(s2, e2->alphabetic);
    return(strcmp(string_downshift(s1), string_downshift(s2)));
}


int rmtoi (romanstr)
char *romanstr;
{
	register char  *p = romanstr;
	register int    w;
	register int    prevw = (-1);
	register int    result = 0;
	int             romanwt();

	while (*p) {
		if ((w = romanwt(*p)) == (-1)) {
			fprintf(stderr, "illegal char in roman string:'%c'\n", (*p));
			return (-1);
		} else {
			if (prevw > 0) {	/* check for subtractive
						 * notation */
				if (w > prevw) { /* e.g., the case "ix" */
					result += (w - prevw) - prevw;
				} else
					result += w;
			} else {
				result += w;
			}
		}
		prevw = w;
		p++;
	}
	return (result);
}

static int romanwt (c)
register char c;
{
	static char     romanlett[7] = {'m', 'd', 'c', 'l', 'x', 'v', 'i'};
	static int      weight[7] = {1000, 500, 100, 50, 10, 5, 1};

	register char  *pt;

	if (isupper(c))
		c = tolower(c);
	pt = romanlett;
	while (*pt) {
		if (*pt == c) {
			return (weight[(int) (pt - romanlett)]);
		} else
			pt++;
	}
	return (-1);		/* roman letter not found */
}

