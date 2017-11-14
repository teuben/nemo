/*  burststring, burst2string, freestrings, splitstring(not tested!!!)
 *
 *  25-dec-94 added burst2string()
 *  20-may-96 added burst0string() for checker
 *  28-nov-00 changed some names to make it compile for C++
 *  20-jun-01 prototypes to help gcc3
 *
 * BURSTSTRING: break a string of the form "word1, word2, ..." into
 * seperate strings "word1", "word2", ... and return them in an
 * extended-string (ie, NULL-terminated sequence of pointers).
 */

#include <stdinc.h>
#include <extstring.h>

#define MWRD 2048	/* max words in list */
#define MSTR 256	/* max chars per word */

string *burststring(string lst, string sep)
{
    string wrdbuf[MWRD], *wp;
    char strbuf[MSTR], *sp, *lp;

    wp = wrdbuf;
    sp = strbuf;
    lp = lst;
    do {						/* scan over list */
	if (*lp == 0 || strchr(sep, *lp) != NULL) {	/*   is this a sep? */
	    if (sp > strbuf) {				/*     got a word? */
		*sp = 0;
		*wp++ = (string) copxstr(strbuf, sizeof(char));
		if (wp == &wrdbuf[MWRD])		/*       no room? */
		    error("burststring: too many words (%d)",MWRD);
		sp = strbuf;				/*       for next 1 */
	    }
	} else {					/*   part of word */
	    *sp++ = *lp;				/*     so copy it */
	    if (sp == &strbuf[MSTR])			/*     no room? */
		error("burststring: word too long (%d)",MSTR);
	}
    } while (*lp++ != 0);				/* until list ends */
    *wp = NULL;
    return ((string *) copxstr((char *)wrdbuf, sizeof(string)));	/*PPAP*/
}

/*
 * BURST1STRING: similar to burststring, but now each separator
 * triggers a new word for output so you can have zero-length
 * words coming out. Use with care.
 */

string *burst0string(string lst, string sep)
{
    string wrdbuf[MWRD], *wp;
    char strbuf[MSTR], *sp, *lp;

    wp = wrdbuf;
    sp = strbuf;
    lp = lst;
    do {						/* scan over list */
	if (*lp == 0 || strchr(sep, *lp) != NULL) {	/*   is this a sep? */

		*sp = 0;
		*wp++ = (string) copxstr(strbuf, sizeof(char));
		if (wp == &wrdbuf[MWRD])		/*       no room? */
		    error("burststring: too many words\n");
		sp = strbuf;				/*       for next 1 */

	} else {					/*   part of word */
	    *sp++ = *lp;				/*     so copy it */
	    if (sp == &strbuf[MSTR])			/*     no room? */
		error("burststring: word too long\n");
	}
    } while (*lp++ != 0);				/* until list ends */
    *wp = NULL;
    return ((string *) copxstr((char *)wrdbuf, sizeof(string)));	/*PPAP*/
}

/*
 * BURST2STRING: like burstring() but also copies the separaters
 * in the alternating parts of the extneded string.
 *  To find out if the first 'word' is a word, or a separator 'word'
 *  you can try strchr(sep,*b2s[0]): should return NULL is the first
 *  word is a true word.
 *  This routine is useful if you need to paste the 'words' back together,
 *  perhaps with some editing, into the original string (lst)
 *
 */

string *burst2string(string lst, string sep)
{
    string wrdbuf[MWRD], *wp;
    char strbuf[MSTR], *sp, *lp;
    int olds, news;

    wp = wrdbuf;
    sp = strbuf;
    lp = lst;
    if (*lp)
    	olds = (strchr(sep,*lp) != NULL);
    else
        olds = 0;    /* doesn't matter what this is */
    do {						/* scan over list */
        if (*lp)
            news = strchr(sep,*lp) != NULL;
        else
            news = !olds; /* force update */

        if (news != olds) {
            *sp = 0;
            *wp++ = (string) copxstr(strbuf, sizeof(char));
            if (wp == &wrdbuf[MWRD])		/*       no room? */
	        error("burst2string: too many words");
	    sp = strbuf;				/*       for next 1 */
	}
        *sp++ = *lp;				/*     so copy it */
	if (sp == &strbuf[MSTR])			/*     no room? */
	    error("burst2string: word too long");
        olds = news;
    } while (*lp++ != 0);				/* until list ends */
    *wp = NULL;
    return ((string *) copxstr((char *)wrdbuf, sizeof(string)));	/*PPAP*/
}


void freestrings(              /* free a previously allocated */
	string *strptr)                  /* array of strings */
{
    string *s;

    for (s=strptr; *s != NULL; s++)             /* loop over all strings */
        free(*s);                               /* free each one */
    free(strptr);                               /* and free the 'parent' */
}


int splitstring(
int maxout,		/* declared length of out[] array */
string out[],		/* array for pointers to words that remain in lst */
string lst,		/* list of words to patch and seperate */
string sep)		/* chars which seperate words */
{
    string wrdbuf[MWRD], *wp;
    char strbuf[MSTR], *sp, *lp;

    error("splitstring: not implemented yet");
    wp = out;	/* point to first word */
    sp = lst;	/* point to first char */
    lp = lst;
    do {						/* scan over list */
	if (*lp == 0 || strchr(sep, *lp) != NULL) {	/*   is this a sep? */
	    if (sp > strbuf) {				/*     got a word? */
		*sp = 0;
		*wp++ = (string) copxstr(strbuf, sizeof(char));
		if (wp == &wrdbuf[MWRD])		/*       no room? */
		    error("splitstring: too many words\n");
		sp = strbuf;				/*       for next 1 */
	    }
	} 
    } while (*lp++ != 0);				/* until list ends */

    return wp-out;
}


#ifdef TESTBED

#include <getparam.h>

string defv[] = {
    "lst=foo, bar,waldo ,",
    "sep= ,",
    NULL,
};

nemo_main()
{
    string lst, sep, *wrds;

    lst = getparam("lst");
    sep = getparam("sep");
    printf("lst=%s\n",lst);
    printf("sep=%s\n",sep);

    wrds = burststring(lst, sep);
    printf("(%d items): ",xstrlen(wrds,sizeof(string))-1);
    while (*wrds != NULL)
	printf("\"%s\"  ", *wrds++);
    printf("\n");
    freestrings(wrds);

    wrds = burst0string(lst, sep);
    printf("(%d items): ",xstrlen(wrds,sizeof(string))-1);
    while (*wrds != NULL)
	printf("\"%s\"  ", *wrds++);
    printf("\n");
    freestrings(wrds);
    

    
}
#endif
