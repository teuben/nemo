/*
 * HISTORY.C: handle history and headline items on the toplevel 
 *            of a structured binary file.
 *
 *	 9-Mar-88	V1.0	Created					PJT
 *	 1-Jun-88	V1.1	renamed read/write to get/put		PJT
 *	 7-jun-88	V1.2	installed in libT with getparam();
 *                              dynamic HISTORY				PJT
 *      14-Jun-88       V1.3    integrated with IAS NEMO		JEB
 *	15-Dec-88	V1.4	allow creation from scratch		PJT
 *	15-jul-89	V1.5	added reset_history(), fixed order	PJT
 *				of Headline's and History's
 *	29-sep-89	V1.6	added item= keyword                     PJT
 *	 8-dec-89	V1.7	moved some stuff to history.h, no error PJT
 *	10-sep-90       V1.8    << moved history_level around >>	PJT
 *	 9-oct-90	V1.8a   documentation upgraded			PJT
 *      17-jul-91       V1.9    extra option skip=                      PJT
 *
 *	 9-apr-92	V2.0    add_history -> app_history		PJT
 *				because of conflicting libreadline.a
 *	20-feb-94	V2.1    ansi
 *	13-mar-95       V2.2	app_history() now makes new copy of history 
 *				string -- hmm, this bug didn't show up until
 *				linux! - also did the prototypesm
 *
 *	 6-may-95       V2.3    (hisf) allow item=headline         	PJT
 *	21-may-01	V2.4    history <old history_level> back here   PJT
 *       1-aug-02       V2.5    increased memory a bit 256->1024        
 *                              and free up memory in reset_history     PJT
 *      28-dec-04       V2.5a   tag headline properly when free'd       PJT
 *      20-jun-08       V2.5b   history -> nemo_history (MAC library conflict)  PJT
 *
 *  ToDo: not all local variables free up memory, despite that some
 *        have clearly come from allocate'd memory. 
 */

#include <stdinc.h>
#include <strlib.h>
#include <filestruct.h>
#include <history.h>

/* 
 * The datastructure to store history is a simple array of pointers to
 * history strings (in safe memory we may hope).
 * The array is declared static, but you may also think of
 * a dynamically allocated linked list in the future.
 *				-- sure, and just how will that work?  JEB
 */

#ifndef MAXHIST
# define MAXHIST 1024			/* max size of history array        */
#endif

local int nhist = 0;			/* count history data stored so far */
local string histbuf[MAXHIST+1];	/* history string array             */
local string headline = NULL;		/* last headline read in            */
local bool freeup[MAXHIST+1];           /* if space should be free'd        */

int nemo_history = 1;          /* 1=history is auto-mode  0=no history done */
                               /* this item should eventually disappear     */
                               /* it is defined by the user interface,      */
                               /* see getparam.c                            */

/*
 * GET_HISTORY: read history and headline info into local buffer;
 * returns current number of history items
 */

int get_history(stream instr)
{
    for(;;) {                   		/* loop reading input data  */
	if (get_tag_ok(instr, HeadlineTag)) {
	    headline = get_string(instr, HeadlineTag);
	    dprintf(5, "get_history: headline = %s\n", headline);
	} else if (get_tag_ok(instr, HistoryTag)) {
	    if (nhist > MAXHIST) {
		warning("get_history: no more history saved; MAXHIST=%d",
                        MAXHIST);
		return MAXHIST;
	    }
	    histbuf[nhist] = get_string(instr, HistoryTag);
	    dprintf(5, "get_history: histbuf[%d] = %s\n",
		    nhist, histbuf[nhist]);
	    freeup[nhist] = FALSE;
	    nhist++;
	} else
            break;                          /* done with reading loop */
    }
    return nhist;
}

/*
 * PUT_HISTORY:  write current history and headline data to output.
 *               always returns 0 ...
 */

int put_history(stream outstr)
{				
    int i;

    if (!nemo_history) {
    	dprintf(5, "put_history: history data suppressed\n");    
        return 0;
    }
    if (headline != NULL) {
	dprintf(5, "put_history: headline = %s\n", headline);
	put_string(outstr, HeadlineTag, headline);
    }
    dprintf(5, "put_history: writing %d history items\n", nhist);
    for (i = 0; i < nhist; i++) {
	dprintf(5, "             histbuf[%d] = %s\n", i, histbuf[i]);
    	put_string(outstr, HistoryTag, histbuf[i]);
    }
/*  nhist = 0;			// reset counter -- why? JEB */
    return 0;
}

/* 
 * APP_HISTORY: add to history
 *              returns current length of history list
 *
 */

int app_history(string s)
{
    permanent bool warned = FALSE;

    if (nhist > MAXHIST) {
	if (!warned) warning("app_history: too much history");
        warned = TRUE;
	return nhist;
    } else if (s == NULL || streq(s, "")) {
	dprintf(1, "app_history: null history string\n");
	return nhist;
    }
    histbuf[nhist] = scopy(s);
    freeup[nhist] = TRUE;
    dprintf(9,"app_history: histbuf[%d] = %s\n", nhist, s);
    nhist++;
    return nhist;
}

/*
 * RESET_HISTORY: reset the history list, and also free up all memory
 */

void reset_history()
{
  int i;
  for (i=0; i<nhist; i++)
    if (freeup[i]) free(histbuf[i]);
  if (headline) {
    free(headline);
    headline=0;
  }
  nhist = 0;
}

/*
 * SET_HEADLINE: set new value for headline string.
 */

void set_headline(string s)
{
    headline = scopy(s);
}

/*
 * ASK_HEADLINE: enquire about value of headline.
 */

string ask_headline()
{
    return headline;
}

/*
 * ASK_HISTORY: enquire about history data.
 */

string *ask_history()
{
    if (nhist > MAXHIST)
	error("ask_history: too much history");
    histbuf[nhist] = NULL;			/* NULL-term the hist array */
    return histbuf;
}

#ifdef TOOLBOX
			/* PROGRAM: hisf */
#include <getparam.h>

string defv[] = {
    "in=\n		Input file name",
    "out=\n		If given, output history",
    "item=all\n		Show which items",
    "history=\n		If given, add this text to history",
    "skip=1\n           Number to skip before display (0 or 1)",
    "VERSION=2.3\n	6-may-95 PJT",
    NULL,
};

string usage = "list history items from a binary structured file";

string iname,oname;			/* input/output file  name */
stream instr,outstr;

void nemo_main()
{
    int i, skip;
    int item = -1;
    bool Qall = FALSE, Qhead=FALSE;

    iname = getparam("in");
    oname = getparam("out");
    if (streq(getparam("item"),"all")) 
        Qall = TRUE;
    else if (streq(getparam("item"),"headline"))
        Qhead = TRUE;
    else
        item = getiparam("item");       /* get a particular history item */
    skip = getiparam("skip");
    if (skip<0) error("skip=%d cannot be negative",skip);

    if (hasvalue("in")) {
        instr = stropen(iname, "r");                /* open input file */
        get_history(instr);                         /* get old history */
        if (hasvalue("history"))
	    app_history(getparam("history"));       /* add new */
        for (i = skip; i < nhist; i++) 
            if (Qall || i==item) {
                dprintf(1,"%d: ",i);
                printf("%s\n", histbuf[i]);     /* list old history */
            }
        if (Qhead && headline) printf("%s\n",headline);     /* list headline */
        strclose(instr);
    } else if (hasvalue("history")) {
        nhist--;				/* delete first one */
        app_history(getparam("history"));           /* add new from scratch */
    }

    if (hasvalue("out")) {                       /* see if output needed */
        outstr = stropen(oname, "w");               /*  open file */
        put_history(outstr);                        /*  and dump history */
        strclose(outstr);                           
    }
}

#endif
