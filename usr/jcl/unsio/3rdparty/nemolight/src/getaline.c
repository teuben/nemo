/*
 *
 *  GETALINE:  return pointer to next line from a file using dynamic memory
 *  GETSLINE:  --same-- using shared pointer            (not implemented yet)
 * 
From wirzeniu@kruuna.Helsinki.FI Tue Sep 28 17:18:44 1993
Newsgroups: comp.lang.c,alt.sources
Subject: getaline (was: Re: fgets() a text file line)
Followup-To: comp.lang.c
Date: 27 Sep 1993 22:47:29 +0200
Organization: University of Helsinki
NNTP-Posting-Host: kruuna.helsinki.fi
 *
 *	3-nov-93 ported to NEMO, but using K&R notation (no prototype) PJT
 *		 also using NEMO's (re)allocate
 *               idea about getsline()
 * 
 *	20-jun-01	gcc3
 */

#include <stdinc.h>
#include <stdlib.h>

char *getaline(stream f)
{
	char *buf;		/* buffer for line */
	size_t size;		/* size of buffer */
	size_t inc;		/* how much to enlarge buffer */
	size_t len;		/* # of chars stored into buf before '\0' */
	char *p;
	const size_t thres = 128; /* initial buffer size (most lines should
				     fit into this size, so think of this as
				     the "long line threshold").  */
	const size_t mucho = 128; /* if there is at least this much wasted
				     space when the whole buffer has been
				     read, try to reclaim it.  Don't make
				     this too small, else there is too much
				     time wasted trying to reclaim a couple
				     of bytes.  */

	len = 0;
	size = thres;
	buf = (char *) allocate(size);
	if (buf == NULL)
		return NULL;

	while (fgets(buf+len, size-len, f) != NULL) {
		len += strlen(buf+len);
		if (len > 0 && buf[len-1] == '\n')
			break;		/* the whole line has been read */

		for (inc = size, p = NULL; p == NULL && inc > 0; inc /= 2)
			p = (char *) reallocate(buf, size + inc);

		if (inc <= 0 || p == NULL) {
			free(buf);
			return NULL;	/* couldn't get more memory */
		}

		size += inc;
		buf = p;
	}

	if (len == 0) {
		free(buf);
		return NULL;	/* nothing read (eof or error) */
	}

	if (buf[len-1] == '\n')	/* go back on top of the newline */
		--len;
	buf[len] = '\0';	/* unconditionally terminate string,
				   possibly overwriting newline */

	if (size - len > mucho) { /* a plenitude of unused memory? */
		p = (char *) reallocate(buf, len+1);
		if (p != NULL) {
			buf = p;
			size = len+1;
		}
	}

	return buf;
}

/*
 *
 * Example usage:
 *
 *      string s = NULL;
 *      char *cp;
 *
 *      while( getsline(f,&s) != NULL) {
 *          printf("%s\n",s);
 *      }
 *      free(s);
 */
char *getsline(stream f, string *s)
{
    error("Not implemented yet");
    return NULL;
}

#ifdef TESTBED

string defv[] = {
    "in=???\n       input file",
    "VERSION=1\n    pjt",
    NULL,
};

string usage="test getaline";

nemo_main()
{
    char *cp;
    stream f = stropen(getparam("in"),"r");

    while ( (cp=getaline(f)) != NULL )
        printf("%s\n",cp);
}

#endif
