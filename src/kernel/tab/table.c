/* some utility routines for table manipulations */
/*
 * get_line (instr, line)		get a line from a stream and return len
 *					(dangerous: max line len unpecified)
 *
 * parse (linenr, line,dat,ndat)	parse a line into data items
 *
 * strinsert (a, b, n)			insert string 'b' into 'a' for 'n' chars
 *
 * iscomment(line)			is this line a blank or comment line?
 * 
 *    1-jan-04      get_line::  changed EOF to return -1, and empty line to 0
 */
 
#include <stdinc.h>
#include <ctype.h>

#if !defined(HUGE)
#define HUGE 1e20
#endif

#ifndef MAX_LINELEN
#define MAX_LINELEN  16384
#endif

/*
 * insert a string 'b' into 'a' replacing the first 'n' positions into 'a'
 *      (see als hsh.h in ..../hermes/lib  --  PJT)
 */

void strinsert(char *a, char *b, int n)
{
	int    i, idiff, alen, blen;

	alen = strlen(a);
	blen = strlen(b);
	idiff = blen - n;

	dprintf(8,"strinsert  in: (%s)-(%s)-%d\n", a, b, n);
	if (idiff < 0)		/* shift a to the left */
		for (i = n; i <= alen; i++)
			a[i + idiff] = a[i];
	else if (idiff > 0)	/* shift a to the right */
		for (i = alen; i >= n; i--)
			a[i + idiff] = a[i];

	for (i = 0; i < blen; i++)
		a[i] = b[i];
	dprintf(8,"strinsert out: (%s)-(%s)-%d\n", a, b, n);
}


int get_line(stream instr, char *line)
{
	int  c, i=0;
	
	for(;;) {
		c=getc(instr);
		if (c==EOF)
			return -1;	/* error code: EOF */
		else if (c=='\n')
			break;
		line[i++]=c;
		if (i>MAX_LINELEN) {
			warning("get_line: max linelen (%d) exceeded; return 0",MAX_LINELEN);
			return 0;
		}
	}
	line[i]=0;
	return strlen(line);
}

#define PARESC  '%'

	/* linenr: current line number for %0 parameter */
	/* line:   ascii string */
	/* dat:    space where to store parsed items */
	/* ndat:   maximum that can be filled */

void parse(int linenr, char *line, double *dat, int ndat)
{
    int i, k, idat;
    char number[100];

    i=0;			/* position in line */
    while (line[i] != '\0') {
        if (line[i]==PARESC) {  /* referenced a column ? */
            k = i+1;            /* point to start of column index */
            while (isdigit(line[k]))
                k++;
            dprintf (2,"%s PARESC: i=%d k=%d\n",line,i,k);
            k -= i;                 /* number of digits + 1 of column index */
            idat = atoi(&line[i+1]);     /* index to column in dat[] */
            if (idat==0)
                sprintf (number,"(%d)",linenr);
            else if (idat <= ndat)
                sprintf (number,"(%16.10e)",dat[idat-1]); 
            else
                error("tabmath: referencing column %d\n",idat);
            strinsert(&line[i],number,k);
            i += strlen(number);
        } else
            i++;
    }
}


/*
 * iscomment: is a line a comment or blank line ?
 *		comment lines must start with '#', ';' or '!'
 *		or a NULL or newline '\n'
 */

int iscomment(char *line)
{
    char *cp;

    /* first do a quick check on the first character */    
    if (*line=='#' || *line==';' || *line=='!' || *line=='\n' || *line=='\0')
        return 1;
    /* then loop over the line until first non-whitespace is seen */
    for (cp=line; *cp; cp++)
        if (!isspace(*cp)) return 0;
    /* if none found, it bumped into the end of line and it's a comment */
    return 1;
}


/*
 * parse_select, with dynamic reallocation
 *
 */

#if 0

void parse_select (select_str, select, nselect, Qsel, nQsel)
char *select_str;
int  **select, *nselect, nQsel;
bool **Qsel;
{
    error ("parse_select: Not implemented yet");
}

#endif

