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
#include <table.h>

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
#if 1	
	static int _first = TRUE;

	if (_first) {
	  warning("old get_line() is still used - it is being deprecated");
	  _first = FALSE;
	}
#endif
	
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

table *table_open(stream instr, int mode)
{
  tableptr tptr = (tableptr) allocate(sizeof(table));

  tptr->str   = instr;
  tptr->mode  = mode;
  tptr->lines = NULL;


  // in full buffering mode, the whole file is read into memory
  // using the *tptr->lines (or linked list?)

  return tptr;
}

void table_close(tableptr tptr)
{
  
}

ssize_t table_line(tableptr tptr, char **line, size_t *linelen)
{
  // in simple (streaming) mode, just get the next line
  if (tptr->mode == 0)
    // return getdelim(line, linelen, ' ', tptr->str);
    return getline(line, linelen, tptr->str);

  
  error("Unsupported table mode=%d", tptr->mode);
  return -1;
}




#ifdef TESTBED

#include <getparam.h>

string defv[] = {
    "file=???\n           Input or Output file",
    "mode=r\n             Read (r) or Write (w)",
    "VERSION=2.0\n        10-feb-2022 PJT",
    NULL,
};

string usage = "testing tables";


void nemo_main()
{
    tableptr tp1;
    stream instr, outstr;
#if 1
    size_t linelen = 0;                  // getline() is allowed to start from 0
    char *line = NULL;
#else
    size_t linelen = MAX_LINELEN;        // does we need to worry about the extra newline
    char *line = allocate(linelen);      // allocate formally has the wrong argument type
    //char *line = malloc(linelen);
#endif

    if (strcmp(getparam("mode"),"w") == 0) {
      dprintf(0,"write mode\n");
      
      outstr = stropen(getparam("file"),"w");

#if 0
      // special test to get an embedded 0
      char *buffer = "A\n\0\nHello\n";
      fwrite(buffer,10,1,outstr);
#endif
      
      strclose(outstr);
      
    } else {
      dprintf(0,"read mode\n");
      
      instr = stropen(getparam("file"),"r");
      dprintf(0,"linelen=%ld\n", linelen);
#if 0
      // original C
      while (getline(&line, &linelen, instr) >= 0)
	printf("line[%ld] = %s\n", linelen, line);
#else
      // new style table2 
      tp1 = table_open(instr, 0);
      while (table_line(tp1, &line, &linelen) >= 0)
	printf("line[%ld] = %s\n", linelen, line);
      table_close(tp1);
#endif
      
      free(line);
    }

}


#endif
