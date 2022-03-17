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

bool ispipe(stream instr)
{
  off_t try = lseek(fileno(instr), 0, SEEK_CUR);
  if (try < 0) return TRUE;
  return FALSE;
}

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
  tptr->nr    = 0;
  tptr->nc    = 0;

#if 0  
  tptr->linelen = 0;
  tptr->line    = NULL;
#else
  tptr->linelen = 20;
  tptr->line    = malloc(tptr->linelen);
#endif
  dprintf(0,"table_open - got %d chars allocated at the start\n", tptr->linelen);

  // in full buffering mode, the whole file is read into memory
  // using the *tptr->lines (linked list?)

  if (mode == 1) {
    //
    warning("new mode=1");
    // we could cheat for now and find the number of lines 
    // allocate tptr->lines and stuff them there
    // for flexibility this should become a linked list
  }

  return tptr;
}

struct LineNode {
  char *val; // strdup
  struct LineNode* next;
};


// non-seekable file mode=1
table *table_open0(stream instr, int mode)
{
  tableptr tptr = (tableptr) allocate(sizeof(table));

  tptr->str   = instr;
  tptr->mode  = mode;
  tptr->lines = NULL;
  tptr->nr    = 0;
  tptr->nc    = 0;

#if 0  
  tptr->linelen = 0;
  tptr->line    = NULL;
#else
  tptr->linelen = 20;
  tptr->line    = malloc(tptr->linelen);
#endif
  dprintf(0,"table_open - got %d chars allocated at the start\n", tptr->linelen);


  // we could cheat for now and find the number of lines 
  // allocate tptr->lines and stuff them there
  // for flexibility this should become a linked list

  struct LineNode* first = (struct LineNode*) allocate(sizeof(struct LineNode));
  struct LineNode* curr = first;

  int nLines = 0;
  string line;

  while (1) {
    line = table_line0(tptr);

    if (line == NULL)
      break;
      
    line = strdup(line);
    nLines++;
    curr->val = line;
    curr->next = (struct LineNode*) allocate(sizeof(struct LineNode));
    curr = curr->next;
    curr->next = NULL;
  }
  tptr->nr = nLines;

  return tptr;
}

// seekable file mode=1
table *table_open1(stream instr, int mode, int nlines)
{
  tableptr tptr = (tableptr) allocate(sizeof(table));

  tptr->str   = instr;
  tptr->mode  = mode;
  tptr->lines = NULL;
  tptr->nr    = 0;
  tptr->nc    = 0;

#if 0  
  tptr->linelen = 0;
  tptr->line    = NULL;
#else
  tptr->linelen = 20;
  tptr->line    = malloc(tptr->linelen);
#endif
  dprintf(0,"table_open - got %d chars allocated at the start\n", tptr->linelen);


  warning("table_open1: mode=%d nlines=%d",mode,nlines);
  tptr->nr = nlines;
  tptr->lines = (string *) allocate(nlines*sizeof(string));
  for (int i=0; i<nlines; i++)
    tptr->lines[i] = strdup(table_line0(tptr));

  return tptr;
}

size_t table_nrows(tableptr tptr)
{
  return tptr->nr;
}


size_t table_ncols(tableptr tptr)
{
  return tptr->nc;
}


void table_reset(tableptr tptr)
{
  // reset the table 
}


void table_close(tableptr tptr)
{
  // free that memory
}

// deprecate
ssize_t table_line1(tableptr tptr, char **line, size_t *linelen)
{
  // in simple (streaming) mode, just get the next line
  if (tptr->mode == 0)
    // return getdelim(line, linelen, ' ', tptr->str);
    return getline(line, linelen, tptr->str);

  error("Unsupported table mode=%d", tptr->mode);
  return -1;
}

string table_line0(tableptr tptr)
{
  // in simple (streaming) mode, just get the next line
  if (tptr->mode == 0) {
    // return getdelim(line, linelen, ' ', tptr->str);
    ssize_t ret = getline(&(tptr->line), &(tptr->linelen), tptr->str);
    if (ret >= 0)
      return tptr->line;
    return NULL;
  }

  error("Unsupported table mode=%d", tptr->mode);
  return NULL;
}

string table_row(tableptr tptr, int row)
{
  return tptr->lines[row];
}

#ifdef TESTBED

#include <getparam.h>

string defv[] = {
    "file=???\n           Input or Output file",
    "mode=r\n             Read (r) or Write (w)",
    "test=0\n             testmode",
    "VERSION=2.0\n        10-mar-2022 PJT",
    NULL,
};

string usage = "testing tables";

void testmode1();
void testmode2();

void nemo_main()
{
    int testmode = getiparam("test");
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

    if (testmode == 1) {
      testmode1();
      return;
    }

    if (testmode == 2) {
      testmode2();
      return;
    }

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

#if 0
      // using your own allocation
      while (table_line1(tp1, &line, &linelen) >= 0)
	printf("line[%ld] = %s", linelen, line);
#else
      // depending on table internals
      while ( table_line0(tp1) )
	printf("line[%ld] = %s", tp1->linelen, tp1->line);
#endif
      
      table_close(tp1);

      
#endif
      
      free(line);
    }

}


#endif


void testmode1()
{
  string input = getparam("file");
  int nlines = nemo_file_lines(input,0);
  stream instr = stropen(input,"r");
  tableptr tp1 = table_open1(instr, 0, nlines);     // read the whole file in memory

  dprintf(0,"nlines: %d\n",tp1->nr);

  printf("first line: %s",table_row(tp1,0));
  printf("last  line: %s",table_row(tp1,tp1->nr - 1));
	 
}

void testmode2()
{
  string input = getparam("file");
  stream instr = stropen(input,"r");
  tableptr tp1 = table_open0(instr, 0);     // read the whole file in memory - linked list

  dprintf(0,"nlines: %d\n",tp1->nr);

  // printf("first line: %s",table_row(tp1,0));
  // printf("last  line: %s",table_row(tp1,tp1->nr - 1));
	 
}
