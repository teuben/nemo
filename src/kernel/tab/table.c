/* some utility routines for table manipulations.
 *
 * New table V2 I/O routines: - test2 -
 *
 *    table_open
 *    table_line
 *    table_line1
 *    table_close
 *
 *
 * ==== some deprecatable functions ====
 *
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
#include <extstring.h>
#include <mdarray.h>

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
 *		comment lines must start with '#' ';' '!' '/'
 *		or a NULL or newline '\n' (blank line)
 */

int iscomment(char *line)
{
    char *cp;

    /* first do a quick check on the first character */    
    if (*line=='#' || *line==';' || *line=='!' || *line=='/' || *line=='\n' || *line=='\0')
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


struct LineNode {        // simple linked list while reading lines
  char *val;
  struct LineNode* next;
};


table *table_open(stream instr, int mode)
{
  tableptr tptr = (tableptr) allocate(sizeof(table));
  dprintf(1,"table_open: mode=%d\n",mode);

  tptr->str     = instr;
  tptr->mode    = mode;
  tptr->lines   = NULL;
  tptr->nr      = 0;
  tptr->nc      = 0;
  tptr->linelen = 0;
  tptr->line    = NULL;
  dprintf(1,"table_open - got %d chars allocated at the start\n", tptr->linelen);

  if (mode <= 0) {   //  read table in memory, also separate header (comments) from body of table
    // note:    mode<0 treats all lines the same
    //          mode=0 should split comments out @todo
    dprintf(1,"linked list reading of table\n");

    // use a linked list to read all lines
    struct LineNode* first = (struct LineNode*) allocate(sizeof(struct LineNode));
    struct LineNode* curr = first;

    int nLines = 0;
    string line;

    while (1) {
      line = table_line(tptr);
      if (line == NULL)
	break;
      if (iscomment(line)) {
	// for now, skip them
	continue;
      }
      nLines++;
      curr->val = strdup(line);
      curr->next = (struct LineNode*) allocate(sizeof(struct LineNode));
      curr = curr->next;
      curr->next = NULL;
    }
    // allocate the array of nLines strings, point each line there
    tptr->nr = nLines;
    tptr->lines = (string *) allocate(nLines*sizeof(string));
    curr = first;
    for (int i=0; i<nLines; i++) {
      tptr->lines[i] = curr->val;
      curr = curr->next;
    }
    // @todo free up space of the linked list itself
  }
  
  // done!

  dprintf(1,"Read %d lines so far\n",tptr->nr);
  return tptr;
}

// seekable file mode=1 where you know number of lines
table *table_open1(stream instr, int mode, int nlines)
{
  tableptr tptr = (tableptr) allocate(sizeof(table));

  tptr->str   = instr;
  tptr->mode  = mode;
  tptr->lines = NULL;
  tptr->nr    = 0;
  tptr->nc    = 0;

  tptr->linelen = 0;
  tptr->line    = NULL;
  dprintf(0,"table_open - got %d chars allocated at the start\n", tptr->linelen);

  if (mode == 1)
    return tptr;

  // @todo what if nlines was not enough.
  // @todo skip comments here too?
  warning("table_open1: mode=%d nlines=%d",mode,nlines);
  tptr->nr = nlines;
  tptr->lines = (string *) allocate(nlines*sizeof(string));
  for (int i=0; i<nlines; i++)
    tptr->lines[i] = strdup(table_line(tptr));

  return tptr;
}

size_t table_nrows(tableptr tptr)
{
  return tptr->nr;
}


size_t table_ncols(tableptr tptr)
{
  // if tptr->nr > 0 and nc==0, force to read a line in mode=0 and set nc
  if (tptr->nr > 0 && tptr->nc == 0) {
    if (tptr->mode == 0) {
      string *sp = table_rowsp(tptr,0);
      dprintf(1,"table_ncols: processed first line to get nc -> %d\n",tptr->nc);
      tptr->nc = xstrlen(sp,sizeof(string)) - 1;
      // free
    } else {
      warning("mode=1 ... does not have ncols yet");
    }
  }
  return tptr->nc;
}


void table_reset(tableptr tptr)
{
  // reset the table   --    TBD
}


void table_close(tableptr tptr)
{
  // free that memory
  free(tptr->line);
  tptr->linelen = 0;
  // @todo - free more
}

string table_line(tableptr tptr)
{
  ssize_t ret = getline(&(tptr->line), &(tptr->linelen), tptr->str);
  if (ret >= 0) {
    if (tptr->line[ret-1] == '\n')
      tptr->line[ret-1] = '\0';
    return tptr->line;
  }
  // end of file
  return NULL;
}

table *table_cat(table* t1, table* t2, int mode)
{
  tableptr tptr = (tableptr) allocate(sizeof(table));
  
  tptr->lines = NULL;
  tptr->nr = 0;
  tptr->nc = 0;

#if 0
    tptr->linelen = 0;
    tptr->line = NULL;
#else
    tptr->linelen = 20;
    tptr->line = malloc(tptr->linelen);
#endif 
  dprintf(0,"table_cat - got %d chars allocated at the start\n", tptr->linelen);

  //concatenation(above and below)
  if(mode == 0) {
    //updating # of rows in the new table
    int totalLines = t1->nr + t2->nr;
    tptr->nr = totalLines;

    //updating # of columns to be the larger of the original two tables
    if(t1->nc > t2->nc){
      tptr->nc = t1->nc;
    } else {
      tptr->nc = t2->nc;
    }

    tptr->lines = (string*) allocate((tptr->nr) * sizeof(string));
    for(int i = 0; i < t1->nr; i++){
      tptr->lines[i] = strdup(table_line(t1));
    }
    for(int i = 0; i < t2->nr; i++){
      tptr->lines[i + t1->nr] = strdup(table_line(t2));
    }

  } else { //paste mode
    //assume number of rows are the same
    tptr->nr = t1->nr;
    //How do I update #of column in this case
    tptr->nc = t1->nc + t2->nc;
    tptr->lines = (string*) allocate((tptr->nr) * sizeof(string));

    for(int i = 0; i < t1->nr; i++){
      //@TODO: add white spaces between two lines
      tptr->lines[i] = strcat(strdup(table_line(t1)), strdup(table_line(t2)));
    }

  }

  return tptr;
} 


ssize_t table_line1(tableptr tptr, char **line, size_t *linelen, int newline)
{
  // in simple (streaming) mode, just get the next line
  if (tptr->mode == 1) {
    ssize_t len1 = getline(line, linelen, tptr->str);
    if (newline) return len1;
    if (*line[len1-1] == '\n') {
      line[len1-1] = '\0';
      len1--;
      return len1;
    }
  }

  error("line1: Unsupported table mode=%d", tptr->mode);
  return -1;
}


string table_row(tableptr tptr, int row)
{
  return tptr->lines[row];
}

// linked list placeholder
typedef struct lls {
  char *val;
  struct lls *next;
} lls;

// return list of zero terminated (extstring) pointers to the words in a row
string *table_rowsp(table *t, int row)
{
  char *line = strdup(t->lines[row]);
  int ntok = 0;
  char *token = strtok(line," ,");
  lls *first = (lls *) allocate(sizeof(lls));
  lls *curr = first;
  string *sp;
  int i;

  // find tokens
  while (token != NULL) {
    dprintf(2,"token: %d  %s\n",ntok, token);
    curr->val = token;
    curr->next = (lls *) allocate(sizeof(lls));
    curr = curr->next;
    curr->next = NULL;

    ++ntok;
    token = strtok(NULL, " ,");
  }

  // allocate the array to hold the tokens
  sp = (string *) allocate((ntok+1) * sizeof(string));
  curr = first;
  for (i=0; i<ntok; i++) {
    sp[i] = curr->val;
    curr = curr->next;
  }
  sp[ntok] = NULL;

  // consistency check on # columns
  if (t->nc == 0) {
    dprintf(1,"table_rowsp[%d]: setting ncols=%d\n",row,ntok);
    t->nc = ntok;
  } else if (ntok != t->nc)
    error("column number changed:  %d -> %d\n",t->nc, ntok);
  
  // free memory of linked list
  curr = first;
  while (curr->next) {
    dprintf(2,"free up LLS with %s\n",curr->val);
    first = curr->next;
    free(curr);
    curr = first;
  }
  // list of strings
  return sp;
}

/*
 *    a_i_j   row-major if a11 and a12 follow in memory    C,python:    data[0][0], data[0][1], ...
 *            col-major if a11 and a21 follow in memory    Fortran:     data(1,1), data(2,1)
 */

//
// nrow=0 and/or ncol=0 are allowed. It simply means all rows/cols are
// used in order
// Returns data[row][col]
//
mdarray2 table_md2rc(table *t, int nrow, int *rows, int ncol, int *cols)
{
  int nr = table_nrows(t);
  int nc = table_ncols(t);
  dprintf(1,"table_md2rc: table %d x %d \n",nr,nc);
  dprintf(1,"table_md2rc: data2 ncol=%d nrow=%d\n",ncol,nrow);
  mdarray2 a = allocate_mdarray2(nr,nc);    // a[nr][nc]
  int i,j;

  for (i=0; i<nr; i++) {
    string s = table_row(t,i);
    dprintf(1,"%d: %s\n",i,s);
    string *sp = table_rowsp(t,i);
    for (j=0; j<nc; j++) {
      a[i][j] = atof(sp[j]);
    }
  }

  return a;
}

// return a data[col][row] based table
// this is the more practical
// note column 0 has a special meaning, it's the row number
mdarray2 table_md2cr(table *t, int ncol, int *cols, int nrow, int *rows)
{
  int i,j,jidx;
  int nr = table_nrows(t);
  int nc = table_ncols(t);
  dprintf(1,"table_md2cr: table %d x %d \n",nr,nc);
  dprintf(1,"table_md2cr: data2 ncol=%d nrow=%d\n",ncol,nrow);
  if (ncol > 0) {
    for (j=0; j<ncol; j++) {
      jidx = cols[j];
      if (jidx < 0)  error("illegal column reference %d < 0", jidx);
      if (jidx > nc) error("illegal column reference %d > %d", jidx,nc);
    }
  }
  
  if (ncol>0) nc=ncol;
  if (nrow>0) nr=nrow;   // not supported yet  
  mdarray2 a = allocate_mdarray2(nc,nr);                  // a[nc][nr]

  for (i=0; i<nr; i++) {
    string s = table_row(t,i);
    dprintf(1,"%d: %s\n",i,s);
    string *sp = table_rowsp(t,i);
    for (j=0; j<nc; j++) {
      jidx = (ncol == 0 ?  j  :  cols[j]-1);
      if (jidx < 0)
	a[j][i] = i+1;
      else
	a[j][i] = atof(sp[jidx]);
    }
  }

  return a;
}


#ifdef TESTBED

#include <getparam.h>

string defv[] = {
    "file=???\n           Input or Output file",
    "write=f\n            Write mode? or else Write (w)",
    "test=0\n             testmode",
    "mode=0\n             buffering mode for table_open",
    "VERSION=2.0\n        29-apr-2022 PJT",
    NULL,
};

string usage = "testing tables";

void testmode1();
void testmode2();
void testmode3();
void testmode4();
void testmode5();

void nemo_main()
{
    int testmode = getiparam("test");
    int mode = getiparam("mode");
    tableptr tp1;
    stream instr, outstr;

    if (testmode == 1) return testmode1();
    if (testmode == 2) return testmode2();
    if (testmode == 3) return testmode3();
    if (testmode == 4) return testmode4();
    if (testmode == 5) return testmode5();

    warning("Falling through testmodes");

    if (getbparam("write")) {
      dprintf(0,"write mode\n");
      outstr = stropen(getparam("file"),"w");
      strclose(outstr);
      return;
    } 


    dprintf(0,"read mode\n");
    instr = stropen(getparam("file"),"r");
#if 0
    // original C
    size_t linelen = 0;                  // getline() is allowed to start from 0
    char *line = NULL;
    while (getline(&line, &linelen, instr) >= 0)
      printf("line[%ld] = %s\n", linelen, line);
    free(line);
#else
    // new style table2 : mode picks the reading by line, or using the
    //                    slightly slower linked-list reader
    char *cp;
    int nl = 0;
    tp1 = table_open(instr, mode);
    while ( (cp=table_line(tp1)) )
      nl++;
    table_close(tp1);
    dprintf(0,"Read %d extra in mode %s\n",nl,mode);
#endif

}



// linear list, but not supported in a pipe
void testmode1()
{
  string input = getparam("file");
  int nlines = nemo_file_lines(input,0);            // this will not work on pipes
  dprintf(0,"nemo_file_lines: %d\n",nlines);
  stream instr = stropen(input,"r");
  if (ispipe(instr))
    error("This mode is not supported with a pipe");
  rewind(instr);
    
  tableptr tp1 = table_open1(instr, 0, nlines);     // read the whole file in memory
  dprintf(0,"found nlines: %d\n",tp1->nr);
  printf("last  line: %s\n",table_row(tp1,tp1->nr - 1));
  printf("first line: %s\n",table_row(tp1,0));
}

// streaming
void testmode3()
{
  string input = getparam("file");
  string s;
  int nlines = nemo_file_lines(input,0);            // this will not work on pipes
  dprintf(0,"nemo_file_lines: %d\n",nlines);
  stream instr = stropen(input,"r");
  if (ispipe(instr))
    error("This mode is not supported with a pipe");
  rewind(instr);
    
  tableptr tp1 = table_open1(instr, 1, nlines);     // read the whole file in memory
  int nl=0;
  while ((s=table_line(tp1)))
    nl++;
    
  dprintf(0,"found nlines: %d\n",nl);
  // note cannot use table_row() since lines have not been saved

}

// linked list read, supports pipes
void testmode2()
{
  string input = getparam("file");
  stream instr = stropen(input,"r");
  tableptr tp1 = table_open(instr, 0);     // read the whole file in memory - linked list

  dprintf(0,"nlines: %d\n",tp1->nr);

  // assuming the linked list has been converted....
  printf("first line: %s\n",table_row(tp1,0));
  printf("last  line: %s\n",table_row(tp1,tp1->nr - 1));
	 
}

// linked list, streaming one line, and not keeping a buffer
void testmode4()
{
  string input = getparam("file");
  string s;
  stream instr = stropen(input,"r");
    
  tableptr tp1 = table_open(instr, 1);
  int nl=0;
  while ((s=table_line(tp1)))
    nl++;
    
  dprintf(0,"found nlines: %d\n",nl);
  // note cannot use table_row() since lines have not been saved
  dprintf(0,"has buffered nlines: %d\n",tp1->nr);
}

void testmode5()
{
  string input = getparam("file");
  stream instr = stropen(input,"r");
  tableptr tp1 = table_open(instr, 0);
  
#if 1
  int nl = table_nrows(tp1);
  int nc = table_ncols(tp1);
  dprintf(0,"testmode5: table %d x %d\n",nl,nc);

  for (int j=0; j<nl; j++) {
    string *sp = table_rowsp(tp1,j);
    int nsp = xstrlen(sp, sizeof(string))-1;
    if (nc == 0)
      nc = nsp;
    else if (nsp != nc)
      warning("Number of columns changed in row %d.... %d -> %d",j+1,nc,nsp);
    dprintf(1,"%d",j);
    for (int i=0; i<nsp; i++)
      dprintf(1," %s",sp[i]);
    dprintf(1,"\n");
    // free this memory is painful, so table should keep it.....
    char *cp = sp[0];
    free(cp);
  }
#endif
  
  mdarray2 a = table_md2rc(tp1,0,0,0,0);    //   a[row][col]
  
  dprintf(0,"after table_md2rc: table %d x %d\n",tp1->nr,tp1->nc);  
  for (int j=0; j<tp1->nr; j++) {
    dprintf(0,"%d",j);      
    for (int i=0; i<tp1->nc; i++)
      dprintf(0," %g",a[j][i]);
    dprintf(0,"\n");
  }

  mdarray2 b = table_md2cr(tp1,0,0,0,0);    //   a[col][row]
  
  dprintf(0,"after table_md2cr: table %d x %d\n",tp1->nr,tp1->nc);  
  for (int j=0; j<tp1->nr; j++) {
    dprintf(0,"%d",j);      
    for (int i=0; i<tp1->nc; i++)
      dprintf(0," %g",b[i][j]);
    dprintf(0,"\n");
  }
  
    
}

#endif
