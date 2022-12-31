/*
 *  various support for table I/O
 *
 *  feb-2022     Table I/O - V2.0
 *
 *  Additional support is given via burststring.c and extstring.c
 *  Deprecation messages added to old routine
 */

#include <mdarray.h>

#ifndef _h_table
#define _h_table

//  code can use #if defined(TABLE2) to test out new features.
#define TABLE2

/* getaline.c */
char *getaline(stream);              // deprecated - table_line0() does this
char *getsline(stream, string *);    // deprecated - was never used, getline() does this

/* gettab.c */
int get_atable(stream , int , int *, real **, int);
int get_itable(stream , int , int *, int **, int);
int get_ftable(stream , int , int *, string *, real **, int);

/* table.c */
int get_line(stream, string);		/* should be deprecated in V2 -> getline() */
void parse(int, string, double *, int);
void strinsert(string, string, int);
int iscomment(string);
void sanitize(string);

//  For a new table system we need a table struct, and each table
//  will have some columns, with properties that we hide in another
//  struct.

typedef struct {
  
  string *name;    // name of the column
  string *unit;    // units 
  int type;        // type (integer, real, string)

} column, *columnptr;
   
typedef struct {
  
  int  mode;        // I/O mode  (0=streaming, 1=all-in-memory, ...)
  int  type;        // type of table (SSV, TSV, CSV, ECSV, ipac, ....)
  size_t  nh;       // number of header/comment rows
  size_t  nr;       // number of rows
  size_t  nc;       // number of columns

  columnptr *cols;  // optional column designators

  string name;      // filename, if used
  stream str;       // stream, if used

  string *comments; // pointer to 'nh' (header/comment) lines
  string *lines;    // pointer to 'nr' lines (depends on mode)   lines[0], lines[1], ....

  size_t linelen;   // see Posix getline(3)
  char  *line;      // see Posix getline(3)
  
} table, *tableptr;

//  this API is not final yet; see also table.3 for a proposed API
table  *table_open(stream instr, int mode);
table  *table_open1(stream instr, int mode, int nlines);
table  *table_cat(table *tp1, table *tp2, int mode);
void    table_close(tableptr tptr);
string  table_line(tableptr tptr);
ssize_t table_line1(tableptr tptr, char **line, size_t *linelen, int newline);
size_t  table_nrows(tableptr tprt);
size_t  table_ncols(tableptr tprt);
string  table_row(tableptr tptr, int row);
string *table_rowsp(tableptr tptr, int row);
mdarray2 table_md2rc(table *t, int nrow, int *rows, int ncol, int *cols);  // a[row][col]
mdarray2 table_md2cr(table *t, int ncol, int *cols, int nrow, int *rows);  // a[col][row]





#endif



