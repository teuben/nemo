/*
 *  various support for table I/O
 *
 *  feb-2022     Table I/O - V2.0
 *
 *  Additional support is given via burststring.c and extstring.c
 *
 */

#ifndef _h_table
#define _h_table

//  code can use #if defined(TABLE2) to test out new features.
#define TABLE2
 
/* getaline.c */
char *getaline(stream);
char *getsline(stream, string *);

/* gettab.c */
int get_atable(stream , int , int *, real **, int);
int get_itable(stream , int , int *, int **, int);
int get_ftable(stream , int , int *, string *, real **, int);

/* table.c */
int get_line(stream, string);		/* should be deprecated in V2 -> getline() */
void parse(int, string, double *, int);
void strinsert(string, string, int);
int iscomment(string);

//  For a new table system we need a table struct, and each table
//  will have some columns, with properties that we hide in another
//  struct.

typedef struct {
  string *name;    // name of the column
  string *unit;    // units 
  int type;        // type (integer, real)

} column, *columnptr;
   
typedef struct {
  int  mode;        // I/O mode  (streaming, all-in-memory, ...)
  int  type;        // type of table (SSV, TSV, CSV, ....)
  int    nr;        // number of rows
  int    nc;        // number of columns

  columnptr *cols;  // optional column designators

  string name;      // filename, if used
  stream str;       // stream, if used

  string *lines;    // pointer to 'nr' lines (depending on mode)
} table, *tableptr;

//  this API is not final yet
table  *table_open(stream instr, int mode);
void    table_close(tableptr tptr);
ssize_t table_line(tableptr tptr, char **line, size_t *linelen);

#endif



