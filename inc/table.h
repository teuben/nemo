/*
 *  various support for table I/O
 *
 *  Additional support is given via burststring.c and extstring.c
 *
 */

#ifndef _h_table
#define _h_table
 
/* getaline.c */
char *getaline(stream);
char *getsline(stream, string *);

/* gettab.c */
int get_atable(stream , int , int *, real **, int);
int get_itable(stream , int , int *, int **, int);
int get_ftable(stream , int , int *, string *, real **, int);

/* table.c */
int get_line(stream, string);		/* should be deprecated */
void parse(int, string, double *, int);
void strinsert(string, string, int);
int iscomment(string);


/* table, *tableptr:
   a structure containing a string table - new 2020 style

*/
   
typedef struct {
  int  mode;        // I/O mode  (streaming, all-in-memory, ...)
  int  type;        // type of table (SSV, TSV, CSV, ....)
  int    nr;        // number of rows
  int    nc;        // number of columns

  string name;      // filename, if used
  stream str;       // stream, if used

  string *lines;    // pointer to 'nr' lines (depending on mode)
} table, *tableptr;

#endif
