/*
 * Various ansi prototypes for tables
 */
 
/* getaline.c */
char *getaline(stream);
char *getsline(stream, string *);

/* gettab.c */
int get_atable(stream , int , int *, real **, int);
int get_ftable(stream , int , int *, string *, real **, int);

/* table.c */
int get_line(stream, string);		/* should be deprecated */
void parse(int, string, double *, int);
void strinsert(string, string, int);
int iscomment(string);

