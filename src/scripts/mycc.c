/***************************************************************/
/* File: mycc.c                                                */
/* Last modified: Fri May  8 1987 16:00:50, Josh Barnes, IAS.  */
/* ----------------------------------------------------------- */
/*     This program augments the C compiler and allows the     */
/* user to specify (1) a new search path to be used for        */
/* #include files marked as library inclusions (i.e., using    */
/* the #include <...> syntax), and (2) a directory containing  */
/* additional archives of library routines which are to be     */
/* searched at load time.  Moreover, each of these functions   */
/* is specified through the setting of an environment variable */
/* (INCLUDE and LIBRARY, respectively) so that the inclusion   */
/* of the new libraries is automatic, requiring no additional  */
/* user specification.                                         */
/*                                                             */
/*     When the compiler encounters an #include line, the      */
/* INCLUDE directory path, if specified, is searched before    */
/* the system directory /usr/include.  This makes it possible  */
/* to supersede header files in the standard library with      */
/* private copies.  If the environment variable INCLUDE is     */
/* not set, the standard search path is used.  However, the    */
/* syntax of the -I switch is extended so that search paths    */
/* and the ~ notation may also be included here.               */
/*                                                             */
/*     If the environment variable LIBRARY is set, cc looks    */
/* for library files (specified by -l<file>) in the LIBRARY    */
/* directory before checking the standard library directories  */
/* such as /usr/lib.					       */
/***************************************************************/

#include <stdio.h>
#include <ctype.h>
#include <pwd.h>
#include <strings.h>

/***************************************************************/
/* Include some definitions from defs.h                        */
/***************************************************************/

#define TRUE  1
#define FALSE 0
#define void int

typedef char bool;
typedef char *string;
typedef void (*proc)();

char *malloc();
char *getenv();

/***************************************************************/
/* Package definitions                                         */
/***************************************************************/

#define MAXARG 1000
#define CCOM "/bin/cc"

/***************************************************************/
/* Package variables                                           */
/***************************************************************/

static string newargv[MAXARG];
static int newargc;

/***************************************************************/
/* Local function declarations                                 */
/***************************************************************/

static void mappath(/* fn, path */);
static void includepath(/* dir */);
static void includelib(/* dir */);
static string expandtilde(/* name */);
static string sconc(/* s1, s2 */);
static string substr(/* s, p1, p2 */);



main(argc, argv)
int argc;
string argv[];
{
    register int i;
    string incpath, libpath;
    bool libflag;

    incpath = getenv("INCLUDE");
    libpath = getenv("LIBRARY");
    newargc = 0;
    newargv[newargc++] = argv[0];
    if (incpath != NULL) mappath(includepath, incpath);
    libflag = (libpath != NULL);
    for (i = 1; i < argc; i++) {
	if (strcmp(argv[i], "-c") == 0) libflag = FALSE;
	if (strcmp(argv[i], "-E") == 0) libflag = FALSE;
	if (strcmp(argv[i], "-S") == 0) libflag = FALSE;
	if (strncmp(argv[i], "-I", 2) == 0)
	    mappath(includepath, &argv[i][2]);
	else
	    newargv[newargc++] = argv[i];
    }
    if (libflag) mappath(includelib, libpath);
    newargv[newargc] = NULL;
    execv(CCOM, newargv);
    fprintf(stderr, "Can't execute cc\n");
    exit(1);
}



/***************************************************************/
/* mappath(fn, path)                                           */
/*                                                             */
/*     Maps a function over each element in the search path,   */
/* which has the form of files or directories separated by     */
/* colons.  Leading and trailing white space is removed        */
/* prior to the internal call to fn.                           */
/***************************************************************/

static void mappath(fn, path)
proc fn;
string path;
{
    register char *cp;
    string start, finish, localpath;

    localpath = sconc(path, ":");
    for (start = localpath; finish = index(start, ':'); start = finish+1) {
	while (isspace(*start)) start++;
        for (cp = finish-1; cp > start && isspace(*cp); cp--);
	*++cp = '\0';
	if (start != finish) (*fn)(expandtilde(start));
    }
}

/***************************************************************/
/* These are the functions actually mapped in this program     */
/*                                                             */
/*     includepath  -- adds an -Iname switch to the arguments  */
/*     includelib   -- adds an -Lname switch to the arguments  */
/***************************************************************/

static void includepath(dir)
string dir;
{
    newargv[newargc++] = sconc("-I", dir);
}

static void includelib(dir)
string dir;
{
    newargv[newargc++] = sconc("-L", dir);
}



/***************************************************************/
/* newname = expandtilde(oldname);                             */
/*                                                             */
/*     Given a filename, returns it after performing tilde     */
/* expansion.                                                  */
/***************************************************************/

static string expandtilde(name)
string name;
{
    string slashpos, homedir, newname;
    struct passwd *pw;

    if (*name != '~') return (name);
    slashpos = index(name, '/');
    if (slashpos == NULL) slashpos = name + strlen(name);
    if (slashpos - name == 1) {
	homedir = getenv("HOME");
	if (homedir == NULL) homedir = getpwuid(getuid())->pw_dir;
    } else {
	homedir = substr(name, 1, slashpos-name-1);
	pw = getpwnam(homedir);
        if (pw == NULL) {
	    fprintf(stderr, "Error [cc]: No such user -- %s\n", homedir);
	    exit(1);
	}
	homedir = pw->pw_dir;
    }
    newname = sconc(homedir, slashpos);
    return (newname);
}



/***************************************************************/
/* s = sconc(s1, s2);                                          */
/*                                                             */
/*     Concatenates two strings and returns the result in      */
/* dynamically-allocated storage.                              */
/***************************************************************/

static string sconc(s1, s2)
string s1, s2;
{
    int l;
    string result;

    result = (string) malloc((unsigned) (l = strlen(s1)) + strlen(s2) + 1);
    strcpy(result, s1);
    strcpy(result+l, s2);
    return (result);
}



/***************************************************************/
/* s = substr(s, p1, p2);                                      */
/*                                                             */
/*     Returns the substring of s extending from the integer   */
/* indices p1 and p2 (inclusive).  The following edge cases    */
/* apply:                                                      */
/*                                                             */
/*      if p1 < 0 then p1 <- 0;                                */
/*      if p2 > strlen(s) then p2 <- strlen(s);                */
/*      if p1 > p2 then return "";                             */
/***************************************************************/

static string substr(s, p1, p2)
string s;
int p1, p2;
{
    int l, i;
    string result;

    l = strlen(s);
    if (p1 < 0) p1 = 0;
    if (p2 >= l) p2 = l - 1;
    if ((l = p2 - p1 + 1) <= 0) return ("");
    result = (string) malloc ((unsigned) l + 1);
    for (i = 0; i < l; i++)
	result[i] = s[p1+i];
    result[l] = 0;
    return (result);
}
