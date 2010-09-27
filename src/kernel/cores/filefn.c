/***************************************************************/
/* File: filefn.c                                              */
/* Last modified on Sat Jan 11 09:26:56 1986 by roberts        */
/* Last modified on Sat Dec 06          1986 by josh           */
/* 		    Sun Oct 07 13:21:10 1990 by Peter          */
/*          ANSI: string.h is now used - index,rindex replaced */
/* Various difficulties in compiling with Turbo C              */
/* and made gcc2.0 silent by declaring more 		       */
/*                  20-nov-94 added () to shutup gcc/lint      */
/*                  12-apr-95 no more ARGS                     */
/*                  20-jun-01 gcc3 - removed old BORLAND code  */
/*                  17-mar-06 added fullname                   */
/*                  27-Sep-10 MINGW32/WINDOWS support (JCL)    */
/***************************************************************/

#include <stdinc.h>
#include <strlib.h>
#include <filefn.h>

#include <ctype.h>
#ifndef __MINGW32__
#include <pwd.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>

#include <unistd.h>     
#include <limits.h>

#ifndef MAXPATHLEN
#define MAXPATHLEN    PATH_MAX
#endif

#define DIR_SEP       '/'
/***************************************************************/
/* Local forward function declarations                         */
/***************************************************************/

local string expandtilde ( string );
local string checkexists ( string, string );




/***************************************************************/
/* str = root(filename);                                       */
/* str = extension(filename);                                  */
/* str = head(filename);                                       */
/* str = tail(filename);                                       */
/*                                                             */
/*     Each of these routines returns a component of a file    */
/* name and are equivalent to the C shell substitution         */
/* characters r, e, h, and t, respectively.                    */
/*							       */
/* root(), extension() fixed to behave consistently when given */
/* a file name with no extension at all.  JEB  06 Dec 1986     */
/***************************************************************/

string root(string filename)
{
    char *dotpos;

    dotpos = strrchr(filename, '.');
    if (dotpos != NULL && strchr(dotpos, DIR_SEP) != NULL)
	dotpos = NULL;
    if (dotpos == NULL)
	return (scopy(filename));
    else
        return (substr(filename, 0, (dotpos - filename) - 1));
}

string extension(string filename)
{
    char *dotpos;

    dotpos = strrchr(filename, '.');
    if (dotpos != NULL && strchr(dotpos, DIR_SEP) != NULL)
	dotpos = NULL;
    if (dotpos == NULL)
	return ("");
    else
        return (scopy(dotpos + 1));
}

string head(string filename)
{
    char *slashpos;

    slashpos = strrchr(filename, DIR_SEP);
    if (slashpos == NULL)
	return ("");
    else
        return (substr(filename, 0, (slashpos - filename) - 1));
}

string tail(string filename)
{
    char *slashpos;

    slashpos = strrchr(filename, DIR_SEP);
    if (slashpos == NULL)
	return (scopy(filename));
    else
        return (scopy(slashpos + 1));
}


/*
 * fullname - return the full name of a file
 * as to make it easier to open it when
 * a program uses chdir(2)
 *
 * TODO:    
 */

string fullname(string filename)
{
  char pathname[MAXPATHLEN];
  char pathsep[2];
  int n1,n2;
  char *outname;

  if (*filename == DIR_SEP)
    return scopy(filename);
  if (getcwd(pathname,MAXPATHLEN) == 0)
    error("Directory name too long (MAXPATHLEN=%d)",MAXPATHLEN);
  n1 = strlen(pathname) + 1;
  n2 = strlen(filename) + 1;
  sprintf(pathsep,"%c",DIR_SEP);
  if (n1+n2 < MAXPATHLEN) {
    strcat(pathname,pathsep);
    strcat(pathname,filename);
    return scopy(pathname);
  } else {
    outname = (char *) allocate(n1+n2+1);
    sprintf(outname,"%s%c%s",pathname,DIR_SEP,filename);
    return outname;
  }
  /* actually never reaches here */
  return 0;
}



/***************************************************************/
/* newname = defext(oldname, ".xxx")                           */
/*                                                             */
/*     The defext routine adds an extension to a file name     */
/* if none already exists.  Alternatively, if the extension    */
/* field begins with a *, any old extension in the first       */
/* filename is replaced with the given extension.              */
/*                                                             */
/*     defext(filename, ".xxx")   --  add .xxx if no ext       */
/*     defext(filename, "*.xxx")  --  force .xxx as ext        */
/*                                                             */
/* Note:  defext returns a pointer to dynamically-allocated    */
/* string storage which is never freed.  This is necessary     */
/* to ensure safety on multiple calls in a single statement.   */
/***************************************************************/

string defext(string filename, string ext)
{
    register char c, *cp;
    char *xp;
    bool forceext;

    if ((forceext = (ext[0] == '*'))) ext++;  /* LINT: questionable assignment */
    xp = NULL;
    for (cp = filename; (c = *cp); cp++) {    /* LINT: questionable assignment */
	switch (c) {
	    case DIR_SEP : 
            case ':'     : xp = NULL; break;
	    case '.'     : xp = cp;
	}
    }
    if (xp == NULL) {
	forceext = TRUE;
	xp = cp;
    }
    if (forceext)
	return (sconc(substr(filename, 0, xp-filename-1), ext));
    else
	return (scopy(filename));
}



/***************************************************************/
/* stream = pathopen(path, filename, mode);                    */
/*                                                             */
/*     The pathopen routine is used to open files using a      */
/* search path similar to that used, for example, by csh       */
/* in searching for a command.  The pathopen routine has       */
/* the same structure as fopen in the standard library and     */
/* the filename and mode arguments are the same as in that     */
/* call.  The path argument consists of a list of directories  */
/* which are prepended to the filename, unless the filename    */
/* begins with either a / or a ~.  The directories in the      */
/* list are separated by colons as in the definition of the    */
/* PATH environment variable.  White space and empty fields    */
/* are ignored to simplify formatting of paths in a definition */
/* file.                                                       */
/*                                                             */
/*     After each directory name has been added, the           */
/* pathopen function performs ~ expansion in the same form as  */
/* csh.  The path argument may be NULL, in which case no       */
/* directories are prepended.  This is useful if ~ expansion   */
/* is the only required function.                              */
/*                                                             */
/*     The pathopen function returns an open stream to         */
/* the indicated file, or NULL, if no existing file is         */
/* found.                                                      */
/***************************************************************/

stream pathopen(string path, string filename, string mode)
{
    return      /* PPAP */
        (stream) _mappath((strfn) fopen, path, filename, mode);
}



/***************************************************************/
/* newname = pathfind(path, filename);                         */
/*                                                             */
/*     The pathfind routine is similar to pathopen, except     */
/* that it does not try to open the file but instead returns   */
/* the full name of the first file that exists along the path, */
/* or NULL if none exist.                                      */
/***************************************************************/

string pathfind(string path, string filename)
{
    return
        (string) _mappath((strfn) checkexists, path, filename, (string) NULL);
}

/***************************************************************/
/* Internal routine to determine if a file exists              */
/***************************************************************/

local string checkexists(string name, string dummy)
{
    permanent struct stat statbuf;

    return ((stat(name, &statbuf) == 0) ? name : NULL);
}



/***************************************************************/
/* s = _mappath(fn, path, filename, arg);                      */
/*                                                             */
/*     Maps the string function fn over each path/filename     */
/* combination, stopping when fn(name, arg) returns non-null.  */
/* That value is returned.  If the function always returns a   */
/* NULL value for each element in the path, NULL is returned   */
/* from _mappath.                                              */
/***************************************************************/

string _mappath(strfn fn, string path, string filename, string arg)
{
    register char *cp;
    char c;
    string start, finish, localpath, tempname, result;

    if (filename == NULL) filename = "";
    c = filename[0];
    if (path == NULL || path[0] == '\0' || c == '~' || c == DIR_SEP)
	return (*fn)(expandtilde(filename), arg);
    localpath = sconc(path, ":");
    result = NULL;
    start = localpath;
    while (result == NULL && (finish = strchr(start, ':')) != NULL) {
	while (isspace(*start)) start++;
        for (cp = finish-1; cp > start && isspace(*cp); cp--); /*LINT: empty for*/
	*++cp = '\0';
	if (start != finish || strlen(localpath) == 1) {
	    tempname = sconc(start, sconc("/", filename));
	    result = (*fn)(expandtilde(tempname), arg);
	}
	start = finish+1;
    }
    return result;
}



/***************************************************************/
/* newname = expandtilde(oldname);                             */
/*                                                             */
/*     Given a filename, returns it after performing tilde     */
/* expansion.                                                  */
/* For MSDOS this always returns the oldname                   */
/***************************************************************/

local string expandtilde(string name)
{
#if defined(__BORLANDC__ ) || defined(__MINGW32__) 
    if (*name != '~') return (name);
    warning("Cannot parse a tilde (~) in a filename, returning %s",name);
    return(name);
#else
    string slashpos, homedir, newname;
    struct passwd *pw;

    if (*name != '~') return (name);

    slashpos = strchr(name, DIR_SEP);
    if (slashpos == NULL) slashpos = name + strlen(name);
    if (slashpos - name == 1) {
	homedir = getenv("HOME");
	if (homedir == NULL) homedir = getpwuid(getuid())->pw_dir;  /*3b1*/
    } else {
	homedir = substr(name, 1, slashpos-name-1);
	pw = getpwnam(homedir);					    /*3b1*/
        if (pw == NULL)
	    error("expandtilde: no such user: %s\n", homedir);
	homedir = pw->pw_dir;
    }
    newname = sconc(homedir, slashpos);
    return (newname);
#endif
}

#ifdef TESTBED
#include <getparam.h>
string defv[] = {
    "path=\n		full path to test",
    "name=\n		name",
    "fullname=\n        find fullname",
    "VERSION=1.2\n      17-mar-06 PJT",
    NULL,
};

string path, name;

#ifndef MAXLINE
#define MAXLINE 256
#endif

char line[MAXLINE];

nemo_main()
{
    stream infile;

    if (hasvalue("fullname")) {
      printf("%s\n",fullname(getparam("fullname")));
      stop(0);
    }

    path = getparam("path");
    name = getparam("name");
    printf("head(\"%s\") = \"%s\"\n", path, head(path));
    printf("tail(\"%s\") = \"%s\"\n", path, tail(path));
    printf("root(\"%s\") = \"%s\"\n", name, root(name));
    printf("extension(\"%s\") = \"%s\"\n", name, extension(name));
    printf("defext(\"%s\", \".xxx\") = \"%s\"\n", name, defext(name, ".xxx"));
    printf("defext(\"%s\", \"*.xx\") = \"%s\"\n", name, defext(name, "*.xx"));
    infile = pathopen(path, name, "r");
    if (infile == NULL)
	error("Can't pathopen file %s\n", name);
    printf("----------------\n");
    while (fgets(line, MAXLINE, infile) != NULL)
	fputs(line, stdout);
}

#endif
