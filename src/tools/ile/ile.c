/*
                          COPYRIGHT 1988
              Evans & Sutherland Computer Corporation
                        Salt Lake City, Utah
                        All Rights Reserved.

     THE INFORMATION  IN  THIS  SOFTWARE  IS  SUBJECT  TO  CHANGE
     WITHOUT  NOTICE  AND SHOULD NOT BE CONSTRUED AS A COMMITMENT
     BY  EVANS  &  SUTHERLAND.   EVANS  &  SUTHERLAND   MAKES  NO
     REPRESENTATIONS  ABOUT  THE SUITABILITY OF THIS SOFTWARE FOR
     ANY PURPOSE.  IT IS SUPPLIED  "AS  IS"  WITHOUT  EXPRESS  OR
     IMPLIED WARRANTY.

     IF THE SOFTWARE IS MODIFIED IN A MANNER CREATING  DERIVATIVE
     COPYRIGHT  RIGHTS,  APPROPRIATE LEGENDS MAY BE PLACED ON THE
     DERIVATIVE WORK IN ADDITION TO THAT SET FORTH ABOVE.

     Permission  to  use,  copy,  modify,  and  distribute   this
     software  and  its documentation for any purpose and without
     fee is hereby granted, provided  that  the  above  copyright
     notice  appear  in  all  copies  and that both the copyright
     notice and this permission notice appear in supporting docu-
     mentation,  and  that  the name of Evans & Sutherland not be
     used in advertising or publicity pertaining to  distribution
     of the software without specific, written prior permission.

Written by:

                        Robert C. Pendleton <bobp@hal.com>

Grateful acknowledgement is made of code and ideas contributed by

 	Ian Donaldson,
 	Department of Communications & Electronic Engineering,
 	Royal Melbourne Institute of Technology,
 	Melbourne, Australia.

$Header$
*/

/*
ile is compiled using:

cc ile.c -o ile -ltermcap
*/

#include <stdio.h>
#include <fcntl.h>
#include <sgtty.h>
#include <signal.h>
#include <string.h>
#include <strings.h>
#include <pwd.h>
#include <utmp.h>
#include <errno.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <sys/dir.h>
#include <sys/file.h>
#include <sys/time.h>

/*------------------------------------------------------------------*/
/*
Definitions of system stuff.
*/

extern int errno;

long lseek();
char *malloc();
char *realloc();
time_t time();

/*------------------------------------------------------------------*/

#define FALSE 0
#define TRUE  1

#define READ  0
#define WRITE 1
#define ERROR 2

#define BUFFER_SIZE 255

#define HISTORY_SIZE 101

#define USER_NAME_SIZE 8

#define EOL (-2)
/*------------------------------------------------------------------*/

/* special characters used by ile */

#define del '\177'

#define CA '\1'
#define CB '\2'
#define CC '\3'
#define CD '\4'
#define CE '\5'
#define CF '\6'
#define bel '\7'
#define bs '\10'
#define CI '\11'
#define nl '\12'
#define CK '\13'
#define CL '\14'
#define cr '\15'
#define CN '\16'
#define CO '\17'
#define CP '\20'
#define CQ '\21'
#define CR '\22'
#define CS '\23'
#define CT '\24'
#define CU '\25'
#define CV '\26'
#define CW '\27'
#define CX '\30'
#define CY '\31'
#define CZ '\32'

#define esc '\33'

/*------------------------------------------------------------------*/
/* areas and varaibles used to get termcap information */

char *getenv();

char *tgetnum();
char *tgetflag();
char *tgetstr();

char termcap_entry[1024];	/* termcap entry for the users terminal */
char term_seqs[1024];		/* area to store control sequences in */
char *where = term_seqs;

char *cle;			/* move cursor left one space */

char *cce;			/* clear to end of line */

char *cbl;			/* audible bell */

char *cnl;			/* new line character */

char *ccr;			/* carriage return */

/*------------------------------------------------------------------*/
/*
The value of HOME
*/

char *homedir = NULL;

/*------------------------------------------------------------------*/
/*
The current working directory as set by query_path.
initialized to PWD
*/

char currentdir[MAXNAMLEN + 1] = "";

/*------------------------------------------------------------------*/

/* tty status flags */
/*
  The original tty status flags are stored so that they can be
  restored when ile exits.
*/

struct sgttyb tty_sgttyb;
struct tchars tty_tchars;
struct ltchars tty_ltchars;
struct winsize tty_winsize;
int windowchanged;
int tty_ldisc;
int tty_mode;

/*------------------------------------------------------------------*/

/* file descriptors for tty and pty */

int master_pty;
int slave_tty;

/* the names of the tty and pty opened by getpty */

char ttydev[] = "/dev/ttyxx";
char ptydev[] = "/dev/ptyxx";

/* path and name of the lock file */

char lock[] = "/tmp/ile.lock";

/*------------------------------------------------------------------*/
/*
  getpty opens a pty, storing file descriptors in pty and tty.
  It trys pairs in order until it finds a pair that is not in use.
*/

getpty(pty, tty)
    int *pty;
    int *tty;

{
    int devindex;
    int letter;

    static char ptychar1[] = "pqrstuvwxyz";
    static char ptychar2[] = "0123456789abcdef";

    letter = 0;
    while (letter < 11)
    {
	ttydev[strlen(ttydev) - 2] = ptychar1[letter];
	ptydev[strlen(ptydev) - 2] = ptychar1[letter];
	letter++;

	devindex = 0;
	while (devindex < 16)
	{
	    ttydev[strlen(ttydev) - 1] = ptychar2[devindex];
	    ptydev[strlen(ptydev) - 1] = ptychar2[devindex];
	    devindex++;

	    if ((*pty = open(ptydev, O_RDWR)) >= 0)
	    {
		if ((*tty = open(ttydev, O_RDWR)) >= 0)
		{
		    return;
		}
		else
		{
		    (void) close(*pty);
		}
	    }
	}
    }

    (void) fprintf(stderr, "ile: unable to allocate pty/tty pair\n");
    exit(1);
    /* NOTREACHED */
}

/*------------------------------------------------------------------*/
/*
Termcap entries may have a sequences of digits optionally followed
by a '*' in front of the actual sequence. This routine increments
the pointer past this information.
*/
void
strip(ptr)
    char **ptr;
{
    while (('0' <= **ptr) && (**ptr <= '9'))
    {
	(*ptr)++;
    }

    if (**ptr == '*')
    {
	(*ptr)++;
    }
}
/*------------------------------------------------------------------*/
/*
Set up everything needed to use the control sequences from the
termcap entry for the terminal.
*/
void
get_termcap()
{
    char *terminal_type;	/* type of terminal */

    /* get the terminal name */

    terminal_type = getenv("TERM");

    /* get termcap entry */

    if (tgetent(termcap_entry, terminal_type) < 1)
    {
	(void) fprintf(stderr, "ile: can't find %s\n", terminal_type);
	exit(1);
	/* NOTREACHED */
    }

    /* get the control sequences ile needs */

    if ((cbl = tgetstr("bl", &where)) == NULL)
    {
	cbl = "\7";
    }

    if ((cnl = tgetstr("nl", &where)) == NULL)
    {
	cnl = "\n";
    }

    if ((ccr = tgetstr("cr", &where)) == NULL)
    {
	ccr = "\r";
    }

    if ((cle = tgetstr("le", &where)) == NULL)
    {
	if (tgetflag("bs"))
	{
	    cle = "\b";
	}
    }

    if ((cle == NULL) ||
	((cce = tgetstr("ce", &where)) == NULL))
    {
	(void) fprintf(stderr,
	    "ile: can't run on %s (need capabilities \"le\" and \"ce\")\n",
	    terminal_type);
	exit(1);
	/* NOTREACHED */
    }

    /* strip timing info from strings */

    strip(&cle);
    strip(&cce);
    strip(&cbl);
    strip(&cnl);
    strip(&ccr);
}
/*------------------------------------------------------------------*/
/*
If the window changes size, tell the slave_tty about it.
*/
void
change_window()
{
    int pgrp;

    (void) ioctl(READ, TIOCGWINSZ, &tty_winsize);
    (void) ioctl(slave_tty, TIOCSWINSZ, &tty_winsize);

    (void) ioctl(slave_tty, TIOCGPGRP, (char *) &pgrp);
    (void) killpg(pgrp, SIGWINCH);

    /* note the change so that we don't die after select */

    windowchanged = TRUE;
}
/*------------------------------------------------------------------*/
/*
 * set/clear the utmp slot for the pty
 */
int
setutmp(fd, set)
{
    int old0;
    int old1;
    int old2;
    int slot;
    int f;
    struct utmp ut;

    /* Must make fd's 0,1,2 correspond to slave_tty for ttyslot() to
     * function.  Ugh!  Why doesn't ttyslot() accept a fd argument?
     * 
     * save fd's */

    old0 = dup(0);
    old1 = dup(1);
    old2 = dup(2);

    if (old0 == -1 || old1 == -1 || old2 == -1)
    {
	perror("ile: dup");
	return (-1);		/* file table full ? */
    }

    /* set fd's 0,1,2 for ttyslot() */

    (void) dup2(fd, 0);
    (void) dup2(fd, 1);
    (void) dup2(fd, 2);

    slot = ttyslot();

    /* put the fd's back */

    (void) dup2(old0, 0);
    (void) dup2(old1, 1);
    (void) dup2(old2, 2);

    (void) close(old0);
    (void) close(old1);
    (void) close(old2);

    if (slot < 0)
    {
	(void) fprintf(stderr, "ile: don't know where you are\n");
	return (-1);
    }

    f = open("/etc/utmp", O_WRONLY);
    if (f == -1)
    {
	return (-1);
    }

    bzero((char *) &ut, sizeof(ut));

    if (set)
    {
	struct passwd *pw;
	char *cp;

	pw = getpwuid(getuid());
	if (pw == 0)
	{
	    (void) fprintf(stderr, "ile: who are you?\n");
	    (void) close(f);
	    return (-1);
	}

	/* skip "/dev/" */

	cp = rindex(ttydev, '/');
	if (cp == 0)
	{
	    cp = ttydev;
	}
	else
	{
	    cp++;
	}

	(void) strncpy(ut.ut_line, cp, sizeof(ut.ut_line));
	(void) strncpy(ut.ut_name, pw->pw_name, sizeof(ut.ut_line));
	(void) time(&ut.ut_time);
    }
    (void) lseek(f, (long) (slot * sizeof(struct utmp)), L_SET);
    (void) write(f, (char *) &ut, sizeof(ut));
    (void) close(f);

    return (0);
}
/*------------------------------------------------------------------*/
/*
clean up and leave.

This function is bound to the SIGCHLD signal so that when the
child process exits, so does ile. It is also called when an exception
is detected by select() in ile().
*/
void
clean_up()
{
    int pgrp;

    /* kill off the child process */

    (void) ioctl(slave_tty, TIOCGPGRP, (char *) &pgrp);
    (void) killpg(pgrp, SIGTERM);

    /* restore terminal status */

    (void) ioctl(READ, TIOCSETP, &tty_sgttyb);
    (void) ioctl(READ, TIOCSETC, &tty_tchars);
    (void) ioctl(READ, TIOCSLTC, &tty_ltchars);
    (void) ioctl(READ, TIOCLSET, &tty_mode);

    /* "logout" the user */

    (void) setutmp(slave_tty, FALSE);

    /* clean up the tty/pty pair */

    (void) close(master_pty);
    (void) close(slave_tty);

    /* make things look nice */

    fputs(cnl, stdout);

    exit(0);
    /* NOTREACHED */
}
/*------------------------------------------------------------------*/
/*
Write a line to the slave_tty.
Get the slave_tty parameters, turn off echo, send the line to the
slave_tty, restore the slave_tty paramters to the way they were
before. If echo was already off, this will have no effect.
*/
write_line(line, length)
    char *line;
    int length;
{
    struct sgttyb params;
    struct sgttyb new_params;

    /* get the current parameters */

    (void) ioctl(slave_tty, TIOCGETP, &params);
    new_params = params;
    new_params.sg_flags &= ~ECHO;

    /* turn off echo so we don't see the characters twice */

    (void) ioctl(slave_tty, TIOCSETP, &new_params);

    (void) write(master_pty, line, length);

    /* set the parameters back the way they were */

    (void) ioctl(slave_tty, TIOCSETN, &params);
}

/*------------------------------------------------------------------*/
/*
The editing routines are called through the edit variable. This allows
the quote and escape commands to be implemented as a straight forward
state machine instead of requiring state flags and complex switch
statements.
*/
/*------------------------------------------------------------------*/

/* line edit buffer */

static char line[BUFFER_SIZE];

static int point;		/* insertion point */
static int length;		/* total chars in buffer */

/* procedure to edit next character */

void (*edit) ();

/* history buffer */

struct
{
    int length;
    char *line;
} hist[HISTORY_SIZE];

int head;			/* insertion point */
int here;			/* current displayed line */

/*------------------------------------------------------------------*/
/*
The delimiter vector is used by the forward, backward, and delete
word operations to decide that a character is a delimiter.
*/
/*------------------------------------------------------------------*/

#define CHAR_SET_SIZE 127
#define CHAR_MASK 0177

char delimit[CHAR_SET_SIZE];

/*------------------------------------------------------------------*/
/*
The action_table is used to bind sequences of keys to operations or strings.
*/
/*------------------------------------------------------------------*/

typedef enum
{
    is_action, is_string
} action_type;

struct
{
    action_type flag;
    union
    {
	void (*action) ();
	char *string;
    } aors;
} action_table[4][CHAR_SET_SIZE];

/*------------------------------------------------------------------*/

void echo();
void echoline();
void cleartoend();
void clearline();
void backspace();
void quote_edit();
void edit_0();
void edit_1();
void edit_2();
void edit_3();
void bell();
void insert();

/*------------------------------------------------------------------*/
/*
The following routines are action routines that are executed by the
editor to carry out commands. Each routine has a single character
argument. Each routine is invoked with the character that caused it
to be invoked as its argument.

The argument isn't always useful, but it is included to provide a
consistent interface for the routines.
*/
/*------------------------------------------------------------------*/
/*
Given a specific directory and the starting string of a file name,
find the longest partial file name that starts with the substring.
*/
void
complete_file_name(dir, name)
    char *dir;
    char *name;
{
    DIR *dirp;
    struct direct *dp;

    int len;
    int maxlen;
    int oldlen;

    char oldname[MAXNAMLEN + 1];
    char newname[MAXNAMLEN + 1];

    if ((dir != NULL) &&
	(name != NULL) &&
	((oldlen = strlen(name)) > 0) &&
	((dirp = opendir(dir)) != NULL))
    {
	maxlen = oldlen;
	(void) strcpy(oldname, name);
	(void) strcpy(newname, name);

	/* find the longest name starting with name */

	for (dp = readdir(dirp); dp != NULL; dp = readdir(dirp))
	{
	    if (dp->d_name != NULL)
	    {
		len = strlen(dp->d_name);
		if ((maxlen < len) &&
		    (strncmp(oldname, dp->d_name, oldlen) == 0))
		{
		    maxlen = len;
		    (void) strcpy(newname, dp->d_name);
		}
	    }
	}

	rewinddir(dirp);

	/* find the longest common sub string */

	for (dp = readdir(dirp); dp != NULL; dp = readdir(dirp))
	{
	    if (dp->d_name != NULL)
	    {
		len = strlen(dp->d_name);
		if ((len <= maxlen) &&
		    (strncmp(oldname, dp->d_name, oldlen) == 0))
		{
		    for (;
			(oldlen < len) &&
			(strncmp(newname, dp->d_name, len) != 0);
			len--);

		    maxlen = len;
		    newname[maxlen] = '\0';
		}
	    }
	}

	if (strlen(name) != strlen(newname))
	{
	    /* return the extended name */

	    (void) strcpy(name, newname);
	}
	else
	{
	    /* no difference so beep */

	    bell('\0');
	}

	(void) closedir(dirp);
    }
}
/*------------------------------------------------------------------*/
/*
Hidden parameters to dirselect. They must be hidden because dirselect
is passed as an argument to scandir.
*/

static char *namep;
static int namelen;

/*------------------------------------------------------------------*/
/*
Passed to scandir. It is used to decide which files to display.
*/

static
dirselect(dp)
    struct direct *dp;
{
    return (strncmp(dp->d_name, namep, namelen) == 0);
}

/*------------------------------------------------------------------*/
/*
List all the files in a given directory that start with the string
passed as name.
*/
void
list_file_names(dir, name)
    char *dir;
    char *name;
{
    struct direct **dlist;
    int i;
    int nfiles;
    int colwidth;
    int cols;
    int ncols;
    int nlines;
    int alphasort();

    if (dir == NULL || name == NULL)
    {
	return;
    }

    cols = tty_winsize.ws_col;
    if (cols <= 0)
    {
	cols = 80;
    }

    namelen = strlen(name);
    namep = name;

    nfiles = scandir(dir, &dlist, dirselect, alphasort);

    /* determine the longest file name length */

    colwidth = 8;		/* minimum width */
    for (i = 0; i < nfiles; i++)
    {
	struct direct *dp;

	dp = dlist[i];

	if (dp->d_namlen > colwidth)
	{
	    colwidth = dp->d_namlen;
	}
    }

    colwidth++;			/* at least 1 space between them */

    /* print the names, sorted vertically per column */
    ncols = cols / colwidth;
    if (ncols == 0)
    {
	/* longest filename is wider than the screen */
	ncols = 1;
    }

    nlines = (nfiles + ncols - 1) / ncols;

    if (nfiles > 0)
    {
	for (i = 0; i < nlines; i++)
	{
	    int j;
	    int l;

	    l = 0;
	    for (j = 0; j < ncols; j++)
	    {
		int m;
		struct direct *dp;

		m = l + i;
		if (m >= nfiles)
		{
		    break;
		}

		dp = dlist[m];
		fputs(dp->d_name, stdout);

		if (j < (ncols - 1))
		{
		    int k;

		    for (k = dp->d_namlen; k < colwidth; k++)
		    {
			(void) fputc(' ', stdout);
		    }
		    l += nlines;
		}
	    }
	    fputs(ccr, stdout);
	    fputs(cnl, stdout);
	}
    }
    free((char *) dlist);
}
/*------------------------------------------------------------------*/
/*
Assuming that there is a file name under the cursor, return a path
and a file name. If there is no path name return "."
*/
int
get_dir_and_name(dir, name, username, userend, start, middle, tail)
    char *dir;
    char *name;
    char *username;
    int *userend;
    int *start;
    int *middle;
    int *tail;
{
    int dirlen;

    int newstart;

    int punlen;
    char pun[USER_NAME_SIZE + 1];
    struct passwd *userpwd;

    int i;

    /* set the default path and file name */

    dir[0] = '\0';
    name[0] = '\0';
    username[0] = '\0';

    /* search for the start of the file name */
    /* start will be left pointing to the first character of the path  */

    for ((*start) = point;
	((0 < (*start)) && (line[(*start) - 1] != ' '));
	(*start)--
	);

    /* search for the end of the file name */
    /* tail will be left pointing at the last character of the path */

    for ((*tail) = point - 1;
	(((*tail) < (length - 1)) && (line[(*tail) + 1] != ' '));
	(*tail)++
	);

    /* search for the middle of the file name */
    /* middle will be left pointing at the first character of the last
     * element of the path */

    for ((*middle) = (*tail) + 1;
	((0 < (*middle)) &&
	    (line[(*middle) - 1] != '/') &&
	    (line[(*middle) - 1] != ' '));
	(*middle)--
	);

    /* copy path from line to dir */

    /* what base path */

    newstart = (*start);

    if ((line[newstart] == '~') &&
	((newstart + 1) < length) &&
	(line[newstart + 1] == '/'))
    {
	/* "~/" means use the value of HOME */

	newstart++;
	(void) strcpy(dir, homedir);
    }
    else if (line[newstart] == '~')
    {
	/* "~username" means use the users login directory */

	/* search for the end of the user name */

	for ((*userend) = newstart,
	    punlen = 0;
	    (((*userend) < (length - 1)) &&
		(line[(*userend) + 1] != ' ') &&
		(line[(*userend) + 1] != '/'));
	    (*userend)++,
	    punlen++);

	/* make middle point to middle */

	if ((*start) == (*middle))
	{
	    (*middle) = (*start) + punlen + 1;
	}

	/* extract partial user name from line */

	(void) strncpy(pun, &line[newstart + 1], punlen);
	pun[punlen] = '\0';

	/* search passwd file for partial match */

	for (userpwd = getpwent();
	    userpwd != NULL;
	    userpwd = getpwent())
	{
	    if ((punlen <= strlen(userpwd->pw_name)) &&
		(strncmp(pun, userpwd->pw_name, punlen) == 0))
	    {

		/* we have a partial match, record it */

		if (strlen(dir) == 0)
		{
		    newstart = (*userend) + 1;
		    (void) strcpy(dir, userpwd->pw_dir);
		    (void) strcpy(username, userpwd->pw_name);
		}
		else
		{
		    /* second partial match, forget the first one. */

		    newstart = (*start);
		    dir[0] = '\0';
		    username[0] = '\0';
		    return (FALSE);
		}

	    }
	}
	(void) setpwent();
    }
    else if ((line[newstart] == '.') &&
	    ((newstart + 1) < length) &&
	(line[newstart + 1] == '/'))
    {
	/* if it's "./" use current dir */

	newstart++;
	(void) strcpy(dir, currentdir);
    }
    else if ((line[newstart] == '.') &&
	    ((newstart + 1) < length) &&
	    (line[newstart + 1] == '.') &&
	    ((newstart + 2) < length) &&
	(line[newstart + 2] == '/'))
    {
	/* if it's "../" strip off one name from currentdir and use that */

	newstart += 2;
	(void) strcpy(dir, currentdir);
	for (i = strlen(dir); (i > 0) && (dir[i] != '/'); i--)
	{
	    /* nothing */
	}
	dir[i] = '\0';
    }
    else if (line[newstart] != '/')
    {
	/* doesn't start with a "/"? use currentdir */

	(void) strcpy(dir, currentdir);
	(void) strcat(dir, "/");
    }

    /* add on the rest of the path */

    dirlen = strlen(dir);
    for (i = 0; i < ((*middle) - newstart); i++)
    {
	dir[dirlen + i] = line[newstart + i];
    }
    dir[dirlen + i] = '\0';

    /* copy file name from line to name */

    for (i = 0; i < ((*tail) - (*middle) + 1); i++)
    {
	name[i] = line[(*middle) + i];
    }
    name[i] = '\0';

    return (TRUE);
}
/*------------------------------------------------------------------*/
/*
Perform file name completion. Put the full path and file name in the
line.
*/
/*ARGSUSED*/
void
complete_file_full(ch)
    char ch;
{
    char dir[10 * (MAXNAMLEN + 1)];
    char name[MAXNAMLEN + 1];
    char username[USER_NAME_SIZE + 1];

    char newline[BUFFER_SIZE];
    int newlength;
    int newpoint;

    int userend;
    int start;
    int middle;
    int tail;

    int i;

    /* get the path and file name in the line */

    if (get_dir_and_name(dir,
	    name,
	    username,
	    &userend,
	    &start,
	    &middle,
	    &tail))
    {

	/* complete the file name if possible */

	complete_file_name(dir, name);

	/* create a new line */

	/* start with the line prefix */

	(void) strncpy(newline, line, start);
	newline[start] = '\0';

	/* add in the new path */

	(void) strcat(newline, dir);

	/* stick in the new file name */

	(void) strcat(newline, name);
	newpoint = strlen(newline);

	/* finish with the line postfix */

	(void) strncat(newline, &line[tail + 1], (length - tail - 1));
	newlength = strlen(newline);

	/* display the new line */

	clearline('\0');

	point = newpoint;
	length = newlength;
	(void) strncpy(line, newline, newlength);

	echoline(line, length);

	for (i = point; i < length; i++)
	{
	    backspace(line[i]);
	}
    }
    else
    {
	bell('\0');
    }
}
/*------------------------------------------------------------------*/
/*
Perform file name completion in much the same style as csh.
*/
/*ARGSUSED*/
void
complete_file(ch)
    char ch;
{
    char dir[10 * (MAXNAMLEN + 1)];
    char name[MAXNAMLEN + 1];
    char username[USER_NAME_SIZE + 1];

    int userend;
    int start;
    int middle;
    int tail;

    char newline[BUFFER_SIZE];
    int newlength;
    int newpoint;

    int userlen;
    int len;

    int i;

    /* get the path and file name in the line */

    if (get_dir_and_name(dir,
	    name,
	    username,
	    &userend,
	    &start,
	    &middle,
	    &tail))
    {
	/* how long is the user name */

	userlen = strlen(username);

	/* complete the file name if possible */

	complete_file_name(dir, name);
	/* create a new line */

	/* start with the line prefix */

	(void) strncpy(newline, line, start);
	newline[start] = '\0';

	/* add in the new username */

	if (userlen != 0)
	{
	    /* put in new user name */

	    (void) strcat(newline, "~");
	    (void) strcat(newline, username);
	    len = strlen(newline);

	    /* put in the existing path */

	    (void) strncat(newline, &line[userend + 1], middle - userend - 1);
	    newline[len + (middle - userend - 1)] = '\0';
	}
	else
	{
	    /* put in the existing path */

	    len = strlen(newline);
	    (void) strncat(newline, &line[start], middle - start);
	    newline[len + (middle - start)] = '\0';
	}

	/* stick in the new file name */

	(void) strcat(newline, name);
	newpoint = strlen(newline);

	/* finish with the line postfix */

	(void) strncat(newline, &line[tail + 1], (length - tail - 1));
	newlength = strlen(newline);

	/* display the new line */

	clearline('\0');

	point = newpoint;
	length = newlength;
	(void) strncpy(line, newline, newlength);

	echoline(line, length);

	for (i = point; i < length; i++)
	{
	    backspace(line[i]);
	}
    }
    else
    {
	bell('\0');
    }
}
/*------------------------------------------------------------------*/
/*
List the names of files that start with the directory path and
file name under the cursor.
*/
/*ARGSUSED*/
void
show_files(ch)
    char ch;
{
    static char divider[] = "----------";

    void retype_line();

    char dir[10 * (MAXNAMLEN + 1)];
    char name[MAXNAMLEN + 1];
    char username[USER_NAME_SIZE + 1];

    int userend;
    int start;
    int middle;
    int tail;

    if (get_dir_and_name(dir,
	    name,
	    username,
	    &userend,
	    &start,
	    &middle,
	    &tail))
    {

	fputs(ccr, stdout);
	fputs(cnl, stdout);

	fputs(divider, stdout);

	fputs(ccr, stdout);
	fputs(cnl, stdout);

	list_file_names(dir, name);

	fputs(divider, stdout);

	fputs(ccr, stdout);
	fputs(cnl, stdout);

	retype_line('\0');
    }
    else
    {
	bell('\0');
    }
}
/*------------------------------------------------------------------*/
/*
Make the gross assumption that the program we are talking to is
a shell and send "pwd\n" to it. Whatever comes back is saved as the
value of currentdir.
*/
/*ARGSUSED*/
void
query_path(ch)
    char ch;
{
    static char command[] = "pwd\n";
    char buffer[BUFFER_SIZE];
    int readfd;
    struct timeval timeout;
    int status;
    int cc;
    int i;

    /* send the command to the shell, we hope. */

    write_line(command, strlen(command));

    /* read the directory path back */

    readfd = 1 << master_pty;
    timeout.tv_sec = 2;
    timeout.tv_usec = 0;

    do
    {
	status = select(32,
	    (fd_set *) & readfd,
	    (fd_set *) NULL,
	    (fd_set *) NULL,
	    &timeout);

	if (0 < status)
	{
	    cc = read(master_pty, currentdir, sizeof(currentdir));

	    /* strip off trailing control chars and blanks */

	    for (i = cc - 1; (currentdir[i] <= ' ') && (i > 0); i--)
	    {
		currentdir[i] = '\0';
	    }

	    /* read the prompt so it can be ignored */

	    readfd = 1 << master_pty;
	    timeout.tv_sec = 2;
	    timeout.tv_usec = 0;

	    do
	    {
		status = select(32,
		    (fd_set *) & readfd,
		    (fd_set *) NULL,
		    (fd_set *) NULL,
		    &timeout);

		if (0 < status)
		{
		    cc = read(master_pty, buffer, sizeof(buffer));
		}
		else if ((-1 == status) && windowchanged)
		{
		    windowchanged = FALSE;
		}
	    } while (status == -1);
	}
	else if ((-1 == status) && windowchanged)
	{
	    windowchanged = FALSE;
	}
    } while (status == -1);
}
/*------------------------------------------------------------------*/
/*
Ring the bell on the terminal.
*/
/*ARGSUSED*/
void
bell(ch)
    char ch;
{
    fputs(cbl, stdout);
}
/*------------------------------------------------------------------*/
/*
Pass characters to the slave. Don't mess with them at all.
*/
void
pass(ch)
    char ch;
{
    (void) write(master_pty, &ch, 1);
}
/*------------------------------------------------------------------*/
/*
Insert a character at point in the line buffer. While we are at it
update the display to show the insertion.
*/
void
insert(ch)
    char ch;
{
    int i;

    if (length < (BUFFER_SIZE - 2))
    {

	/* display the character */

	echo(ch);

	/* redisplay the rest of the line */

	echoline(&line[point], (length - point));

	/* move the characters in the line buffer */
	/* and put the cursor back at point */

	for (i = length; i > point; i--)
	{
	    line[i] = line[i - 1];
	    backspace(line[i]);
	}

	/* add the character to the line buffer */
	/* and increment point and length */

	line[point] = ch;
	length++;
	point++;
    }
    else
    {
	bell('\0');
    }
}
/*------------------------------------------------------------------*/
/*
Transpose the letter under the cursor and the letter immediately to
the left of the cursor.
*/
/*ARGSUSED*/
void
transpose_chars(ch)
    char ch;
{
    char tch;

    if ((0 < point) && (point < length))
    {
	/* first, update the display */

	backspace(line[point]);

	echo(line[point]);
	echo(line[point - 1]);

	/* now swap the chars in the line buffer */

	tch = line[point];
	line[point] = line[point - 1];
	line[point - 1] = tch;

	/* point moved forward one char */

	point++;
    }
}
/*------------------------------------------------------------------*/
/*
Delete a character at point in the line buffer. While we are at it
update the display to reflect the deletion.
*/
/*ARGSUSED*/
void
delete_char_under(ch)
    char ch;
{
    int i;

    if (point < length)
    {

	/* clear to the end of the line */

	cleartoend();

	/* retype the rest of the line */

	echoline(&line[point + 1], (length - point - 1));

	/* build the new line */

	for (i = point + 1; i < length; i++)
	{
	    line[i - 1] = line[i];
	    backspace(line[i]);
	}

	length--;

	if (point > length)
	{
	    point = length;
	}
    }

}
/*------------------------------------------------------------------*/
/*
Delete the character to the left of point in the line buffer. While we
are at it update the display to reflect the deletion.
*/
/*ARGSUSED*/
void
delete_char(ch)
    char ch;
{
    int i;

    if (point > 0)
    {
	/* move the cursor left one character */

	backspace(line[point - 1]);

	/* clear to the end of the line */

	cleartoend();

	/* retype the rest of the line */

	echoline(&line[point], (length - point));

	/* build the new line */

	for (i = point; i < length; i++)
	{
	    line[i - 1] = line[i];
	    backspace(line[i]);
	}

	length--;
	point--;
    }

}
/*------------------------------------------------------------------*/
/*
Bind the edit vector to quote_edit so that the next character
will be placed in the line buffer.
*/
/*ARGSUSED*/
void
quote(ch)
    char ch;
{
    edit = quote_edit;
}
/*------------------------------------------------------------------*/
/*
The next character will select an action from action_table[1]
*/
/*ARGSUSED*/
void
escape_1(ch)
    char ch;
{
    edit = edit_1;
}
/*------------------------------------------------------------------*/
/*
The next character will select an action from action_table[2]
*/
/*ARGSUSED*/
void
escape_2(ch)
    char ch;
{
    edit = edit_2;
}
/*------------------------------------------------------------------*/
/*
The next character will select an action from action_table[3]
*/
/*ARGSUSED*/
void
escape_3(ch)
    char ch;
{
    edit = edit_3;
}
/*------------------------------------------------------------------*/
/*
Delete the word to the left of the cursor.
*/
/*ARGSUSED*/
void
delete_word(ch)
    char ch;
{
    int i;
    int old;

    if (length > 0)
    {
	/* find the new deletion point */

	old = point;

	/* first skip over any delimiters */

	for (; (point > 0) && (delimit[line[point - 1]]); point--)
	{
	    backspace(line[point - 1]);
	}

	/* now delete until we find a delimiter */

	for (; (point > 0) && (!delimit[line[point - 1]]); point--)
	{
	    backspace(line[point - 1]);
	}

	/* clear to the end of the line */

	cleartoend();

	/* retype the rest of the line */

	echoline(&line[old], (length - old));

	/* construct the new line */

	for (i = 0; i < (length - old); i++)
	{
	    line[point + i] = line[old + i];
	    backspace(line[point + i]);
	}

	/* update the length */

	length = length - (old - point);
    }
}
/*------------------------------------------------------------------*/
/*
Go forward one word.
*/
/*ARGSUSED*/
void
forward_word(ch)
    char ch;
{
    if (length > 0)
    {
	/* first skip any delimiters */

	for (; (point < length) && (delimit[line[point]]); point++)
	{
	    echo(line[point]);
	}

	/* now skip until we find a delimiter */

	for (; (point < length) && (!delimit[line[point]]); point++)
	{
	    echo(line[point]);
	}
    }

}
/*------------------------------------------------------------------*/
/*
Lower case the word.
*/
/*ARGSUSED*/
void
lower_word(ch)
    char ch;
{
    if (length > 0)
    {
	/* first skip any delimiters */

	for (; (point < length) && (delimit[line[point]]); point++)
	{
	    echo(line[point]);
	}

	/* now skip until we find a delimiter */

	for (; (point < length) && (!delimit[line[point]]); point++)
	{
	    if ((line[point] >= 'A') && (line[point] <= 'Z'))
	    {
		line[point] = line[point] - 'A' + 'a';
		echo(line[point]);
	    }
	    else
	    {
		echo(line[point]);
	    }
	}
    }

}
/*------------------------------------------------------------------*/
/*
Upper case the word.
*/
/*ARGSUSED*/
void
upper_word(ch)
    char ch;
{
    if (length > 0)
    {
	/* first skip any delimiters */

	for (; (point < length) && (delimit[line[point]]); point++)
	{
	    echo(line[point]);
	}

	/* now skip until we find a delimiter */

	for (; (point < length) && (!delimit[line[point]]); point++)
	{
	    if ((line[point] >= 'a') && (line[point] <= 'z'))
	    {
		line[point] = line[point] - 'a' + 'A';
		echo(line[point]);
	    }
	    else
	    {
		echo(line[point]);
	    }
	}
    }

}
/*------------------------------------------------------------------*/
/*
Capitalize the word.
*/
/*ARGSUSED*/
void
capitalize_word(ch)
    char ch;
{
    if (length > 0)
    {
	/* first skip any delimiters */

	for (; (point < length) && (delimit[line[point]]); point++)
	{
	    echo(line[point]);
	}

	/* now skip until we find a delimiter */

	if ((point < length) && (!delimit[line[point]]))
	{
	    if ((line[point] >= 'a') && (line[point] <= 'z'))
	    {
		line[point] = line[point] - 'a' + 'A';
		echo(line[point]);
	    }
	    else
	    {
		echo(line[point]);
	    }
	}
	point++;

	for (; (point < length) && (!delimit[line[point]]); point++)
	{
	    if ((line[point] >= 'A') && (line[point] <= 'Z'))
	    {
		line[point] = line[point] - 'A' + 'a';
		echo(line[point]);
	    }
	    else
	    {
		echo(line[point]);
	    }
	}
    }

}
/*------------------------------------------------------------------*/
/*
Go backward one word.
*/
/*ARGSUSED*/
void
backward_word(ch)
    char ch;
{
    if (length > 0)
    {
	/* first backspace over any delimiters */

	for (; (point > 0) && (delimit[line[point - 1]]); point--)
	{
	    backspace(line[point - 1]);
	}

	/* now backspace until we find a delimiter */

	for (; (point > 0) && (!delimit[line[point - 1]]); point--)
	{
	    backspace(line[point - 1]);
	}
    }

}
/*------------------------------------------------------------------*/
/*
Move the cursor to the start of the line.
*/
/*ARGSUSED*/
void
start_of_line(ch)
    char ch;
{
    int i;

    if (length > 0)
    {
	for (i = 0; i < point; i++)
	{
	    backspace(line[i]);
	}
	point = 0;
    }
}
/*------------------------------------------------------------------*/
/*
Move the cursor one character to the left.
*/
/*ARGSUSED*/
void
backward_char(ch)
    char ch;
{
    if ((length > 0) && (point > 0))
    {
	backspace(line[point - 1]);
	point--;
    }
}
/*------------------------------------------------------------------*/
/*
Move the cursor to the right of the last character on the line.
*/
/*ARGSUSED*/
void
end_of_line(ch)
    char ch;
{
    if ((length > 0) && (point < length))
    {
	echoline(&line[point], (length - point));
	point = length;
    }
}
/*------------------------------------------------------------------*/
/*
Move the cursor one character to the right.
*/
/*ARGSUSED*/
void
forward_char(ch)
    char ch;
{
    if ((length > 0) && (point < length))
    {
	echo(line[point]);
	point++;
    }
}
/*------------------------------------------------------------------*/
/*
Add a line to the history buffer and pass it to the child process
as input.
*/
/*ARGSUSED*/
void
add_to_history(ch)
    char ch;
{
    /* Put the line in the history buffer. Make here point to the current
     * line. And increment head to point to the next history slot. */

    /* If the current line is identical to the current history line, don't
     * add it. */

    /* don't save blank lines */

    int prev;

    if ((head - 1) < 0)
    {
	prev = HISTORY_SIZE - 1;
    }
    else
    {
	prev = head - 1;
    }

    if ((length != 0) &&
	((length != hist[prev].length) ||
	    (strncmp(hist[prev].line, line, length) != 0)))
    {
	/* set the length of the entry */

	hist[head].length = length;

	/* make sure there is enough storage for the new line */

	if (hist[head].line == NULL)
	{
	    if ((hist[head].line = (char *) malloc((unsigned) length)) == NULL)
	    {
		perror("ile");
	    }
	}
	else
	{
	    if ((hist[head].line =
		    (char *) realloc(hist[head].line, (unsigned) length))
		== NULL)
	    {
		perror("ile");
	    }
	}

	(void) strncpy(hist[head].line, line, length);

	head = (head + 1) % HISTORY_SIZE;

	if (hist[head].line != NULL)
	{
	    free(hist[head].line);
	    hist[head].length = 0;
	    hist[head].line = NULL;
	}
    }

    /* reset here */

    here = head;

    /* Echo a carriage return or a newline as a cr-nl sequence. Then send the
     * line to the child process. Finally, clear the buffer for reuse. */

    fputs(ccr, stdout);
    fputs(cnl, stdout);

    line[length] = nl;
    length++;

    write_line(line, length);

    point = 0;
    length = 0;

}
/*------------------------------------------------------------------*/
/*
Erase the entire line.
*/
/*ARGSUSED*/
void
erase_line(ch)
    char ch;
{
    /* remove any text from the display */

    clearline(ch);

    /* nothing in the line buffer */

    point = 0;
    length = 0;

    /* reset here */

    here = head;

}
/*------------------------------------------------------------------*/
/*
Erase from the current cursor position to the end of the line.
*/
/*ARGSUSED*/
void
erase_to_end_of_line(ch)
    char ch;
{
    if ((length > 0) && (point < length))
    {
	cleartoend();
	length = point;
    }

}
/*------------------------------------------------------------------*/
/*
Retype the current contents of the edit buffer.
*/
/*ARGSUSED*/
void
retype_line(ch)
    char ch;
{
    int i;

    fputs(ccr, stdout);
    fputs(cnl, stdout);

    echoline(line, length);

    for (i = point; i < length; i++)
    {
	backspace(line[i]);
    }
}
/*------------------------------------------------------------------*/
/*
Go to the the next entry in the history buffer and display it.
If we are past the last history entry, then beep.
*/
void
forward_history(ch)
    char ch;
{
    if (here != head)
    {
	clearline(ch);

	here = (here + 1) % HISTORY_SIZE;
	length = hist[here].length;
	point = length;

	(void) strncpy(line, hist[here].line, length);
	echoline(line, length);
    }
    else
    {
	bell('\0');
    }
}
/*------------------------------------------------------------------*/
/*
Search backward in the history list for a line that starts with
the characters left of the cursor. If it is found make it the
current line.
*/
void
search_backward_history(ch)
    char ch;
{
    int prev;
    int i;

    /* search backward in the history */

    prev = here;

    do
    {
	prev--;

	if (prev < 0)
	{
	    prev = HISTORY_SIZE - 1;
	}
    }
    while ((hist[prev].line != NULL) &&
	(strncmp(line, hist[prev].line, point) != 0));

    /* if something was found, make it the current line */

    if (hist[prev].line != NULL)
    {
	/* remember the position in the history */

	here = prev;

	/* set the length and point correctly */

	length = hist[here].length;
	if (point > length)
	{
	    point = length;
	}

	/* redraw the line */

	clearline(ch);
	(void) strncpy(line, hist[here].line, hist[here].length);
	echoline(line, length);

	for (i = point; i < length; i++)
	{
	    backspace(line[i]);
	}
    }
    else
    {
	bell('\0');
    }
}
/*------------------------------------------------------------------*/
/*
Go back one entry in the history buffer and display it. If we are
already at the last entry, then beep.
*/
void
backward_history(ch)
    char ch;
{
    int prev;

    prev = here - 1;

    if (prev < 0)
    {
	prev = HISTORY_SIZE - 1;
    }

    if (hist[prev].line != NULL)
    {
	clearline(ch);

	here = prev;
	length = hist[here].length;
	point = length;

	(void) strncpy(line, hist[here].line, length);
	echoline(line, length);
    }
    else
    {
	bell('\0');
    }
}
/*------------------------------------------------------------------*/
/*
The following routines are utility routines used by the editing
routines.
*/
/*------------------------------------------------------------------*/
/*
Clear to the end of the current input line.
*/
void
cleartoend()
{
    /* send the clear character */

    fputs(cce, stdout);

    /* send somes nulls for padding */

    fputs("\0\0\0\0", stdout);
}
/*------------------------------------------------------------------*/
/*
Clear the input line. Backspace to the start of the line. Then clear
to the end of the line.
*/
/*ARGSUSED*/
void
clearline(ch)
    char ch;
{
    int i;

    for (i = 0; i < point; i++)
    {
	backspace(line[i]);
    }

    cleartoend();
}
/*------------------------------------------------------------------*/
/*
Echo a character. Not all characters are created equal. Control characters
are echoed in ^X form. So they take up two character positions instead of
the normal 1 character position.
*/
void
echo(ch)
    char ch;
{
    /* how should we echo the char? */

    if (ch < ' ')
    {
	(void) fputc('^', stdout);
	(void) fputc('@' + ch, stdout);
    }
    else
    {
	(void) fputc(ch, stdout);
    }
}
/*------------------------------------------------------------------*/
/*
Echo a line. Print a whole line with control characters printed in
^X form.
*/
void
echoline(line, length)
    char *line;
    int length;
{
    int i;

    for (i = 0; i < length; i++)
    {
	echo(*line++);
    }

}
/*------------------------------------------------------------------*/
/*
Backspace over a character. Generate enough bs characters to backspace
over any character.
*/
void
backspace(ch)
    char ch;
{
    if (ch < ' ')
    {
	fputs(cle, stdout);
	fputs(cle, stdout);
    }
    else
    {
	fputs(cle, stdout);
    }
}
/*------------------------------------------------------------------*/
/*
Add any character to the line buffer.
*/
void
quote_edit(ch)
    char ch;
{
    insert(ch);

    edit = edit_0;
}
/*------------------------------------------------------------------*/
/*
Given a character and an action table number either execute the
action or pass the string to (*edit)(ch)
*/
void
dispatch(table, ch)
    int table;
    char ch;
{
    char *cptr;

    switch (action_table[table][ch].flag)
    {
    case is_action:

	(*(action_table[table][ch].aors.action)) (ch);

	break;

    case is_string:

	cptr = action_table[table][ch].aors.string;
	while ((*cptr) != '\0')
	{
	    (*edit) (*cptr);
	    cptr++;
	}

	break;
    }
}
/*------------------------------------------------------------------*/
/*
Select an action from action_table[3] and execute it.
*/
void
edit_3(ch)
    char ch;
{
    /* reset so that next input is handled by edit_0 unless over ridden by
     * the action. */

    edit = edit_0;
    dispatch(3, ch);
    (void) fflush(stdout);
}
/*------------------------------------------------------------------*/
/*
Select an action from action_table[2] and execute it.
*/
void
edit_2(ch)
    char ch;
{
    /* reset so that next input is handled by edit_0 unless over ridden by
     * the action. */

    edit = edit_0;
    dispatch(2, ch);
    (void) fflush(stdout);
}
/*------------------------------------------------------------------*/
/*
Select an action from action_table[1] and execute it.
*/
void
edit_1(ch)
    char ch;
{
    /* reset so that next input is handled by edit_0 unless over ridden by
     * the action. */

    edit = edit_0;
    dispatch(1, ch);
    (void) fflush(stdout);
}
/*------------------------------------------------------------------*/
/*
Select an action from action_table[0] and execute it.
*/
void
edit_0(ch)
    char ch;
{
    dispatch(0, ch);
    (void) fflush(stdout);
}
/*------------------------------------------------------------------*/
/*
Input line editor.

Initialize the world. Then loop forever using select to wait for
characters to be available from either stdin or from master_pty.
When characters are available, pass them on after doing any needed
editing.
*/
void
ile()
{
    /* general purpose integer variable */

    int i;

    /* arguments for read and write calls */

    char buffer[BUFFER_SIZE];
    int cc;

    /* current slave_tty parameters */

    struct sgttyb slave_params;

    /* Arguments for select call */

    int nfds;
    int width;
    int readfds;

    /* what to do if the child or parent dies */

    (void) signal(SIGCHLD, clean_up);
    (void) signal(SIGSEGV, clean_up);
    (void) signal(SIGBUS, clean_up);
    (void) signal(SIGTERM, clean_up);
    (void) signal(SIGHUP, clean_up);
    (void) signal(SIGINT, clean_up);
    (void) signal(SIGQUIT, clean_up);

    /* what to do it the window changes size */

    (void) signal(SIGWINCH, change_window);

    /* Get all the different pieces of the current ttys' state and copy them
     * to the slave_tty. */

    /* tty sgttyb */

    (void) ioctl(READ, TIOCGETP, &tty_sgttyb);
    (void) ioctl(slave_tty, TIOCSETP, &tty_sgttyb);

    /* tty line discipline */

    (void) ioctl(READ, TIOCGETD, &tty_ldisc);
    (void) ioctl(slave_tty, TIOCSETD, &tty_ldisc);

    /* tty tchars */

    (void) ioctl(READ, TIOCGETC, &tty_tchars);
    (void) ioctl(slave_tty, TIOCSETC, &tty_tchars);

    /* tty mode */

    (void) ioctl(READ, TIOCLGET, &tty_mode);
    (void) ioctl(slave_tty, TIOCLSET, &tty_mode);

    /* tty ltchars */

    (void) ioctl(READ, TIOCGLTC, &tty_ltchars);
    (void) ioctl(slave_tty, TIOCSLTC, &tty_ltchars);

    /* tty windsize */

    (void) ioctl(READ, TIOCGWINSZ, &tty_winsize);
    (void) ioctl(slave_tty, TIOCSWINSZ, &tty_winsize);
    windowchanged = FALSE;

    /* "login" the user */

    (void) setutmp(slave_tty, TRUE);

    /* set raw mode on tty */
    {
	struct sgttyb params;
	struct tchars tparams;
	struct ltchars ltparams;

	/* Simulate RAW but allow original parity to work.  Thus we use
	 * CBREAK with all the options turned off. */

	params = tty_sgttyb;
	params.sg_flags = CBREAK;
	(void) ioctl(READ, TIOCSETP, &params);

	tparams = tty_tchars;
	tparams.t_intrc = -1;
	tparams.t_quitc = -1;
	tparams.t_startc = -1;
	tparams.t_stopc = -1;
	tparams.t_eofc = -1;
	tparams.t_brkc = -1;
	(void) ioctl(READ, TIOCSETC, &tparams);

	ltparams = tty_ltchars;
	ltparams.t_suspc = -1;
	ltparams.t_dsuspc = -1;
	ltparams.t_rprntc = -1;
	ltparams.t_flushc = -1;
	ltparams.t_lnextc = -1;
	(void) ioctl(READ, TIOCSLTC, &ltparams);
    }

    /* set new mode on tty */

    {
	int mode;
	mode = LNOFLSH | LDECCTQ | LLITOUT;
	(void) ioctl(READ, TIOCLSET, &mode);
    }

    /* get descriptor table size */

    width = getdtablesize();
    if (width > 32)
    {
	width = 32;
    }

    /* set initial edit function */

    edit = edit_0;

    /* initialize line buffer */

    point = 0;
    length = 0;

    /* initialize history buffer */

    head = 0;
    here = 0;

    for (i = 0; i < HISTORY_SIZE; i++)
    {
	hist[i].length = 0;
	hist[i].line = NULL;
    }

    for (;;)
    {
	readfds = (1 << READ) | (1 << master_pty);

	/* wait for input from stdin or master_pty */

	nfds = select(width,
	    (fd_set *) & readfds,
	    (fd_set *) NULL,
	    (fd_set *) NULL,
	    (struct timeval *) NULL);

	if (nfds == -1)		/* an exception has occured */
	{
	    if (windowchanged)
	    {
		/* nothing serious, the window changed size */

		windowchanged = FALSE;
	    }
	    else
	    {
		perror("ile");
		clean_up();
	    }
	}
	else if ((nfds > 0) && (readfds != 0))	/* something to read */
	{
	    if ((readfds & (1 << master_pty)) != 0)
	    {
		/* read the pending characters. */

		cc = read(master_pty, buffer, BUFFER_SIZE);

		/* display the characters. */

		(void) write(WRITE, buffer, cc);
	    }

	    if ((readfds & (1 << READ)) != 0)
	    {
		/* read the pending characters. */

		cc = read(READ, buffer, BUFFER_SIZE);

		/* if the slave is in RAW or CBREAK mode, or has turned off
		 * ECHO then we should not mess with its input characters */

		(void) ioctl(slave_tty, TIOCGETP, &slave_params);

		if (((slave_params.sg_flags & (RAW | CBREAK)) != 0) ||
		    (slave_params.sg_flags & ECHO) == 0)
		{
		    edit = pass;
		}
		else if (edit == (void (*) ()) pass)
		{
		    edit = edit_0;
		}

		/* decide what to do with the characters. */

		for (i = 0; i < cc; i++)
		{
		    (*edit) (CHAR_MASK & buffer[i]);
		}
	    }

	}
    }

}
/*------------------------------------------------------------------*/
/*
The child process.

Make the pty the processes controling terminal. Bind the pty to
stdin, stdout, and stderr. Then exec the users program.
*/
void
child(argv)
    char *argv[];
{
    /* shell name pointers */

    char *shellname;
    char *shellpath;
    char *dashshellname;

    /* close all file descriptors */

    (void) close(READ);
    (void) close(WRITE);
    (void) close(ERROR);
    (void) close(slave_tty);

    /* get rid of controlling terminal */

    {
	int tty;

	if ((tty = open("/dev/tty", O_RDWR) == -1) ||
	    (ioctl(0, TIOCNOTTY, 0) == -1) ||
	    (close(tty) == -1))
	{
	    perror("ile");
	}
    }

    /* open the tty again */
    /* this makes the pty the controlling terminal */

    if ((slave_tty = open(ttydev, O_RDWR)) == -1)
    {
	perror("ile");
    }

    /* slave_tty is now stdin */

    /* bind slave_tty to stdout */

    (void) dup2(slave_tty, WRITE);

    /* bind slave_tty to stderr */

    (void) dup2(slave_tty, ERROR);

    /* close master_pty descriptor */

    (void) close(master_pty);

    /* Fire up application program. If no program name is given then fire up
     * users favorite shell. */

    /* get the name of the users shell. default to /bin/csh */

    if ((shellpath = getenv("SHELL")) == NULL ||
	(*shellpath == '\0'))
    {
	shellpath = "/bin/csh";
    }

    /* get just the name */

    if ((shellname = strrchr(shellpath, '/')) != NULL)
    {
	shellname += sizeof(char);
    }

    /* if the current argv[0] starts with -, then the new argv[0] */
    /* should start with - */

    if (*(argv[0]) == '-')
    {
	dashshellname = (char *) malloc((unsigned) strlen(shellname) + 2);
	(void) strcpy(dashshellname, "-");
	(void) strcat(dashshellname, shellpath);
    }
    else
    {
	dashshellname = shellname;
    }

    /* execute the shell or the specified program */

    if (argv[1] == NULL)
    {
	/* execute default shell */

	execlp(shellpath, dashshellname, 0);
    }
    else if (*argv[1] == '-')
    {
	/* there is an initialization file */

	if (argv[2] == NULL)
	{
	    /* execute default shell */

	    execlp(shellpath, dashshellname, 0);
	}
	else
	{
	    /* execute specified program */

	    execvp(argv[2], &argv[2]);
	}
    }
    else
    {
	/* execute specified program */

	execvp(argv[1], &argv[1]);
    }

    /* this executes if exec fails */

    perror("ile");
    exit(1);
    /* NOTREACHED */
}
/*------------------------------------------------------------------*/
/*
Set up default key bindings and delimeters.
*/
void
default_bindings()
{
    int i;

    /* clear delimiter vector and the action table */

    for (i = 0; i < CHAR_SET_SIZE; i++)
    {
	delimit[i] = FALSE;

	action_table[0][i].aors.action = insert;
	action_table[1][i].aors.action = bell;
	action_table[2][i].aors.action = bell;
	action_table[3][i].aors.action = bell;

	action_table[0][i].flag = is_action;
	action_table[1][i].flag = is_action;
	action_table[2][i].flag = is_action;
	action_table[3][i].flag = is_action;
    }

    /* default delimiters */

    delimit[' '] = TRUE;	/* blank */
    delimit['/'] = TRUE;	/* slash */
    delimit['.'] = TRUE;	/* dot */
    delimit['-'] = TRUE;	/* dash */

    /* default action_table[0] */

    action_table[0][CA].aors.action = start_of_line;
    action_table[0][CB].aors.action = backward_char;
    action_table[0][CE].aors.action = end_of_line;
    action_table[0][CF].aors.action = forward_char;
    action_table[0][CK].aors.action = erase_to_end_of_line;
    action_table[0][CU].aors.action = erase_line;
    action_table[0][CL].aors.action = retype_line;
    action_table[0][CN].aors.action = forward_history;
    action_table[0][CP].aors.action = backward_history;
    action_table[0][CR].aors.action = search_backward_history;
    action_table[0][CT].aors.action = transpose_chars;
    action_table[0][CV].aors.action = quote;
    action_table[0][del].aors.action = delete_char;
    action_table[0][esc].aors.action = escape_1;
    action_table[0][cr].aors.action = add_to_history;
    action_table[0][nl].aors.action = add_to_history;
    action_table[0][CX].aors.action = delete_char_under;

    action_table[0][CC].aors.action = pass;
    action_table[0][CD].aors.action = pass;
    action_table[0][CQ].aors.action = pass;
    action_table[0][CS].aors.action = pass;
    action_table[0][CZ].aors.action = pass;

    /* default action_table[1] ^[ c */

    action_table[1]['b'].aors.action = backward_word;
    action_table[1]['f'].aors.action = forward_word;
    action_table[1][del].aors.action = delete_word;
    action_table[1]['u'].aors.action = upper_word;
    action_table[1]['l'].aors.action = lower_word;
    action_table[1]['c'].aors.action = capitalize_word;
    action_table[1]['['].aors.action = escape_2;
    action_table[1][esc].aors.action = complete_file;
    action_table[1]['s'].aors.action = complete_file_full;
    action_table[1]['d'].aors.action = show_files;
    action_table[1]['p'].aors.action = query_path;

    /* default action_table[2] ^[ [ */

    action_table[2]['A'].aors.action = backward_history;
    action_table[2]['B'].aors.action = forward_history;
    action_table[2]['C'].aors.action = forward_char;
    action_table[2]['D'].aors.action = backward_char;

}
/*------------------------------------------------------------------*/
/*
Return a character or EOF. This routine reads characters from input
and converts them into a character using the following rules.

The character may be a single character, a control
character indicated by ^x, an octal number starting with \, or an
escaped character indictated by \x.
*/
int
scan_char(input)
    FILE *input;
{
    int ch;
    int value;

    ch = fgetc(input);
    switch (ch)
    {
    case '^':

	/* it is a control character */

	for (ch = fgetc(input); '@' <= ch; ch = ch - '@');

	break;

    case '\\':

	/* octal or an escaped character? */

	ch = fgetc(input);
	if (('0' <= ch) && (ch <= '7'))
	{

	    /* its an octal number */

	    value = 0;
	    while (('0' <= ch) && (ch <= '7'))
	    {
		value = (value * 8) + (ch - '0');
		ch = fgetc(input);
	    }
	    (void) ungetc(ch, input);

	    ch = value & 0177;	/* make sure it is in range */
	}
	else
	{
	    /* its an escaped character */

	    ch = fgetc(input);
	}

	break;

    case '\n':

	/* the real end of the line */

	ch = EOL;

	break;

    default:

	/* it is just itself */

	break;
    }

    return (ch);

}
/*------------------------------------------------------------------*/
/*
Set key bindings and delimiters from the users file.
*/
void
user_bindings(file)
    FILE *file;
{

#define NAME_SIZE 40

    static struct action_name_table
    {
	char *name;
	void (*action) ();
    } action_name_table[] =
    {
	{
	    "complete_file_full", complete_file_full
	},
	{
	    "complete_file", complete_file
	},
	{
	    "show_files", show_files
	},
	{
	    "query_path", query_path
	},
	{
	    "bell", bell
	},
	{
	    "pass", pass
	},
	{
	    "insert", insert
	},
	{
	    "transpose_chars", transpose_chars
	},
	{
	    "delete_char", delete_char
	},
	{
	    "delete_char_under", delete_char_under
	},
	{
	    "quote", quote
	},
	{
	    "escape_1", escape_1
	},
	{
	    "escape_2", escape_2
	},
	{
	    "escape_3", escape_3
	},
	{
	    "delete_word", delete_word
	},
	{
	    "upper_word", upper_word
	},
	{
	    "lower_word", lower_word
	},
	{
	    "capitalize_word", capitalize_word
	},
	{
	    "forward_word", forward_word
	},
	{
	    "backward_word", backward_word
	},
	{
	    "start_of_line", start_of_line
	},
	{
	    "backward_char", backward_char
	},
	{
	    "end_of_line", end_of_line
	},
	{
	    "forward_char", forward_char
	},
	{
	    "add_to_history", add_to_history
	},
	{
	    "erase_line", erase_line
	},
	{
	    "erase_to_end_of_line", erase_to_end_of_line
	},
	{
	    "retype_line", retype_line
	},
	{
	    "forward_history", forward_history
	},
	{
	    "backward_history", backward_history
	},
	{
	    "search_backward_history", search_backward_history
	},
	{
	    "", NULL
	}
    };

    char name[NAME_SIZE];

    int ch;
    int i;

    int linecount;
    int table;
    int entry;

    /* First clear the default delimiters */

    for (i = 0; i < CHAR_SET_SIZE; i++)
    {
	delimit[i] = FALSE;
    }

    /* Now read the delimiter characters */

    while (((int) (ch = fgetc(file)) != EOF) && (ch != '\n'))
    {
	delimit[ch] = TRUE;
    }

    linecount = 2;

    /* Now read the character binding pairs */

    while ((int) (ch = fgetc(file)) != EOF)
    {
	switch (ch)
	{
	case '\n':

	    /* skipping a blank line */
	    linecount++;

	    break;

	case '0':
	case '1':
	case '2':
	case '3':

	    /* which table is this entry directed to? */

	    table = ch - '0';

	    /* get the character code */

	    entry = scan_char(file);

	    /* make sure the '=' is there */

	    ch = fgetc(file);
	    if (ch != '=')
	    {
		(void) fprintf(stderr,
		    "ile: '=' missing on line %d\n",
		    linecount);
		exit(1);
		/* NOTREACHED */
	    }

	    /* collect the action name or string */

	    for (ch = scan_char(file), i = 0;
		((int) ch != EOL) && (i < (NAME_SIZE - 1));
		ch = scan_char(file), i++)
	    {
		name[i] = ch;
		name[i + 1] = '\0';
	    }

	    /* look it up in the action_name_table */

	    for (i = 0;
		(action_name_table[i].action != NULL) &&
		(strcmp(name, action_name_table[i].name) != 0);
		i++);

	    /* if it was found, put it in the action array */

	    if (action_name_table[i].action == NULL)
	    {
		/* must be a string */

		action_table[table][entry].flag = is_string;
		action_table[table][entry].aors.string =
		    (char *) malloc((unsigned) strlen(name) + 1);
		(void) strcpy(action_table[table][entry].aors.string, name);
	    }
	    else
	    {
		/* its an action */

		action_table[table][entry].flag = is_action;
		action_table[table][entry].aors.action =
		    action_name_table[i].action;
	    }

	    linecount++;	/* count the line */

	    break;

	default:
	    (void) fprintf(stderr,
		"\nile: error in initialization file on line %d\n",
		linecount);
	    exit(1);
	    /* NOTREACHED */
	}
    }

    (void) fclose(file);
}
/*------------------------------------------------------------------*/
/*
Initialize key bindings and delimiters.
*/
void
initialize(argv)
    char *argv[];
{
    FILE *file;
    char name[BUFFER_SIZE];
    char *pwd;

    /* set up the default bindings */

    default_bindings();

    /* Look for an initialization file. If it's there, load it. */

    name[0] = '\0';
    homedir = getenv("HOME");
    if (homedir == NULL)
    {
	/* no home dir, use / instead */
	name[0] = '\0';
    }
    else
    {
	(void) strcpy(name, homedir);
    }
    (void) strcat(name, "/.ilerc");

    /* initialize currentdir */

    pwd = getenv("PWD");
    if (pwd == NULL)
    {
	/* no pwd, use homedir instead */
	(void) strcpy(currentdir, homedir);
    }
    else
    {
	(void) strcpy(currentdir, pwd);
    }

    if ((argv[1] != NULL) &&
	(*argv[1] == '-') &&
	((file = fopen(argv[1] + 1, "r")) != NULL))
    {
	/* load the users bindings */

	user_bindings(file);
    }
    else if (((file = fopen("./.ilerc", "r")) != NULL) ||
	((file = fopen(name, "r")) != NULL))
    {
	user_bindings(file);
    }
}
/*------------------------------------------------------------------*/
/*
*/
/*ARGSUSED*/
main(argc, argv)
    int argc;
    char *argv[];
{

    /* Child process id */

    int childpid;

    /* identify yourself */

    (void) fprintf(stdout, "ile rev.2\n\r");

    /* create the tty/pty pair */

    getpty(&master_pty, &slave_tty);

    /* get control sequences from termcap */

    get_termcap();

    /* initialize the dispatch vectors */

    initialize(argv);

    /* create the child process */

    childpid = fork();

    switch (childpid)
    {
    case 0:			/* child process */

	child(argv);
	break;

    case -1:			/* fork failed */

	perror("ile");
	exit(1);
	/* NOTREACHED */

    default:			/* parent process */
	ile();
	break;
    }

}
