#define VERSION_ID "(version 1.7a [21-jun-91])"

/* The following section will be extracted by the doc program in a .doc file */
/* The miriad package needs such files, although in this particular case     */
/* the miriad.doc file itself it not really needed.			     */

/*= miriad - general program to run miriad programs */
/*& pjt */
/*: tools */
/*+* /
  Miriad is a command-line interface to run Miriad tasks.

  Some important Miriad commands are:

      inp [taskname] ... preview input parmeters to task taskname
      parm=value ....... set task parameter parm to value
      unset parm ....... unset task parameter parm
      help [taskname] .. show on-line docs for task taskname
      go [taskname] .... run task taskname
      view [taskname] .. use editor to set parms for task taskname
      task taskname .... set default for taskname [taskname]
      quit ............. exit immediately
      exit ............. exit, but save all parm settings
      ? ................ overview of all miriad commands
      help ? ........... listing of all miriad tasks

  The directory MIRBIN contains all miriad tasks, whereas the directory
  MIRPDOC must contain the appropriate doc files (help files with directives 
  on available keywords). The invocation 'miriad -?' (or -\? if a shell
  would use the ?) will display all command line options.

  Any unrecognized commands are passed on to the parent shell for execution.
  Thus, the command "ls -l $HOME" will return a listing of the files
  in the user's home directory on UNIX. UNIX aliases are also understood
  for a number of shells. In VMS the command "DIR SYS$LIBRARY" would
  equally well work.
/*--                                                                    */
/************************************************************************/
/* If you need to toy more seriously with miriad.c - a corresponding    */
/* Makefile to it also belongs to your tools, as well as a number of    */
/* test files. These are not delivered with miriad, but can be found n  */
/* within the NEMO environment. See current maintainer for more details.*/
/************************************************************************/
/*  COMPILE DEFINE OPTIONS:                                             */
/*                                                                      */
/*  Use as:  -Doption (unix) or /DEFINE=(option=1) (VMS) in CC          */
/*                                                                      */
/*   <option>      <explanation>                                        */
/*   --------      -------------                                        */
/*  PATHSEARCH     if set, full $PATH is searched instead of $MIRBIN    */
/*                 not usable in VMS                                    */
/*  INTERRUPT      if set, interrupts like ^\, ^C, ^Y are caught        */
/*  GETENV         set this if your OS has no char *getenv()            */
/*  NEWENV         if variables can be added to **menviron               */
/*  EXECL          if foreign commands can use execl() to get at aliases*/
/*  MULTISET       if 'set key val1 val2 val3' is allowed               */
/*		   [not sure if this macro will survive - should be default] */
/*  TCL            use TCL shell? (testing by PJT-nobody should use it) */
/*  READLINE       use the readline GNU command line editor?            */
/*		   this not only gives BASH command line editing, but   */
/*		   also file completion					*/
/*                                                                      */
/************************************************************************/

#if defined(sun)
# define PATHSEARCH 1
# define INTERRUPT  1
# define NEWENV     1
# define EXECL      1
# define MULTISET   1
#endif

#if defined(vms)
# define INTERRUPT  1
# define NEWENV     1
#endif

#if defined(READLINE) && defined(TCL)
Cannot install Tcl and Readline together
#endif

/************************************************************************/
/*                                                                      */
/*  History:                                                            */
/*    rjs  Dark-ages Original version.                                  */
/*    pjt   4dec89   Warning when no write permission in save_vars.     */
/*    rjs  21feb90   Minor enhancements, suggested by Brian Glendenning.*/
/*    pjt  26feb90   Changed environment variables, more messages       */
/*    rjs   5mar90   Merged PJT and RJS versions.                       */
/*    pjt  12mar90   don't write keyfile when not needed                */
/*    pjt  15mar90   'gob' is same as 'go' with backgrounding           */
/*    pjt  16mar90   added save, and help with no options               */
/*    pjt   9apr90   some more help                                     */
/*    rjs  26apr90   Looks in local directory for .doc files. On VMS,   */
/*                   it checks for the foreign command definition,      */
/*                   before overwriting it with its own.                */
/*    pjt   6may90   compile option to search $PATH in Unix (execvp)    */
/*    pjt  13may90   -b BIN -d DEF -p PDOC options                      */
/*    pjt  15jun90   added TASK command - to set default task - setenv  */
/*    pjt  20jun90   quit/exit is now different - load/save have def    */
/*    pjt  10jul90   catch a few signals                                */
/*    pjt   6aug90   new INPUT command                                  */
/*    pjt  21aug90   Testing TCL - for Version 2.x                      */
/*                   added dounset to get_vars                          */
/*    pjt   6sep90   some TCL cleanup - lastexit with -s switch         */
/*                   reset keywords cleans ALL variables                */
/*    pjt  18sep90   cd command                                         */
/*    pjt  15nov90   local shell commands must use ! now                */
/*                   doset() can now have multiple values               */
/*    pjt  19dec90   removed complaint line when shell command used     */
/*    pjt  23jan91   introduced MIRPAGER env.var. (needed for doc)      */
/*    pjt  22feb91   4.1c fixed erroneous switchings of taskname -      */
/*         27feb91   4.1d ALIAS command - added to SET command          */
/*          5apr91   4.2 merged Sanjay's code (unreleased)              */
/*          6may91       -d implemented                                 */
/*    rjs  22may91   Fixed bug in doset()                               */
/*    pjt  22may91   Fixed new VMS problems - set interrupt on in VMS   */
/*         26may91   child signalling from Sanjay added                 */
/*    pjt  28may91   V1.6 added BLANK code for setting variables        */
/*         12jun91   V1.7 added optional GNU's READLINE library         */
/*	   21jun91        renamed environ to menviron for cray		*/
/*                                                                      */
/*    ToDo anyhow:                                                      */
/*      check earlier if lastexit can be written, otherwise complain    */
/*      aliases do not support strings with blanks - they are tokenized */
/*    PJT's ideas:                                                      */
/*      - redirect stdout/stderr to something useful, or:               */
/*      - tee stdout/stderr to a general miriad.log like lastexit       */
/*        [some work just before the execv(P) is needed?                */
/*      - cursor history                                                */
/*      - allow 'go histo in=newfile' and                               */
/*        'go in=newfile' if histo was alread the default????           */
/*      - allow a PATH in -d, -b                                        */
/*    LGM's ideas:                                                      */
/*      - also remember per keyword which program used it, such that    */
/*        they can be used next to each other                           */
/*        (like TGET and TPUT in AIPS)                                  */
/*     Sanjay Bhatnagar idea's                                          */
/*      - expand any regular expression in key values                   */
/*        (for e.g.in=~sanjay/test/test.dat instead of in=/usr2/sanjay  */
/*         /test/test.dat)                                              */
/*                                                                      */
/************************************************************************/

#include <stdio.h>
#if defined(INTERRUPT)
#include <signal.h>
#endif

/*   These come from SYS5 versions which do not always have them def'd */
/*   I need them for 3b1 (mx68k)                                       */
/*   This is for access(2) - enventually POSIX will solve this kludge  */
#if !defined(R_OK)
#define R_OK 4
#endif
#if !defined(X_OK)
#define X_OK 1
#endif

#define MAXARGS   64
#define MAXBUF   256
#define HASHSIZE 127
#define MAXINPUT  10

typedef struct variable {        /* structure def for miriad variables */
    char *name;
    char *value; 
    struct variable *fwd; 
} VARIABLE;

VARIABLE *hashtable[HASHSIZE];

struct {
    char *name;
    char *value;
} args[MAXARGS];

char names[MAXBUF]; 
char buffer[MAXBUF];    /* move over ??? */
char taskname[MAXBUF];          /* current default taskname */
char *lastexit="lastexit";         /* default name of lastexit file */

int Qkeys = 0;      /* 0: keys were not updated     1: were */
int debug_level = 0;        /* 0:  initially no debug output */
int input_level = 0;        /* nesting level of INPUT command */
FILE *fpinput[MAXINPUT];    /* open files for INPUT command */
char *bindir = NULL;        /* -b: come from **menviron */
char *docdir = NULL;        /* -d: come from **menviron */
char *defdir = NULL;        /* -p: come from **menviron */


#if defined(TCL)
# include <tcl.h>
  Tcl_Interp *tcl_interp;
  Tcl_CmdBuf tcl_buffer;
  int cmdCd(), cmdTest(), cmdEcho(), dP();
#endif

/* forward references to make (ansi) compilers happy */

char *getenv();
void get_vars(),save_vars(),doset(),dounset(),doinp(),dogo(),dohelp(),doq(),
     dosetenv(), dotask(), dot(), doversion(), doinput(),dotcl(),doreset(),
     doload(), dosave(), docommand(),doview(),filename(), bug(), docd();
int  getline(), task_args();
#ifdef unicos
  void unicos_check();
#endif
#if defined(INTERRUPT)
  void review();
#endif

#if defined(_trace_)
extern int errno;       /* mostly <errno.h> would do the job */
#else
#include <errno.h>
#endif

char **menviron;  /* point to environment - not used in vms */
char *shell=NULL;   /* remember which shell user uses */

char *malloc(), *strchr();

/************************************************************************/
main(ac,av,ep)
int ac;         /* number of arguments + 1 */
char *av[];     /* progname + arguments */
char *ep[];     /* environment strings ENV=val */
{
  int argc;
  char *argv[MAXARGS];

  ini_miriad(ac,av,ep);                /* initialize lot's of stuff        */

  for(;;) {                         /* infinite loop                    */
    argc = getline(argv,NULL);      /* get a new line; and tokenize it  */
    if (argc)                       /* if any arguments left            */
        do_miriad(argc,argv);       /* ... call dispatcher to execute   */
  }                                 /* end_miriad() is called inside    */
}

ini_miriad(ac,av,ep)
int ac;
char *av[];
char *ep[];
{
  int i;
  char name[MAXBUF], *s, *startup[2];

  printf("Miriad shell %s\n",VERSION_ID);              /* Hello world ! */
#if defined(INTERRUPT)
  if (debug_level) printf("Catching signals SIGTERM, SIGQUIT, SIGINT\n");
  signal(SIGTERM, review);            /* catch interrupts */
  signal(SIGQUIT, review);            /* for review */
  signal(SIGINT,  review);            /* and ^C also */ 
#endif
  menviron = ep;                     /* set environment for global access */
  set_defaults(ac,av);              /* (re)set defaults from arguments */
  taskname[0] = '\0';               /* make sure this is NULL at start */

  for(i=0; i<HASHSIZE; i++) {        /* Initialise the hash table */
        hashtable[i] = NULL;
  }

  get_vars(lastexit);                /* Read "lastexit" file */

                                /* set name for miriad startup aliases */
#ifdef VMS 
  filename(name,"SYS$LOGIN","miriad",".rc");  /* HOME is also in **menviron */
#else
  filename(name,"HOME","",".miriadrc");
#endif
  if (access(name,R_OK)==0) {           /* if startup file indeed found */
    printf("[Reading %s]\n",name);
    startup[0] = "input";               /* build cmdline args */
    startup[1] = name;
    do_miriad(2,startup);               /* and execute */
  }

#if defined(TCL)
  install_tcl();
#endif

#if defined(READLINE)
  ini_readline();
#endif
}


end_miriad()
{
  if(input_level>0) {            /* if exit was from an input file */
      input_level--;                /* decrease stack of input files */
      fclose(fpinput[input_level]); /* close that file */
      return;                       /* and keep on trucking */
  }

  if (Qkeys)
      save_vars(lastexit);       /* Save all the parameters in lastexit. */
  else
      printf("### Warning: Variables not saved in %s\n",lastexit);

  exit(0);
}

/************************************************************************/
int getline(argv,cmdnow)
char *argv[], *cmdnow;
/*
  This prompts and reads a line from STDIN. It breaks it into tokens.
  If the second token is an equals sign, then this makes it into a
  "set" command.
  If TCL is compiled into miriad, assembling the commandline and tokening
  is mostly done by TCL now. Without TCL a commandline MUST be on one line

  If 'cmdnow' is not a NULL pointer, this is taken as a command and parsed
  tokenized and returned. In case TCL is used, no multiple command lines
  are allowed for now.

  Problem:  \n in read line in TCL vs. non-TCL
------------------------------------------------------------------------*/
{
  int n,inter,doset,i,gotPartial, result;
  char *s, *cmd;
  char buffer2[MAXBUF], *line, prompt[40];

  if (cmdnow==NULL) {       /* get new command from stdin or a file */

  gotPartial = 1;  /* needed if we allow multiple input lines before parsing */
  do {                  /* reads from stdin/file until line complete */
    if (input_level==0) {		/* if interactive */
      if (*taskname)
        sprintf(prompt,"%s%% ",taskname);
      else
        sprintf(prompt,"Miriad%% ");
#if !defined(READLINE)
      printf("%s",prompt);
      clearerr(stdin);          /* clear junk from buffer ?? */
      if(fgets(buffer,MAXBUF-1,stdin) == NULL) 
        strcpy(buffer,"exit");
#else
      for(;;) {                 /* read until we have something */
        char *readline();

        line = readline(prompt);
        if (line) {
            stripwhite(line);    /* get rid of leading and trailing blank */ 
            if (*line) {
                strcpy(buffer,line);
                break;
            }
         }
       }
       if (line) free(line);
#endif
    } else {				/* or reading from input file */
      if(fgets(buffer,MAXBUF-1,fpinput[input_level-1])==NULL) {
        input_level--;
        fclose(fpinput[input_level]);
        return(0);
      } 
    }
#if defined(TCL)
                /* TCL needs buffer to be \n terminated */
    cmd = Tcl_AssembleCmd(tcl_buffer,buffer); 
    if (cmd)
        gotPartial = 0;         /* got it, can try and execute */
#else
                /* we'll pass a NULL terminated string back */
    if (buffer[strlen(buffer)-1] == '\n')  /* need to overwrite '\n' ? */
        buffer[strlen(buffer)-1] = '\0';                   /* do it */
    gotPartial = 0;             /* got it, can try and execute */
#endif
  } while (gotPartial);

  } else {              /* input from an argument */
    printf("[Executing: %s]\n",cmdnow);
    cmd = cmdnow;
    strcpy(buffer,cmdnow);
  }

#if defined(TCL)
  result = Tcl_RecordAndEval(tcl_interp, cmd, 0);
  if (result == TCL_OK) {
    if(*tcl_interp->result != 0)
        printf("%s\n",tcl_interp->result);
  } else {
    if (result == TCL_ERROR) {
        printf("Error");
    } else {
        printf("Error %d", result);
    }
    if (*tcl_interp->result != 0) {
        printf(": %s\n", tcl_interp->result);
    } else {
        printf("\n");
    }
  }

  return(0);

/* 
        what needs to be done here, is testing if user wanted
        to use assignment, ala 'x=b', i.e. first non-alphanum is
        a '=', or if first char is an '!', in which case
        the buffer+1 is passed to system(2)
        If none of these are true, TCL mode can return with 0
        to fake that it has to do nothing in do_miriad() and
        nothing is returned in the output argv[] array.
*/
        
#endif

#if defined(READLINE)
  add_history(buffer);
#endif  
  s = buffer;
  inter = 1;
  n = 0;
  doset = 0;
  while(*s){
    if(*s == ' ' || *s == '\t'){
      inter = 1;
      *s = 0;
    } else if(*s == '=' && n == 1){
      doset = 1;
      inter = 1;
      *s = 0;
    } else if(inter){
      argv[n++] = s;
      inter = 0;
    }
    s++;
  }

/* If it was a set command, shift everything down by one, add the "set"
   to the top of the list, and increase the count. */

  if(doset){
    for(i=n; i>0; i--) argv[i] = argv[i-1];
    n++;
    argv[0] = "set";
  }
  return(n);
}
/************************************************************************/
/* Strip whitespace from the start and end of STRING. */

#define whitespace(c)    ((c) == ' ' || (c) == '\t' || (c) == '\n')

stripwhite (string)
char *string; 
{
  register int i = 0;
        
  while (whitespace (string[i]))
    i++;
	           
  if (i)
    strcpy (string, string + i);
	     
  i = strlen (string) - 1;
         
  while (i > 0 && whitespace (string[i]))
    i--;
				           
  string[++i] = '\0';
}

/************************************************************************/
void doset(argc,argv)
int argc;
char *argv[];
/*
------------------------------------------------------------------------*/
{
  char *t, cat[MAXBUF]; 
  int hashval,found,i,len;
  VARIABLE *v;

/* Check the arguments. */

  if(argc == 1) {                      /* print all values - and return */
    for(i=0;i<HASHSIZE;i++) {
        v = hashtable[i];
        while(v) {
            printf("%8s = %s\n",v->name,v->value);
            v = v->fwd;
        }
    }
    return;
  } 

/* Find the value of the parameter, stored in the hash table. */

  hashval = 0;
  t = argv[1];      /* name of the variable */
  while(*t) hashval += *t++;
  v = hashtable[hashval % HASHSIZE];
  while(v != NULL){
    if(!strcmp(argv[1],v->name)) break;
    v = v->fwd;
  }

  if (argc == 2) {          /* print value of one variable - and return */
    if (v)
        printf("%8s = %s\n",v->name,v->value);
    else
        printf("Variable %s not been set yet\n",argv[1]);
    return;
  }

/* Create the variable if needed. Fill in the value. */

  if(v == NULL){
    v = (VARIABLE *)malloc(sizeof(VARIABLE));
    v->fwd = hashtable[hashval % HASHSIZE];
    hashtable[hashval % HASHSIZE] = v;
    v->name = (char *)malloc(strlen(argv[1])+1);
    strcpy(v->name,argv[1]);
  } else {
    free(v->value);
  }
#if defined(MULTISET)
  if (argc>2) {	  		/* assemble all remaining args in the 3rd */
    strcpy(cat,argv[2]);
    for (i=3; i<argc; i++) {
        strcat(cat," ");
        strcat(cat,argv[i]);
    }
    argv[2] = cat;
    argc = 3;
  }
#endif
  len = 1;                          /* figure out how long value is */
  for (i=2; i<argc; i++)
    len += strlen(argv[i])+1;
  v->value = (char *)malloc(len);   /* allocate space for value */
  strcpy(v->value,argv[2]);         /* copy first one */
  for (i=3; i<argc; i++) {          /* and catenate all other ones */
    strcat(v->value," ");
    strcat(v->value,argv[i]);
  }
}


/************************************************************************/
void doreset(argc,argv)
int argc;
char *argv[];
/*
 * This routine needs to be upgraded: when the argument is 'all' it can
 * result the whole hashtable, but when no argument is given, or the
 * name of an existing executable, it should only clear the variables
 * for that executable...
------------------------------------------------------------------------*/
{
  int i;

  save_vars("lastexit.bck");
  printf("[All miriad keywords have been blanked - backup in lastexit.bck]\n");
  printf("[Use 'load <file>' to load them back]\n");

  for(i=0; i<HASHSIZE; i++)         /* Reset HASH table */
        hashtable[i] = NULL;        /* BUG: It doesn't free used memory */
}
/************************************************************************/
void dounset(argc,argv)
int argc;
char *argv[];
/*
------------------------------------------------------------------------*/
{
  int i,hashval;
  char *t;
  VARIABLE *v,*prev;

  for(i=1; i < argc; i++){  /* handle all arguments as keys to be unset */
    t = argv[i];
    hashval = 0;
    while(*t) hashval += *t++;
    prev = NULL;
    v = hashtable[hashval % HASHSIZE];
    while(v != NULL){
      if(!strcmp(argv[i],v->name)) break;
      prev = v;
      v = v->fwd;
    }
    if(v != NULL) {
      printf("[unset %s]\n",argv[i]);
      if(prev == NULL) hashtable[hashval % HASHSIZE] = v->fwd;
      else prev->fwd = v->fwd;
      free(v->name);
      free(v->value);
      free((char *)v);
    } /* else: variable not found, so ignore */
  }
}

/************************************************************************/
void doinp(argc,argv)
int argc;
char *argv[];
/*
------------------------------------------------------------------------*/
{
  int i,n;
  char *testname;

  if(argc > 2) fprintf(stderr,"inp: Extra arguments on line ignored.\n");
  testname = (argc > 1 ? argv[1] : taskname);
  n = task_args(testname);
  if(n < 0){
    fprintf(stderr,"inp: Found no documentation on task %s.\n",testname);
  } else if( n == 0){
    if (argc>1) strcpy(taskname,argv[1]);
    printf("There are no parameters for task %s\n",taskname);
  } else {
    if (argc>1) strcpy(taskname,argv[1]);
    printf("  Task:   %s\n",taskname);
    for(i=0; i<n; i++)
      printf("  %-9s= %s\n",args[i].name,(args[i].value == NULL ?
                                          " " : args[i].value));
  }
}
/************************************************************************/
void doinput(argc,argv)
int argc;
char *argv[];
/*
------------------------------------------------------------------------*/
{
  int i,n;

  if(argc > 2) fprintf(stderr,"input: Extra arguments on line ignored.\n");

  if(argc==1) {
    fprintf(stderr,"input: no filename supplied\n");
    return;               /*  no filename supplied */
  }
  if(input_level+1 > MAXINPUT) {
        fprintf(stderr,"input: Too many nested inputs in %s\n",argv[1]);
        return;
  }
  fpinput[input_level] = fopen(argv[1],"r");
  if (fpinput[input_level] == NULL) {
    fprintf(stderr,"[input: file %s not found]\n",argv[1]);
    return;
  }
  input_level++;
}
/************************************************************************/
void doview(argc,argv)
int argc;
char *argv[];
/*
------------------------------------------------------------------------*/
{
  int i,n;
  FILE *fd;
  char name[MAXBUF],command[MAXBUF],*viewer,*testname;

  if(argc > 2) fprintf(stderr,"view: Extra arguments on line ignored.\n");
  testname = (argc>1 ? argv[1] : taskname );
  n = task_args(testname);
  if(n < 0){
    fprintf(stderr,"view: Found no documenation on task %s.\n",testname);
  } else if( n == 0){
    if (argc>1) strcpy(taskname,argv[1]);
    printf("There are no parameters for task %s\n",taskname);
  } else {
    if (argc>1) strcpy(taskname,argv[1]);
    sprintf(name,"%s.def",taskname);
    fd = fopen(name,"w");
    if(fd == NULL){
      fprintf(stderr,"Failed to open %s\n",name);
      return;
    }
    for(i=0; i<n; i++)
      fprintf(fd,"%-9s= %s\n",args[i].name,(args[i].value == NULL ?
                                          "" : args[i].value));
    fclose(fd);
#ifdef vms
    viewer = "edit";
#else
    {
      char *getenv();
      if((viewer = getenv("VISUAL")) == NULL)
        if((viewer = getenv("EDITOR")) == NULL) viewer = "vi";
    }
#endif
    sprintf(command,"%s %s",viewer,name);
    system(command);
    get_vars(name);
  }
}
#ifdef vms
/************************************************************************/
void dogo(bg,argc,argv)                            /* GO command in VMS */
int bg;                         /* bg still ignored in VMS */
int argc;
char *argv[];
/*  VMS version (no backgrounding/spawning yet)
------------------------------------------------------------------------*/
{
  int i,n,table;
  char line[MAXBUF],parameter[MAXBUF],*testname;
  struct {int length; char *pnt; } name,value;
#define LIB$K_CLI_GLOBAL_SYM 2
#define assign(descriptor,string) descriptor.length = strlen(string);\
                                  descriptor.pnt    = string

  if (bg) fprintf(stderr,"### Warning: SPAWNing ignored - no GOB\n");
  if(argc < 1) return;
  if(argc > 2) fprintf(stderr,"go: Extra arguments on line ignored.\n");
  testname = (argc > 1 ? argv[1] : taskname);
  n = task_args(argv[1]);
  if(n < 0){
    fprintf(stderr,"go: Found no documentation on task %s.\n",testname);
  } else {
    if (argc>1) strcpy(taskname,argv[1]);
    /* Check if the foreign command is defined. If not, define it. */
    assign(name,taskname);
    value.length = MAXBUF; value.pnt = line;
    if(lib$get_symbol(&name,&value) != 1){
      table = LIB$K_CLI_GLOBAL_SYM;
      sprintf(line,"$MIRBIN:%s.exe",taskname);
      assign(name,taskname); assign(value,line);
      lib$set_symbol(&name,&value,&table);
    }

/* Build up the command line. */

    strcpy(line,taskname);
    for(i=0; i<n; i++){           /* CHECK IF THIS STILL WORKS 15-jun-90 PJT */
      if(args[i].value != NULL){
        sprintf(parameter," %s=%s",args[i].name,args[i].value);
        strcat(line,parameter);
      }
    }
    system(line);
    printf("\n");   /* extra newline is apparently needed on VMS */
  }
}
#else
/************************************************************************/
void dogo(bg,argc,argv)                           /* GO command in UNIX */
int bg;
int argc;
char *argv[];
/*      Unix version
------------------------------------------------------------------------*/
{
  int i,n,pid;
  char *arg[MAXARGS],path[MAXBUF],line[MAXBUF],*s,**t,*testname;

  if(argc < 1) return;
  if(argc > 2) fprintf(stderr,"go: Extra arguments on line ignored.\n");
  testname = (argc > 1 ? argv[1] : taskname);
  n = task_args(testname);
  if(n < 0){
    fprintf(stderr,"go: Found no documentation on task %s.\n",testname);
  } else {
    if (argc>1) strcpy(taskname,argv[1]);
#if defined(PATHSEARCH)
    strcpy(path,taskname);       /* let shell search for exe file */
#else
    filename(path,"MIRBIN",taskname,"");     /* full name of exe file */
#endif
    s = line;		      /* start building up cmdline args in here */
    t = arg;		/* but also keep pointers to start of arguments */
    *t++ = taskname;
    for(i=0; i<n; i++){			/* loop foreach known keyword */
      if(args[i].value != NULL){
        if (strlen(s)+strlen(args[i].name)+strlen(args[i].value)>MAXBUF) {
            printf("Out of argument space to execute command\n");
            return;
        }
#if 1
             /* PJT - I think the next SNB section is not good - it is the
	     doset() command which should handle spaces properly ?? */
        sprintf(s,"%s=%s",args[i].name,args[i].value);
#else
             /* SNB */
        sprintf(s,"%s=%s ",args[i].name,args[i].value);
        strcat(s," ");
#endif
        *t++ = s;	           /* fix proper pointer into line[] */
        s += strlen(s) + 1;    /* point just beyond the last NULL in line[] */
      }
    }
    *t = NULL;					/* terminate arg[] properly */
    if (debug_level) {
        printf("Command:\n");
        for (t=arg; *t; t++)
            printf(" %s",*t);
        printf("\n");
    }

#ifdef unicos
    unicos_check(taskname);
#endif

/* Spawn off the command. */

    pid = fork();
    if(pid < 0) bug("go: Failed to fork a child process");
    if(pid == 0) {         /* inside child now */
      if (bg) fprintf(stderr,"[Job %s running in background]",taskname);
#if defined(PATHSEARCH)
      execvp(path,arg);         /* let shell look for exe */
#else
      execv(path,arg);          /* path contains full name of exe */
#endif
      printf("Errno = %d\n",errno);
      perror(" ");
      bug("go: Failed to exec the subprocess");
    } else {
        if (debug_level)
            fprintf(stderr,"[Parent jobid = %d - bg=%d]\n",pid,bg);
        if (!bg) wait(NULL);
    }
  }
}
#endif
#ifdef unicos
/************************************************************************/
void unicos_check(task)
char *task;
/*
  Check if the executable of the task is in the executable directory. If
  not, get it from the appropriate CFS directory.
------------------------------------------------------------------------*/
{
  char path[MAXBUF],cfs[MAXBUF],command[MAXBUF];

  filename(path,"MIRBIN",task,"");
  if(access(path,X_OK)){
    filename(cfs,"BINCFS",task,"");
    sprintf(command,"cfs get %s:%s",path,cfs);
    system(command);
  }
}
#endif
/************************************************************************/
void dohelp(argc,argv)
int argc;
char *argv[];
/*
------------------------------------------------------------------------*/
{
  char path[MAXBUF],command[MAXBUF],*pager, *mirpager, *hp;
  void showbin();

  if( (argc < 2 && taskname[0]==NULL) || (argc > 1 && *argv[1]=='?') ){
    showbin();
    return;
  }

  if(argc > 2) fprintf(stderr,"help: Extra arguments on line ignored.\n");

/* Check if we want help about Miriad itself. */

  if (argc < 2) {   /* take default taskname */
    hp = taskname;
  } else if(strcmp("miriad",argv[1])==0){
    doq();
    return;
  } else
    hp = argv[1];

/* Determine the pathname of the help file. */

  sprintf(path,"%s.doc",hp);
  if(access(path,R_OK)){
    filename(path,"MIRPDOC",hp,".doc");
    if(access(path,R_OK)) path[0] = 0;
  }
  if(!path[0]) {
    fprintf(stderr,"help: Cannot find and/or read help for %s\n",hp);
    fprintf(stderr,"Try: '?' or 'help ?' for other type of help\n");
  }else{
#ifdef vms
    pager = "type/page";
    sprintf(command, "%s %s",pager, path);
#else
    {
      char *getenv();
      if ( (pager = getenv("PAGER")) == NULL) pager = "more";
      if ( (mirpager = getenv("MIRPAGER")) == NULL) mirpager = "cat";
    }
    sprintf(command, "%s %s | %s", mirpager, path, pager);
#endif
    system(command);
  }
}

void showbin()
{
        char path[MAXBUF], command[MAXBUF], *pager;

        filename(path,"MIRBIN","","");
#ifdef vms
        pager = "dir";
#else
        pager = "ls";
#endif        
        sprintf(command,"%s %s",pager,path);
        printf("%s\n",command);
        system(command);
        return;
}

/************************************************************************/
void doq()
/*
------------------------------------------------------------------------*/
{
    printf("miriad is a very simple frontend to run miriad commands\n\n");
    printf("?                     this help\n");
    printf("version               display how miriad was compiled\n");
    printf("input <file>          process commands from input file\n");
    printf("set <key> <value>     set or show a keyword\n");
    printf("<key>=<value>         set a keyword\n");
    printf("unset <key(s)>        unset keyword(s)\n");
    printf("reset                 reset all keywords\n");
    printf("task {<task>]         set/show default taskname\n");
    printf("inp [<task>]          show settings of keywords to [task]\n");
    printf("go [<task>]           run [task]\n");
    printf("gob [<task>]          run [task] in background\n");
    printf("help [<task>]         help (on a task)\n");
    printf("view <task>           edit keyword file for task\n");
    printf("save/load [file]      save/load global keyword file\n");
    printf("tget/tput [<task>]    save/load [task] keyword file\n");
    printf("cd <dir>              change and/or show current directory\n");
    printf("exit,end              exit program, save variables\n");
    printf("quit                  quit program, do not save variables\n");
    printf("<task> [<par(s)>]     any other command is passed to the shell\n");
#if defined(TCL)
    printf("tcl                   Some more help on TCL\n");
#endif
    printf("\n");
}
/************************************************************************/
void dosave(argc,argv)
int argc;
char *argv[];
/*
------------------------------------------------------------------------*/
{
  char path[64];

  if(argc > 1)                          /* if an argument given */
    save_vars(argv[1]);                 /* assume it's the name to save */
  else {                                /* else */
    filename(path,"",taskname,".def");  /* use default task name */
    save_vars(path);
  }
}
/************************************************************************/
void doload(argc,argv)
int argc;
char *argv[];
/*
------------------------------------------------------------------------*/
{
  char path[64];

  if(argc > 1) {
    get_vars(argv[1]);
  } else {
    filename(path,"",taskname,".def");
    get_vars(path);
  }
}
/************************************************************************/
void dot(mode,argc,argv)
int mode,argc;
char *argv[];
/*
------------------------------------------------------------------------*/
{
    printf("tget/tput not yet implemented Lee...\n");
}
/************************************************************************/
void docd(argc,argv)
int argc;
char *argv[];
/*
------------------------------------------------------------------------*/
{
    if (argc == 2) {     /* change directory */
        if (chdir(argv[1]) != 0)
            printf("Failed to change directory %s\n",argv[1]);
        /* printf("Current directory is: ***\n"); */
    } else {            
        if (argc == 1) {    /* if one arg: show current dir */
#ifdef VMS
            system("sho def");
#else
            system("/usr/bin/pwd");
#endif
        } else
            printf("Usage: cd [<directory>]\n");
    }
}
/************************************************************************/
void dotcl(argc,argv)
int argc;
char *argv[];
/*
------------------------------------------------------------------------*/
{
#if !defined(TCL)
    printf("TCL is not compiled into this version of miriad\n");
#else
    printf("is OK\n");
#endif
}
/************************************************************************/
void dotask(argc,argv)
int argc;
char *argv[];
/*
------------------------------------------------------------------------*/
{
  if(argc > 1) {
    strcpy(taskname,argv[1]);
  } else {
    printf("Current default task is: %s\n",taskname);
  }
}
/************************************************************************/
void dosetenv(mode,argc,argv)
int mode, argc;
char *argv[];
/*
------------------------------------------------------------------------*/
{
    if (mode==1 && argc < 3) {
        fprintf(stderr,"Usage: setenv env_var value\n");
        return;
    } else if (mode==0 && argc < 2) {
        fprintf(stderr,"Usage: unsetenv env_var\n");
        return;
    }
    if (mode==1)
        newenv(argv[1],argv[2]);    /* setenv: */
    else if (mode==0)
        newenv(argv[1],"");         /* unsetenv: set to empty string */
    else
        bug("dosetenv: Unknown mode\n");    /* should never happen */
}
/************************************************************************/
void docommand(argc,argv)
int argc;
char *argv[];
/*
------------------------------------------------------------------------*/
{
  char *s, buffer[MAXBUF];
  int i, status;

  buffer[0] = '\0';
  for (i=0; i<argc; i++) {      /* accumulate all stuff into buffer */
    strcat(buffer,argv[i]);
    strcat(buffer," ");
  }
  s = buffer;
  while (*s == ' ' || *s == '\t')   /* skip whitespace */
    s++;
  if (*s && (*s=='#'))          /* treat lines starting with # as comment */
    return;     
#if defined(vms) || !defined(EXECL)

  system(buffer);         /* execute it by the host cmd.interpreter */

#else
                          /* In UNIX: pass it such that aliases are known */
  if (fork()==0) 
    if (shell==NULL)
        execl("/bin/sh","sh","-c",buffer,NULL);
    else
        execl(shell,shell,"-c",buffer,NULL);
#if defined(INTERRUPT)
  signal(SIGTERM, SIG_IGN);            /* ignore interrupts by the parent */
  signal(SIGQUIT, SIG_IGN);            /* till the baby dies */
  signal(SIGINT,  SIG_IGN);
#endif
  wait(&status);
#if defined(INTERRUPT)
  signal(SIGTERM, review);            /* restore status */
  signal(SIGQUIT, review);            /* of signals */
  signal(SIGINT,  review);
#endif

#endif /* vms/EXECL */
}
/************************************************************************/
void get_vars(name)
char *name;
{
  FILE *fd;
  char *argv[3],line[MAXBUF],*s, fname[MAXBUF];

  filename(fname,"MIRDEF",name,"");      /* make filename */
  fd = fopen(fname,"r");
  if(fd == NULL)return;
  argv[0] = "set";
  fprintf(stderr,"[All variables retrieved from %s]\n",fname);
  while(fgets(line,sizeof(line),fd) != NULL){
    s = line;
    while(*s == ' ' || *s == '\t')s++;      /* skip blanks */
    argv[1] = s;                            /* keyword */
    while(*s && *s != ' ' && *s != '\t' && *s != '=')s++;  /* skip key */
    while(*s && (*s == ' ' || *s == '\t' || *s == '='))*s++ = 0;  /* zero it */
    argv[2] = s;                        /* this points to the argument */
    while(*s && *s != '\n')s++;         /* patch newline char */
    *s = 0;
    if(*argv[2]) 
        doset(3,argv);
    else {
        dounset(2,argv);
    }
  }
  fclose(fd);
}
/************************************************************************/
void save_vars(name)
char *name;
{
  int i;
  VARIABLE *v;
  FILE *fd;
  char line[MAXBUF], fname[MAXBUF];


  Qkeys = 0;
  for(i=0; i<HASHSIZE; i++) /* walk through keywords to see if any are there */
    for(v = hashtable[i]; v != NULL; v = v->fwd)
        Qkeys++;
  if (Qkeys==0) {             /* if none found */
    fprintf(stderr,"[Variables not saved - no changes were made]\n");
    return;                 /* don't write them out */
  }

  /* First a check if writing the file is OK, we don't want to overwrite
   * a file that doesn't seem to look like a keyword file, surely */

  filename(fname,"MIRDEF",name,"");      /* make filename */
  fd=fopen(fname,"r");
  if (fd!=NULL) {               /* check that file !! */
    while(fgets(line,sizeof(line),fd) != NULL){
        if (strlen(line)>0) {
          if (strchr(line,'=')==NULL) {
            fprintf(stderr,"[Will not overwrite keywords to a file (%s)",fname);
            fprintf(stderr," which does not look like a keyword file\n");
            return;
          } else
            break;      /* OK, found an '=', assume its OK to overwrite */
        }
    }
    fclose(fd);   /* and close file again */
  } 
    
  fd = fopen(fname,"w");     /* open to write */
  if (fd==NULL) {               /* Perhaps allow a  y/n type-question */
    fprintf(stderr,"Warning: could not write file %s\n",fname);
    return;
  }
  for(i=0; i<HASHSIZE; i++)
    for(v = hashtable[i]; v != NULL; v = v->fwd)
      fprintf(fd,"%s=%s\n",v->name,v->value);
  fclose(fd);
  fprintf(stderr,"[All variables saved in %s]\n",fname);
}
/************************************************************************/
int task_args(task)
char *task;
/*
  This gets the arguments, and their values, for a particular task.
------------------------------------------------------------------------*/
{
  char line[MAXBUF],path[MAXBUF],keyword[MAXBUF],*s,*t;
  FILE *fd;
  VARIABLE *v;
  int n,hashval,found;

/* Check both the local directory, and the standard directory for the .doc
   file. */

  sprintf(path,"%s.doc",task);
  fd = fopen(path,"r");
  if(!fd){
    filename(path,"MIRPDOC",task,".doc");
    fd = fopen(path,"r");
  }
  if(!fd)return(-1);

/* Scan the .doc file for the keywords. */
  s = names;
  n = 0;
  while( fgets(line,sizeof(line),fd) != NULL){
    found = sscanf(line,"%%A %s",keyword) == 1;
    if(found){

/* Find the value of the parameter, stored in the hash table. */

      hashval = 0;
      t = keyword;
      while(*t) hashval += *t++;
      v = hashtable[hashval % HASHSIZE];
      while(v != NULL){
        if(!strcmp(keyword,v->name)) break;
        v = v->fwd;
      }
      args[n].value = (v == NULL ? NULL : v->value);

/* Save the name. */

      args[n].name = s;
      strcpy(s,keyword);
      s += strlen(keyword) + 1;
      n++;
    }
  }
  fclose(fd);
  return(n);
}
/************************************************************************/
void filename(out,envvar,name,type)
char *out,*envvar,*name,*type;
/*
  This makes a filename from the input components.
  In VMS 'envvar' is a logical, rest (Unix) treats it as an environment
  variable
------------------------------------------------------------------------*/
{
#if defined(vms)
  if(envvar && *envvar)sprintf(out,"%s:%s%s",envvar,name,type);
  else       sprintf(out,"%s%s",name,type);
#else
  char *s;
  if(envvar && *envvar){
    s = getenv(envvar);
    if(s == NULL){
      sprintf(out,"### Unable to find environment variable %s.",envvar);
      bug(out);
    } else {
      sprintf(out,"%s/%s%s",s,name,type);
    }
  }else 
    sprintf(out,"%s%s",name,type);
#endif
}
/************************************************************************/
void bug(message)
char *message;
/*
  This prints an error message, then exits.
------------------------------------------------------------------------*/
{
  fprintf(stderr,"### %s\n",message);
  exit(1);
}
/************************************************************************/
set_defaults(ac,av)
int ac;
char *av[];
{
/*
    This gets any command line parameters
    -b dir      take this dir as MIRBIN (envp)
    -p dir      take this dir as MIRPDOC
    -d dir      take this dir as MIRDEF
    -s lastexit default name for lastexit file
    -g          debugging on
    -?/h        inline help
------------------------------------------------------------------------*/
    int i;
    char *cp;

    i=1;
    while(i < ac) {         /* avoid using getopt here */
        cp = av[i];
        if (*cp++=='-') {
            switch(*cp) {
              case 'g':
                    debug_level = 1;
                    printf("Debug turned on\n");
                    break;
              case 'b':
                    if (++i >= ac) bug("missing -b name, try -help");
                    newenv("MIRBIN",av[i]);
                    bindir = av[i];
                    break;
              case 'p':
                    if (++i >= ac) bug("missing -p name, try -help");
                    newenv("MIRPDOC",av[i]);
                    docdir = av[i];
                    break;
              case 'd':
                    if (++i >= ac) bug("missing -d name, try -help");
                    newenv("MIRDEF",av[i]);
                    defdir = av[i];
                    break;
              case 's':
                    if (++i >= ac) bug("missing -s name, try -help");
                    lastexit = av[i];
                    break;
              case 'h':
              case '?':
                    usage(av[0]);                  /* will also exit(0) */
              default:
                    fprintf(stderr,"%s: Illegal option %s\n",av[0],cp);
            }
            i++;
        } else
            fprintf(stderr,"%s: Illegal option %s\n",av[i++],--cp);
    } /* while */
    shell = getenv("SHELL");
}

usage(name)
char *name;
{
    fprintf(stderr,
      "Usage: %s [-b bindir] [-p docdir] [-d defdir] [-g] [-s lastexit]\n",
      name);
    fprintf(stderr,"Options:\n");
    fprintf(stderr,"  -?/h        This help\n");
    fprintf(stderr,"  -g          Turn debugging on (verbose)\n");
    fprintf(stderr,"  -b bindir   Use this to search for commands\n");
    fprintf(stderr,"  -p docdir   Use this to search for documentation\n");
    fprintf(stderr,"  -d defdir   Use this in I/O lastexit/.def files\n");
    fprintf(stderr,"  -s lastexit Name of default lastexit file \n");
    exit(0);
}


void doversion()
{
    printf("Miriad Version=%s\n",VERSION_ID);
    printf("Compile directives used were: \n");
    printf("PATHSEARCH: ");
#if defined(PATHSEARCH)
    printf(" on\n");
#else
    printf(" off\n");
#endif
    printf("INTERRUPT:  ");
#if defined(INTERRUPT)
    printf(" on\n");
#else
    printf(" off\n");
#endif
    printf("TCL:        ");
#if defined(TCL)
    printf(" on\n");
#else
    printf(" off\n");
#endif
    printf("READLINE:   ");
#if defined(READLINE)
    printf(" on\n");
#else
    printf(" off\n");
#endif

    if (bindir)
        printf("BINDIR = %s\n",bindir);
    if (docdir)
        printf("DOCDIR = %s\n",docdir);
    if (defdir)
        printf("DEFDIR = %s\n",defdir);
}

/************************************************************************/
newenv(var,value)
char *var, *value;
{
#if defined(NEWENV)
    char **ep, **newep, **epfrom, **epto, *cp;
    int  elen, vlen, nev, i;
/*
        enter a new environment variable. If 'value' is NULL, erase it.
        Some OS's have putenv(), the opposite to getenv(), but they 
        typically do not work with the third parameter to main().

    input:  var         name of env. variable
            value       value of env. variable
------------------------------------------------------------------------*/
    vlen = strlen(var);

    for (nev=0, ep=menviron; *ep != NULL; ep++) {            /* loop all vars */
        nev++;                                               /* count # vars */
        elen = strlen(*ep);                       /* length of this variable */
        cp = strchr(*ep,'=');                        /* look for an '=' sign */
        if ( (int)(cp - *ep) == vlen && strncmp(*ep,var,vlen)==0) { /* found */
            if (vlen+strlen(value)+2 > elen) {                 /* reallocate */
                cp = malloc(vlen+strlen(value)+2);
                if (cp==NULL) {
                    fprintf(stderr,"No space to expand variable %s\n",var);
                    return;
                }
                *ep = cp;
            }
            sprintf(*ep,"%s=%s",var,value);
            return;                                  /* and return to caller */
        }
    }
    /* if it got here: add it as a new environment variable */

    if (debug_level)
        fprintf(stderr,"### (%d) Adding %s=%s to environment\n",nev,var,value);

    newep = (char **)  malloc( (nev+2) * sizeof(char **) );   /* allocate new */
    for (i=0, epfrom=menviron, epto=newep; i<nev; i++)     /* copy old stuff */
       *epto++ = *epfrom++;   
    cp = malloc(vlen+strlen(value)+2);          /* allocate for new one */
    if (cp==NULL) {
        fprintf(stderr,"No memory to add environment %s=%s\n",var,value);
        return;
    }
    sprintf(cp,"%s=%s",var,value);
    *epto++ = cp;       /* put new thing in array */
    *epto = NULL;       /* terminate array */
    menviron = newep;    /* set new environment */
#else
    printf("### newenv disabled; cannot set %s=%s\n",var,value);
#endif
}

/* getenv: a local version, in case *getenv() is not supplied by your OS */

#if defined(GETENV)
char *getenv(var)
char *var;
{
    char **ep, *cp;
    int vlen, elen;

    vlen = strlen(var);
    for(ep=menviron; *ep != NULL; ep++) {
        elen = strlen(*ep);                       /* length of this variable */
        cp = strchr(*ep,'=');                       /* look for an '=' sign */
        if ( (int)(cp - *ep) == vlen && strncmp(*ep,var,vlen)==0) { /* found */
            fprintf(stderr,"Found %s\n",*ep);
            cp++;
            return(cp);
        }
    }
    return(NULL);
}

#endif

#if defined(INTERRUPT)
void review()
{
  fprintf(stderr,
       "Miriad shell cannot be interrupted, type '?' for info\n");
}
#endif


/* general miriad command dispatcher */
/* Can also be used in TCL environment (see later for TCL) */

do_miriad(argc,argv)
int argc;
char **argv;
{
    int more=1;

    if(!argc);                      /* do nothing if no args were given */
    else if(!strcmp(argv[0],"?"))        {doq();}
    else if(!strcmp(argv[0],"version"))  {doversion();}    
    else if(!strcmp(argv[0],"task"))     {dotask(argc,argv); }
    else if(!strcmp(argv[0],"set"))      {doset(argc,argv); Qkeys++; }
    else if(!strcmp(argv[0],"unset"))    {dounset(argc,argv); Qkeys++; }
    else if(!strcmp(argv[0],"setenv"))   {dosetenv(1,argc,argv); }
    else if(!strcmp(argv[0],"unsetenv")) {dosetenv(0,argc,argv); }
    else if(!strcmp(argv[0],"inp"))      {doinp(argc,argv); }
    else if(!strcmp(argv[0],"go"))       {dogo(0,argc,argv); }
    else if(!strcmp(argv[0],"gob"))      {dogo(1,argc,argv); }
    else if(!strcmp(argv[0],"help"))     {dohelp(argc,argv); }
    else if(!strcmp(argv[0],"view"))     {doview(argc,argv); Qkeys++; }
    else if(!strcmp(argv[0],"save"))     {dosave(argc,argv); }
    else if(!strcmp(argv[0],"load"))     {doload(argc,argv); }
    else if(!strcmp(argv[0],"tget"))     {dot(0,argc,argv); }
    else if(!strcmp(argv[0],"tput"))     {dot(1,argc,argv); }
    else if(!strcmp(argv[0],"input"))    {doinput(argc,argv);}
    else if(!strcmp(argv[0],"reset"))    {doreset(argc,argv);}
    else if(!strcmp(argv[0],"cd"))       {docd(argc,argv);}
    else if(!strcmp(argv[0],"exit"))     {more = 0; }
    else if(!strcmp(argv[0],"quit"))     {Qkeys = 0; more = 0; }
    else if(!strcmp(argv[0],"end"))      {more = 0; }
    else if(!strcmp(argv[0],"tcl"))      {dotcl(argc,argv);}
    else                                 {docommand(argc,argv);}

    if (more==0) end_miriad();
}


/* The following pages are dedicated to some GNU READLINE experiments - 
   they are typically not used in your active version of miriad.c */
#if defined(READLINE)
ini_readline()
{
    printf("*** Experimental GNU READLINE installed ***\n");
}
#endif /* READLINE */


/* The following pages are dedicated to some TCL experiments - they are
   typically not used in your active version of miriad.c */

#if defined(TCL)

cmdEcho(clientData, interp, argc, argv)
    char *clientData;
    Tcl_Interp *interp;
    int argc;
    char **argv;
{
    int i;

    for (i = 1; ; i++) {
        if (argv[i] == NULL) {
            if (i != argc) {
                echoError:
                sprintf(interp->result,
                    "argument list wasn't properly NULL-terminated in \"%s\" command",
                    argv[0]);
            }
            break;
        }
        if (i >= argc) {
            goto echoError;
        }
        fputs(argv[i], stdout);
        if (i < (argc-1)) {
            printf(" ");
        }
    }
    printf("\n");
    return TCL_OK;
}


dP(clientData)
    char *clientData;
{
    printf("Deleting command with clientData \"%s\".\n", clientData);
}

cmdTclHelp(clientData, interp, argc, argv)
    char *clientData;
    Tcl_Interp *interp;
    int argc;
    char **argv;
{
    printf("Some more help on TCL will be available here\n");
    printf("Your friendly programmers was just a bit lazy here late at night\n");
}
cmdCd(clientData, interp, argc, argv)
    char *clientData;
    Tcl_Interp *interp;
    int argc;
    char **argv;
{
    int i;

    if (argc != 2) {
        sprintf(interp->result, "wrong # args:  should be \"%.50s count\"",
                argv[0]);
        return TCL_ERROR;
    }
    i = chdir(argv[1]);
    if (i==0)
      return TCL_OK;
    else {
      sprintf(interp->result, "%s: No such directory",argv[1]);
      return TCL_ERROR;
    }
}

int
cmdTest(clientData, interp, argc, argv)
    char *clientData;
    Tcl_Interp *interp;
    int argc;
    char **argv;
{
    printf("CmdTest: %s\n",argv[0]);
    printf("  Command is loaded as a special TCL command\n");
    printf("  but is not implemented yet; us miriad INP/GO command\n");

    return TCL_OK;
}

int
cmdMiriad(clientData, interp, argc, argv)
    char *clientData;
    Tcl_Interp *interp;
    int argc;
    char **argv;
{
    int do_miriad();

    printf("CmdMiriad: %s\n",argv[0]);
    do_miriad(argc,argv);
    return TCL_OK;
}

install_tcl()
{
  int result;

  printf("[Installing TCL as (additional) command language]\n");

  tcl_interp = Tcl_CreateInterp();
  tcl_buffer = Tcl_CreateCmdBuf();

  Tcl_CreateCommand(tcl_interp,"echo",cmdEcho,(ClientData)"echo",dP);
  Tcl_CreateCommand(tcl_interp,"cd",cmdCd,(ClientData)"cd",  dP);
  Tcl_CreateCommand(tcl_interp,"tcl",cmdTclHelp,(ClientData)"tcl",  dP);

  Tcl_CreateCommand(tcl_interp,"?",cmdMiriad,(ClientData)"?",dP);
  Tcl_CreateCommand(tcl_interp,"version",cmdMiriad,(ClientData)"version",dP);
  Tcl_CreateCommand(tcl_interp,"task",cmdMiriad,(ClientData)"task",dP);
/*  Tcl_CreateCommand(tcl_interp,"set",cmdMiriad,(ClientData)"set",dP);  */
  Tcl_CreateCommand(tcl_interp,"unset",cmdMiriad,(ClientData)"unset",dP);
  Tcl_CreateCommand(tcl_interp,"unsetenv",cmdMiriad,(ClientData)"unsetenv",dP);
  Tcl_CreateCommand(tcl_interp,"inp",cmdMiriad,(ClientData)"inp",dP);
  Tcl_CreateCommand(tcl_interp,"go",cmdMiriad,(ClientData)"go",dP);
  Tcl_CreateCommand(tcl_interp,"help",cmdMiriad,(ClientData)"help",dP);
  Tcl_CreateCommand(tcl_interp,"exit",cmdMiriad,(ClientData)"exit",dP);
  Tcl_CreateCommand(tcl_interp,"quit",cmdMiriad,(ClientData)"quit",dP);
  Tcl_CreateCommand(tcl_interp,"input",cmdMiriad,(ClientData)"input",dP);
  Tcl_CreateCommand(tcl_interp,"save",cmdMiriad,(ClientData)"save",dP);
  Tcl_CreateCommand(tcl_interp,"load",cmdMiriad,(ClientData)"load",dP);
  Tcl_CreateCommand(tcl_interp,"view",cmdMiriad,(ClientData)"view",dP);

  load_environ();
  if (bindir)
    load_commands(bindir);
  
  result = Tcl_Eval(tcl_interp,"source ~/.tclrc",0,0);
  if (*tcl_interp->result != 0)
    printf("TCL_STARTUP: %s\n",tcl_interp->result);
  getline(NULL,"echo Testing cmdnow-mode in getline in TCL mode");
}

load_environ()
{
    char *cp, line[256], *strchr();
    int i, k;

    k = 0;
    for (i=0, cp=menviron[0]; cp != NULL; cp=menviron[++i]) {
        strcpy(line,cp);
        cp = strchr(line,'=');
        if (cp==NULL) {
            printf("### %s: not a proper environment variable\n",line);
            continue;
        }
        *cp = '\0';  /* patch it */
        cp++;
        Tcl_SetVar(tcl_interp,line,cp,1);
        k++;
    }
    printf("Loaded %d environment variables.\n",k);
}

/* 
 *  load all commands from a directory 'dir' as new tcl commands.
 *  this gives some possibility to retrieve command information
 *  back from command. Of course they have to know about it..
 */

load_commands(dir)
char *dir;
{
    FILE *fp;       /* pipe */
    char line[64];
    int k;

    sprintf(line,"ls %s",dir);
    fp = popen(line,"r");
    if (fp==NULL) {
        printf("Could not initialize directory %s\n",dir);
        return(0);
    }
    k=0;
    while (fgets(line,64,fp) != NULL) {
        line[strlen(line)-1] = '\0';
        Tcl_CreateCommand(tcl_interp, line, cmdTest, (ClientData) line,
            dP);
        k++;
    }
    printf("Loaded %d commands from %s.\n",k,dir);
}

#endif /* TCL */

