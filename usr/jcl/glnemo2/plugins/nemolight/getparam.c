/* 
 * GETPARAM.C: command-line processing functions for NEMO and ZENO
 *             this also encompassed most of NEMO's user interface
 *
 *             This package has a long and winded history... 
 *
 * Nov. 1986  Joshua Edward Barnes, Princeton NJ. - original version
 * Oct. 1987  Peter Teuben: added the 'host=' system keyword 
 * Mar. 1988  PJT: added 'debug=' and 'yapp =' system keyword
 *              and used data-history mechanism (HeadlineTag, now HistoryTag)
 * Apr. 1988  PJT: added interactive mode using help=1
 * Jun. 1988  JEB + PJT : --merged IAS and UIUC versions -- 
 * Dec. 1988  PJT: help_levels extended : pick your favorite EDITOR mode 
 * Jan. 1989  PJT: cleared up help-levels and some debug output
 * Mar. 1989  PJT: general calls to nemoinp to parse
 * Jun. 1989  PJT: added setparam - added help=a option for help=8 nemotool
 *                 'key=val' is allowed to have some whitespace, i.e.
 *                 ' key = val' This is for NemoTool
 * Jul. 1989  PJT: interrupting to REVIEW section - help=h installed 
 *                 keyword file is now <prognam>.def (Mirth compatible)
 * Nov. 1989  PJT: date_id() added, fixed annoying save_history() feature
 * Jan  1990  PJT: turning off date_id() [bug??? on sub3]; added USAGELOG
 * Feb  1990  PJT: read include files using key=@file  (no nesting)
 * 14-mar-90  PJT: V2.3 made GCC happy and fixed bug introduces by GCC!
 * 19-mar-90  PJT: V2.3a help=t fixed up a but, helpvec formatting too
 *  9-apr-90  PJT: V2.4 review can be used to chain programs... (NONO not yet)
 * 21-may-90  PJT: V2.5 yapp_string  for yapp_sm now installed
 * 17-jul-90  PJT: V2.5a isolated debug_level through function calls
 *  7-sep-90  PJT: V2.6 - REDIR & system keywords stdout, stderr added
 *                            Won't work the way I did it - disabled
 *                      - help=t also shows the default now
 *  1-oct-90  PJT: V2.6a moved MACROREAD into getparam
 *  7-oct-90           b more ANSI - index -> strchr
 * 11-nov-90           c repaired bug with @include files
 * 21-nov-90           d added 'isaparam'
 * 16-jan-91           e added support for MNEM_ malloc() debug
 * 22-mar-91  PJT V2.7 ERROR environment variable
 * 25-mar-91           a include ieee_handler when debug>0
 * 25-may-91           b appease some gcc warning messages - some Cray tests
 *                       also got rid of MICRO and REDIR stuff
 * 27-jun-91  PJT      c Fixed bug when empty string given to getXparam()
 * 11-jul-91           d fixed bug when empty line in .def file
 *                       blank lines will still crap out
 * 12-oct-91  PJT V2.8 Various small compiler complaints for BorlandC compiler
 *  7-nov-91       2.9 Added hasvalue() function
 * 22-nov-91      2.9a AIX can't handle if(getenv()) if non-existing e.v.
 *                     so use 'ev' as temporary pointer
 *  6-feb-92  PJT 2.9b output KHOROS pane file when help=z is used
 * 13-feb-92         c usage string request via help= or so...; help=u
 * 23-feb-92         d implemeted dmm(3) menu interface via help=8 HELP_DMM
 * 25-feb-92	     e happy gcc2.0
 * 12-mar-92	     f added local nemobin: no usage if no $NEMOBIN
 * 27-apr-92         g parname: space is not a key=val separator
 *  6-jul-92         h help=t more verbose output
 * 18-aug-92         i protected matchname() agains NULL entry
 *                     removed some lint complaints
 *                   j special handling of VERSION keyword to be able to
 *                     warn of old VERSION def file is read.
 * 31-oct-92         k fixed bug when VERSION= absent
 * 18-mar-93	     l default USAGE is now off
 *  6-apr-93         m error_level and debug_level back in their .c routines
 *		       stop() moved to error.c -- still not clean !!!
 * 21-apr-93	     n fixed bug in help=2 reading due to VERSION . ...
 *  6-jun-93         o fixed bugs with badly contructed defv[]'s - HatCreek
 * 15-jun-93         p fixed bugs introduced by previous 'fix' - HatCreek
 * 14-nov-93	V2.10 consolidated two separate mod tracks ...		PJT
 * 18-jan-94    V2.11 solaris (5.3) support				pjt
 *                    -> no more host= remote execution
 *                       ieee disabled (problems finding library)
 * 21-feb-94         a - ansi *** unfinished
 * 25-nov-94    V3.0  big overhaul to use new 'struct keyword'    PJT
 *                    From 2000 to 1600 lines (excluding TCL & READLINE)
 *                    + minimim match keywords
 *                    + special system key processing, added new argv=
 *                    + tcl support
 *                    - no support for DMM, UsageLog, host= (rsh)
 * 12-feb-95          fixed keyword initializer bug
 *                    added reference count (key[].upd)
 *  7-mar-95          no more need for HELPVEC, defined MAXKEYLEN
 *  5-jul-95          show NEMO version in help, fixed help=i bug
 * 11-sep-95    3.0c  allow _ in env.var. in eval_keys
 * 14-aug-96    3.0d  fixed counting bug and introduced maxkeys for this
 * 13-sep-96    3.0e  fixed bug with blank environment variables (SunOS5.5)
 *  6-mar-97    3.0f  added SINGLEPREC output for help=?
 * 18-jan-98    3.0g  added version check for debug>0
 * 14-jan-99    3.0h  support single options -help, --help
 *  3-feb-99    3.1   support for ZENO defv[] vectors with ";.." comments
 * 21-feb-99    3.1a  added -help and -version
 *                    (still no -help and -clue yet)
 *  8-oct-99    3.1b  all gets() replaced with fgets()
 *  6-jan-00       c  fixed bug in fgets() arguments
 * 18-jan-00       d  added cntparam; added [defaults] to help=h; 
 *                    also use 'help' and '-h' to give help.
 *                 e  gcc warnings
 *  8-apr-01    3.2   added getfparam()
 * 21-may-01       a  returned 'int history' to history.c for 'uNEMO'
 * 31-may-01       b  added jJ (dutch/german) to boolean 'true'
 *  7-jun-01       c  using NEMO_VERSION instead of VERSION
 * 23-jun-01       d  PRECISION now in ZENO mode
 *  3-jul-01       e  keybuf always present, not just in INTERACT mode
 * 12-sep-01       f  nemo_file_size now
 * 15-sep-01    3.3   indexed keywords, but oops, can't do MINMATCH
 *                    and does not do keyword file I/O yet, no macroread 
 * 16-sep-01       a  fixed macroread
 * 19-oct-01       b  fixed recursion for those not using nemo_main
 * 24-oct-01       c  some parsing error warning's are now fatal errors
 * 13-nov-01       d  demoted some warnings to dprintf's
 * 25-nov-01       e  fixed hasvalue() for indexed keywords
 * 10-dec-01       f  indexing now based at 0, not 1 (for wasp)
 * 16-jan-02       g  wrote findakey to fix various indexing shortcomings
 * 17-jan-02       h  fix minor indexing, fix khoros output
 * 18-jan-02       i  fix reading indexed from keyword (.def) files
 * 20-jan-02       j  putparam can now create indexed keywords on the fly
 *  5-feb-02       k  findaparam return NULL if illegal indexed
 * 24-apr-02    3.4   added placeholder for new outparam() stuff
 *  1-aug-02       a  try and free() some allocated local data
 *  2-aug-02       b  write kiddy version of outkeys= processing
 *  4-dec-02       c  use MAXBUF for keys - btw, wc is not 2800+ lines :-(
 * 21-may-03       d  document the system keywords in help=\?
 * 12-jul-03       e  allow -- to stop parsing
 * 15-oct-03       f  fix for version.h to give multiple static's in g++ 
 * 23-jan-04       g  fix prototypes for intel compiler
 * 13-may-04       h  add cpu/real time reporting for help=c 
 * 23-jun-04       i  fix bug free's version_  (JC Lambert)
 * 20-sep-04       j  disabled the version checking by default
 * 18-oct-04       k  support for the CVS ID (cvsid)
 * 28-dec-04       l  help/report CVS ID properly (help=I)
 * 30-jun-05       o  parse integers with 0x as hex numbers
 * 14-dec-06    3.5   changed some printf()'s to fprintf(stderr,....)'s
 * 13-sep-07 WD    a  print name of name whose key is not found
 * 23-oct-07    3.6   fix bug when ZENO style defv[] is used, more ZENO compat
 * 22-aug-08       a  add help=m to show memory usage at end
 * 26-sep-08 WD    b  add MPI proc info to some stderr output
 * 27-Sep-10 JCL   c   MINGW32/WINDOWS support

  TODO:
      - what if there is no VERSION=
      - process system keywords in readkeys()
      - allow a (file) keyword to not signify a pipe
        (-) but a popen() based process (POPENIO)

      - instead of using local_error() to catch recursion or
        errors when nemo_main() is not used, perhaps for nkeys==0
	should do a magic init with a dummy defv[] array
        such that getargv0() does not recurse....

  BUGS/FEATURES

      - no filelocking for keyword file I/O
      - gets() fails after somebody has used the scanf function
  	(see promptparam)
  	        allows to process it - and of course dies
      - ESC go key in help=2 prompting mode does not work
      - @macro and $key references get expanded as strings.

      - handle long arguments better. It seems the shell can only
        do 256, but macros (key=@macro instead of key=value) can
	do much more (up to our liking really). getparam should
	warn for this and  suggest a fix for the user

	C-shell variables and command line argument handling all
	have different maxima, however 2000 seems to be ok. Much
	higher, C-shell seems to die around  4000, argv[] around
	10000 on linux-2002.  On solaris they both seem to be 
	around 1000
 
  SOME ALTERNATIVE RESOURCES
	argtable  http://argtable.sourceforge.net/doc/html/index.html
	clig	  http://wsd.iitb.fhg.de/~kir/clighome/
	perl      http://www-106.ibm.com/developerworks/linux/library/l-perl-speak.html
	parseargs
	opt
        getopt
	gengetopt http://www.gnu.org/software/gengetopt/gengetopt.html
	opag      http://www.zero-based.org/software/opag/
 */

#define GETPARAM_VERSION_ID  "3.6a 22-aug-08 PJT"

/*************** BEGIN CONFIGURATION TABLE *********************/

#if !defined(unicos) && !defined(__BORLANDC__) && !defined(mc68k)
# define INTERACT        /* want interactive input ?              */
# define NEMOINP         /* want nemoinp parsing in getXparam ?   */
#endif

#define PUTPARAM        /* allow putparam ?                       */
#define MACROREAD       /* ability to use @include files          */
#define POPENIO         /* ability to use <popenr or >popenw      */
#define HISTORY         /* have history.c implemented ?           */
#define MINMATCH	/* allow keywords to be minimum matched?  */
#define INDEXED         /* allow keywords to be dynamically indexed */
#define OUTKEYS         /* allows writing binary export keyword   */
#define DASHDASH        /* allows stop parsing when -- seen       */

#if 0
#define TCL7		/* TCL (old V7) support?		  */
#define READLINE	/* READLINE support?			  */
#endif

#if 0
# define INTERRUPT       /* reset when you don't want REVIEW interrupt */
#endif


/*************** END CONFIGURATION TABLE ***********************/

#include <stdinc.h>
#include <getparam.h>
#include <extstring.h>
#include <filefn.h>
#include <strlib.h>
#include <history.h>

#ifndef __MINGW32__
#include <sys/types.h>
#include <sys/times.h>
#else
# include <sys/time.h>
#endif
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <unistd.h>
#include <ctype.h>


#if defined(TCL7)
# include <tcl.h>
#endif

#if defined(READLINE)
  extern char *readline(string);
#endif

#if defined(TCL)
  local void myreadline(string,int);
  local void ini_readline(void);
  local void end_readline(void);
  local void stripwhite(char *);
#endif

#if defined(sun)
#   include <floatingpoint.h>               /* special IEEE handler for Sun's */
#endif

#if defined(MIRIAD)
#define KEYVAL 1			/* this would enforce usage of key= */
#endif

#if defined(INTERRUPT)
# include <signal.h>
#endif

#if defined(INTERACT)
#ifndef __MINGW32__
# include <sys/ioctl.h>           /* this is a BSD thing - won't work on SYSV */
#endif
                                  /* fcntl.h doesn't have the right stuff */
# if defined(SOLARIS)
#  include <sys/termio.h>	/* for Solaris - can it hurt on SUN OS 4.x ?  */
# endif

#endif /* INTERACT */

#ifndef MAXBUF
# define MAXBUF 1024
#endif

#define MAXKEYLEN 16


/* keyword is a helper struct containing a keyword. From this the code will
 * build a linear list of all program keywords, and those that are indexed,
 * will follow a linked list. E.g. for 4 program keywords 'a', 'b', 'c' and 'd',
 * of which 'c' is indexed and 3 had been supplied on the command line, the
 * keys[] array will have the following structure:
 *
 *   keys[0]    argv0
 *   keys[1]    a
 *   keys[2]    b
 *   keys[3]    c# ->  c4 -> c0 -> c1
 *   keys[4]    d
 */

typedef struct keyword {      /* holds old all keyword information           */
  string keyval;              /* pointer to original R/O (defv) data         */
  string key;                 /* keyword/name                                */
  string val;                 /* value                                       */
  string help;                /* help/comment -- may be empty                */
  int count;                  /* update count; 0=original default            */
  int upd;                    /* 0=read 1=unread original 2=unread updated   */
  int flag;                   /* zeno flags (unused)                         */
  int indexed;                /* index+1 of keyword, if so; 0 if none        */
  struct keyword *next;       /* pointer to next indexed keyword             */
} keyword;

typedef struct keyword_out {  /* a simple keyword, only meant for outkeys=   */
  string key;                 /* keyword name                                */
  char type;                  /* data type: 'c', 'd', 'f', 'i', 'l', 'b'     */
  int ndat;                   /* number of data elements                     */
  void *dat;                  /* pointer to data (malloc/free)-able          */
} keyword_out;

/*      *       *       *       variables        *       *       */

/* local variables - not visible to outside world */

   local char key_filename[MAXBUF];         /* filename to store keywords */
   local char keybuf[MAXBUF];               /* buffer for keyword i/o to file */
   local string version_e=NULL;  /* extern version (keyfile or commandline)    */
   local string version_i=NULL;  /* internal version (executable)              */

#if defined(INTERRUPT)
   local int initparam_done=0;                 /* trap for too early review() */
#endif

   local string progname=NULL;          /* is really pointer to keys[0].val   */
   local int nkeys, maxkeys;
   local keyword *keys = NULL;          /* point to array of program keywords */
   local int getparam_argc = 0;		/* count commmand line args           */
   local int dd_argc = 0;      /* number of arguments after and including  -- */
   local string *dd_argv = NULL;        /* starts at "--" if present          */

/* global variables - must be visible to the outside world */
/*      some defined here, others defined elsewhere        */

extern int debug_level; /* see dprintf.c   from DEBUG env.var. */
extern bool mpi_proc;   /* dprintf.c */
extern int  mpi_rank;   /* dprintf.c */
extern int error_level; /* see error.c     from ERROR env.var. */
extern string usage;    /* see program.c or usage.c */
extern string cvsid;    /* see program.c or cvsid.c */
extern string defv[];   /* see program.c or defv.c */
extern string outdefv[];/* see program.c or outdefv.c */
extern char **environ;  /* environment variables */

extern int nemo_history;     /* 0=no history written; see history.c */

int yapp_dev = 0;       /* interface hidden keyword yapp to plotting pkg */
int help_level = 0;     /* hidden keyword help for interactive prompting */
int bell_level = 0;     /* noisy terminal when prompted */
int review_flag = 0;    /* review keywords and optional help=8 chaining ? */
int tcl_flag = 0;       /* go into TCL when all parameters set before go */
int report_cpu = 0;     /* report time and cpu usage; activated with help=T */
int report_mem = 0;     /* report memory usage (machine specific) */
string yapp_string = NULL;  /* once only ? */
string help_string = NULL;  /* cumulative ? */
string argv_string = NULL;  /* cumulative ? */
string error_string = NULL;

local clock_t clock1, clock2;
#ifndef __MINGW32__
 local struct tms tms1, tms2;
#else
 local struct tm tms1, tms2;
#endif

#if defined(TCL7)
  Tcl_Interp *tcl_interp;
  int cmdInp(), cmdCd(), cmdLoad(), cmdTest(), cmdGo(), cmdHelp(), cmdEcho();
  int cmdKeys();
  int dP();
#endif

/* local functions */

local void initparam_out(void);
local void finiparam_out(void);
local void writparam_out(string key, char type, int ndat, void *dat);
local string nemokeys(void);

local void scan_environment(void);
local void save_history(string *argv);
local void printhelp(string help);
local void newline(int n);
local void showconfig(void);
local void showsystem(void);
local void printusage(string *defv);
local string get_macro(char *mname);
local void eval_keys(void);
local void setparam(char *par, char *val, char *prompt);
local int readparam(char *buffer, char *text);
local void beep(void);
local string parname(string arg);
local string parvalue(string arg);
local string parhelp(string arg);
local int findkey(string name);     // to be obsoleted
local keyword *findakey(string name);
local int set_indexed(string name, int *idx);
local int addindexed(int j, string keyval, int idx);
local void writekeys(char *mesg);
local void readkeys(string mesg, bool first);
local void set_debug(string);
local void set_error(string);
local void set_help(string);
local void set_review(string);
local void set_tcl(string);
local void set_yapp(string);
local void set_outkeys(string);
local void local_error(string);
local void local_exit(int);
local void report(char);

/* external NEMO functions */

extern int nemo_file_size(string);
extern string date_id(void);
extern double cputime(void);

/* Mask values for help levels */

#define HELP_DEFIO  0x01
#define HELP_PROMPT 0x02
#define HELP_EDIT   0x04
#define HELP_GLOBAL 0x08		/* Note: not used yet */
#define HELP_FUTURE 0x10		/* Note: not used yet */


/*
 * INITPARAM: initalize parameter lists and handle special help and
 *            review requests. This should be the first routine any
 *            NEMO program (main()) calls.
 *            See also $NEMO/src/kernel/misc/nemomain.c[c]
 */

void initparam(string argv[], string defv[])
{
    int i, j, nzeno, idx;
    bool posflag, useflag, defflag;
    string name;
    keyword *kw;

    progname = argv[0];        /* cheat in case local_error called early on */
    scan_environment();              /* scan env.var and set what is needed */
#if defined(sun) && defined(USE_IEEE)
    if(ieee_handler("set","common",ieee_nemo))    /* need to trap FP errors */
        warning("Error setting Sun ieee handler");
#endif
    
    nkeys = xstrlen(defv, sizeof(string));     /* count # params + progname */
    keys = (keyword *) allocate(nkeys * sizeof(keyword));
    maxkeys = nkeys - 1;

    keys[0].keyval = argv[0];
    keys[0].key = scopy("argv0");
    keys[0].val = progname = tail(argv[0]);
    keys[0].help = scopy("Program name");
    keys[0].count = 0;
    keys[0].upd = 0;
    nzeno =  (*defv[0] == ';');

    /* TODO:  if in zeno mode, it needs to skip EACH ; string, since it's a comment/help */

    for (i = 1, j=0; i < nkeys; i++, j++) {		/* loop over key defs */
        if (defv[j] == NULL) break;
	if (nzeno) {
	  if (j==0) {
	    usage = defv[0];
	    dprintf(1,"ZENO Usage: %s\n",usage);
	  } else
	    dprintf(1,"ZENO: %s\n",defv[j]);
	  j++;
	  if (defv[j] == NULL) {
	    nkeys = i;
	    dprintf(2,"ZENO: nkeys=%d %s\n",nkeys,keys[i].keyval);
	    break;
	  }
	}

        keys[i].keyval = defv[j];
        keys[i].key = scopy(parname(defv[j]));
        keys[i].val = scopy(parvalue(defv[j]));
	if (nzeno)
	  keys[i].help = scopy(defv[j+1] + 1);
	else
	  keys[i].help = scopy(parhelp(defv[j]));
        keys[i].count = 0;
        keys[i].upd = 1;
	if (keys[i].key[strlen(keys[i].key)-1] == '#')
	  keys[i].indexed = -1;   /* indexed, but not initialized */
	else
	  keys[i].indexed = -2;   /* not indexed */
	keys[i].next = NULL;
        if (streq(keys[i].key,"VERSION")) {      /* special (last?) keyword */
            version_i = scopy(keys[i].val);            
            keys[i].upd = 0;
            maxkeys--;
        }
    }



#if defined(KEYVAL)
    posflag = FALSE;                            /* force key=val args       */
#else
    posflag = TRUE;                             /* start scan by position   */
#endif
    for (i = 1; argv[i] != NULL; i++) {         /* loop over stuff in argv  */
#if defined(DASHDASH)
        if (streq(argv[i],"--")) {
	  dd_argv = &argv[i];
	  dd_argc = 1;
	  while (argv[++i])
	    dd_argc++;
	  dprintf(2,"DASHDASH: dd_argc=%d\n",dd_argc);
	  break;
        }
#endif
        name = parname(argv[i]);                /*   get param name, if any */
        posflag = posflag && (*name == 0);      /*   see how to match args  */
        if (posflag) {                          /*   match by position?     */
            getparam_argc++;			/* count keywords */
            if (i > maxkeys) error("Too many un-named arguments");
            if (streq(argv[i],"-help") || streq(argv[i],"--help") ||
		streq(argv[i],"-h") || streq(argv[i],"help")) {
	      set_help("");
            } else if (streq(argv[i],"-version") || streq(argv[i],"--version")) {
	      printhelp("V");    /* will also exit */
            } else {
	      if (keys[i].indexed >= -1) 
		error("Cannot match indexed keywords by position");
	      free(keys[i].val);
	      keys[i].val = scopy(argv[i]);
	      keys[i].count++;
            }
        } else {                                /*   match by name?         */
            if (name == NULL)
                error("Arg \"%s\" must be named", argv[i]);
            j = findkey(name);                  /*     find index in keys   */
            if (j >= 0) {                       /*     was name found?      */
	      getparam_argc++;			/* count keywords */
	      if (keys[j].count)
		error("Parameter \"%s\" duplicated", name);
	      free(keys[j].val);
	      keys[j].val = scopy(parvalue(argv[i]));        /* get value */
	      keys[j].count++;
	    } else if (j=set_indexed(name,&idx)) {       /* enter indexed keywords */
	      // process this indexed keyword
#if 1
	      addindexed(j,argv[i],idx);
#else
       	      dprintf(1,"Entering indexed %s, j=%d\n",argv[i],j);
	      kw = &keys[j];       /* start at the indexed (base) keyword */
	      while (kw->next) {        /* link through all previously indexed ones */
		dprintf(1,"Link List Skipping %s\n",kw->key);
		kw = kw->next;
		if (kw->indexed == idx+1) 
		  error("Duplicated indexed keyword %s",name);
	      }
	      kw->next = (keyword *) allocate(sizeof(keyword));
	      kw = kw->next;
	      kw->keyval = scopy(argv[i]);
	      kw->key = scopy(parname(argv[i]));
	      kw->val = scopy(parvalue(argv[i]));
	      kw->help = 0;
	      kw->count = 1;   // indexed keywords are never 'original'
	      kw->upd = 0;
	      kw->indexed = idx+1;
	      kw->next = NULL;
#endif
            } else {                            /*     not listed in defv?  */
	      if (streq(name, "help"))       /* got system help keyword?  */
		set_help(parvalue(argv[i]));
	      else if (streq(name, "debug"))   /*    got debug keyword? */
		set_debug(parvalue(argv[i]));
	      else if (streq(name, "yapp"))   /*     got yapp keyword?  */
		set_yapp(parvalue(argv[i]));
	      else if (streq(name, "review"))   /*    go into review ?  */
		set_review(parvalue(argv[i]));
	      else if (streq(name, "tcl"))      /*    go into tcl ?  */
		set_tcl(parvalue(argv[i]));
	      else if (streq(name, "error")) /*    allowed error count? */
		set_error(parvalue(argv[i]));
	      else if (streq(name, "outkeys")) /*    output keywords ? */
		set_outkeys(parvalue(argv[i]));
	      else
		error("Parameter \"%s\" unknown", name);
            } /* j>0 (i.e. program/system keyword) */
        } /* if(posflag) */
    } /* for (i) */

    /* do we need to check version_i vs. version_e here ??? */

#if defined(INTERRUPT)
    if (review_flag) {
        signal(SIGTERM, review);                /* catch interrupts */
        signal(SIGQUIT, review);                /* for review section */
        signal(SIGINT, review);                 /* and ^C also */
    }
#endif

#if defined(INTERACT)
    sprintf(&key_filename[strlen(key_filename)],"%s.def", progname);
    dprintf(2,"Keyword file will be %s\n",key_filename);

    if (help_level & HELP_DEFIO) {           /* read any existing keys-file */
	readkeys("initparam(1)",TRUE);
    } 
    if (help_level & HELP_GLOBAL) {
        warning("HELP_GLOBAL level not implemented");
    }

#endif /* INTERACT */



    if (help_level & HELP_PROMPT) {
#if defined(INTERACT)
        string key, val;
        int    go=0;            /* set immediate go to false */
        
        if (!isatty(fileno(stdin)) || !isatty(fileno(stdout)))
            error("(initparam) Can only INTERACT with terminal");

        for (i=1; i<nkeys; i++) {               
            if (streq(keys[i].key,"VERSION"))
                continue;   /* could possibly also set version_e here .. */
            do {                                 /*       loop for input     */
                if (streq(keys[i].val, "???")) {  /* still unanswered ? */
                    defflag = TRUE;
                    sprintf(keybuf,"%s=",keys[i].key);
                    go = readparam(keybuf, "No default allowed");  /* get it */
                } else {                  /* keyword with default value */
                    defflag = FALSE;
                    sprintf(keybuf,"%s=%s",keys[i].key,keys[i].val);
                    go = readparam(keybuf, "");     /* get it */
                }
                dprintf(0,"KEYBUF: %s\n",keybuf);
                if (go) break;          /* exit loop when premature done */
                key = parname(keybuf);
                val = parvalue(keybuf);
            } while ((key == NULL) ||
                     ((val == NULL || streq(val, "")) && defflag) );
            if (!streq(keys[i].val,val)) {
                free(keys[i].val);
                keys[i].val = scopy(val);
                keys[i].count++;
            }
            if (go) break;
        } /* for (i) */
#else   
    error("(initparam) INTERACTIVE input not allowed; check help=?");
#endif /* INTERACT */     
    }

#if defined(INTERACT)
    if (help_level&HELP_DEFIO || help_level&HELP_EDIT) {     /* write keyfile */
	writekeys("initparam(1)");
    }
    if (help_level&HELP_EDIT) {    /* into editor mode, edit and read keyfile */
        char *edtcmd;

        edtcmd = getenv("EDITOR");
        if (edtcmd != NULL)
            sprintf(keybuf,"%s %s",edtcmd,key_filename);       /* your editor */
        else
            sprintf(keybuf,"vi %s",key_filename);           /* default editor */

        if (system(keybuf))                         /*  edit the keyword file */
            error("(initparam) Error running editor (%s)",keybuf);
                        
	readkeys("initparam(2)",(bool)FALSE);
    }

    if (help_level&HELP_GLOBAL) {                /* write global keyword file */
        warning("HELP_GLOBAL: Not implemented");
    }
#endif /* INTERACT */

    if (help_string) {
      /* also patch printhelp if you add characters to this strpbrk check */
      if (      strpbrk(help_string,"oiapdqntkvhmc?")!=NULL ||
              ( strpbrk(help_string,"oiapdqntkvhmc?")==NULL && 
                strpbrk(help_string,"0123456789")==NULL
		) )
	printhelp(help_string);     /* give some help and possibly */
    }

    useflag = FALSE;
    for (i=1; i<nkeys; i++)               /* check if any left unresolved */
          useflag = useflag || streq(keys[i].val,"???");

    if (useflag) {
        printusage(defv);         /* give minimum 'usage' to stderr and exit */
        local_exit(0);
        /*NOTREACHED*/
    }
#if defined(INTERRUPT)
    initparam_done = 1;                 /* flag for interruptability */
    review();
    beep();
#endif
    save_history(argv);                 /* save activation record */
    eval_keys();	            	/* resolve $key references */
#if defined(TCL)
    if (tcl_flag) {
        ini_readline();
        run_tcl();
        if (help_level&HELP_DEFIO) writekeys("initparam(2)");
        end_readline();
    }
        
#endif    

#if 0
    if (debug_level > 0) {              /* version check : via external */
        i = findkey("VERSION"); 
        if (i>0) {
            sprintf(keybuf,"nemo.version %s %s\n",keys[0].val,keys[i].val);
            system(keybuf);
        } else
            warning("No VERSION keyword");
    }
#endif
    initparam_out();
#ifndef __MINGW32__
    clock1 = times(&tms1);
#endif
} /* initparam */

local void initparam_out()
{
  int nout = -1;
  if (outdefv && *outdefv)
     nout = xstrlen(outdefv, sizeof(string));     /* count # params + progname */  
  /*  
   * here comes 
   *      - parsing code for outdefv[]         'key=[n]t\n    help'
   *      - build linear list of this
   */

  dprintf(1,"Found %d output keywords\n",nout);
}


/*
 * FINIPARAM: finish up any outstanding requests before the program exists
 *            (users can use stop() to exit and bypass this routine)
 */

void finiparam()
{
    int i, n=0;
#ifndef __MINGW32__
    if (report_cpu) report('c');
    if (report_mem) report('m');
#endif
    for (i=1; i<nkeys; i++)
        n += keys[i].upd ? 1 : 0;

    /* TODO:::   fix this for indexed keywords */

    if (n && debug_level > 0) {
	dprintf(1,"There were %d parameters used on the commandline\n",
		getparam_argc);
        warning("(finiparam) The following %d keywords have never been read:",n);
        for (i=1; i<nkeys; i++)
            if (keys[i].upd) dprintf(1," %s ",keys[i].key);
        dprintf(1,"\n");            
    }
#if defined(INTERACT)

    if (help_level&HELP_DEFIO) {
        dprintf(2,"finiparam: writing keyword file for final time\n");
	writekeys("finiparam");

    }
    if (help_level & HELP_GLOBAL) {
        warning("HELP_GLOBAL: Not implemented yet");
    }

#if defined(INTERRUPT)
    review();
#endif

#endif

    finiparam_out();

    /* free up some junk malloc checkers complain about */

    dprintf(1,"finiparam: now freeup some final memory\n");
    reset_history();
    free(yapp_string);
    for (i=0; i<nkeys; i++) {
      free(keys[i].key); 
      free(keys[i].val);
      free(keys[i].help);
    }
    free(keys);
    if (version_i) free(version_i);
}


local void report(char what)
{
  float ticks;

  if (what == 'c') {
#ifndef __MINGW32__
    clock2 = times(&tms2);
    ticks = (float) sysconf(_SC_CLK_TCK);
    dprintf(0,"CPU_USAGE %s : %.2f    %.2f %.2f  %.2f %.2f  %ld\n",
	    progname, 
	    (clock2-clock1)/ticks,
	    (tms2.tms_utime-tms1.tms_utime)/ticks,
	    (tms2.tms_stime-tms1.tms_stime)/ticks,
	    (tms2.tms_cutime-tms1.tms_cutime)/ticks,
	    (tms2.tms_cstime-tms1.tms_cstime)/ticks,
	    clock1);
#endif
    return;
  }

  if (what == 'm') {
#if defined(linux)
    /* http://www.gnu.org/software/libc/manual/html_node/Statistics-of-Malloc.html */
    struct mallinfo mi = mallinfo();
    dprintf(0,"mallinfo: hblks(d):%d %d uord=%d ford=%d keepcost=%d arena=%d ord=%d\n",
	    mi.hblks,mi.hblkhd,mi.uordblks,mi.fordblks,mi.keepcost,mi.arena,mi.ordblks);
#else
    dprintf(0,"report_mem not implemented for non-linux\n");
#endif      
    return;
  }
  return;
}

/*
 * nemokeys:  set the location of the nemokeys file, since it depends
 *            on the value of NEMODEF
 */

local string nemokeys(void)
{
  static char path[256];
  char *nemodef;
  string fn;

  nemodef = getenv("NEMODEF");
  if (nemodef)
    sprintf(path,"%s/nemokeys.dat",nemodef);
  else
    sprintf(path,"./nemokeys.dat");
  dprintf(1,"NEMOKEYS: using %s\n",path);
  return path;
}

local void finiparam_out(void)
{
  /*
   * this will
   *    - read the $NEMODEF/nemokeys.dat file
   *    - replace values from outdefv[] 
   *    - write the $NEMODEF/nemokeys.dat file
   * Also needs to take care of file locking ..
   */
  string fn = nemokeys();
  int fd, err;
  pid_t pid = getpid();
  char msg[128];
  sprintf(msg,"%s(%d)\n",progname,pid);
#if 0
  dprintf(1,"writing to %s\n",fn);

  /* if it exists */
  fd = open(fn,O_CREAT|O_RDWR|O_APPEND,S_IRUSR|S_IWUSR);
  if (fd<0) local_error("finiparam_out: can't open");
  err = flock(fd,LOCK_EX);
  if (err < 0) local_error("finiparam_out: can't lock");
  err = write(fd,msg,strlen(msg));
  if (err < 0) local_error("finiparam_out: can't write");
  /* should use fcntl() really, for NFS */
  err = flock(fd,LOCK_UN);
  if (err < 0) local_error("finiparam_out: can't unlock");
  err = close(fd);
  if (err < 0) local_error("finiparam_out: can't close");
#endif
}




/*
 * scan the environment variables, and set the ones NEMO wants to
 * know about. This is where future expansion is system and private
 * lookup tables should go. (ala AIPS++)
 */
local void scan_environment()
{
    int i;
    char *ev;

    for (i = 0; environ[i] != NULL; i++) {
    	/* printf("%d: \"%s\"\n",i+1,environ[i]); */
        if (streq("BELL", parname(environ[i])))
            bell_level = atoi(parvalue(environ[i]));
        else if (streq("HISTORY", parname(environ[i])))
            nemo_history = atoi(parvalue(environ[i]));
        else if (streq("DEBUG", parname(environ[i])))
            set_debug(parvalue(environ[i]));
        else if (streq("YAPP", parname(environ[i])))
            set_yapp(parvalue(environ[i]));
        else if (streq("HELP", parname(environ[i]))) 
            set_help(parvalue(environ[i]));
        else if (streq("REVIEW", parname(environ[i])))
            set_review(parvalue(environ[i]));
        else if (streq("ERROR", parname(environ[i])))
            set_error(parvalue(environ[i]));
        else if (streq("TCL", parname(environ[i])))
            set_tcl(parvalue(environ[i]));
    }
    dprintf(5, 
    "scan_environment: debug=%d yapp=%d help=%d history=%d review=%d error=%d\n",
     debug_level, yapp_dev, help_level, nemo_history, review_flag, error_level);
    dprintf(5,"date_id = %s\n",date_id());
#if defined(INTERACT)
    if ((ev=getenv("NEMODEF")) != NULL) {
        strcpy (key_filename,ev);
        strcat (key_filename,"/");
    } else
        key_filename[0] = '\0';
#endif /* INTERACT */
}

local void save_history(string *argv)
{
#if defined(HISTORY)
    int i, histlen;
    string hist;

    if (!nemo_history) return;

    if (help_level) {         /* if external help, argv[] is not appropriate */
        histlen = 0;
        for (i=0; i<nkeys; i++)            /* count how much needed */
            histlen += strlen(keys[i].key) + strlen(keys[i].val) + 2;
        hist = (char *) allocate(histlen+20);   /* extra 20 for help=.... */
        strcpy(hist,keys[0].val);
        for (i=1; i<nkeys; i++) {
            strcat(hist," ");
            strcat(hist,keys[i].key);
            strcat(hist,"=");
            strcat(hist,keys[i].val);
        }
        sprintf(&hist[strlen(hist)]," help=%d",help_level);     /* add help */
    } else {                /* save the argv[] command line */
        histlen = strlen(progname) + 11;
	if (version_i) histlen += strlen(version_i);
        for (i = 1; argv[i] != NULL; i++)
            histlen += strlen(argv[i]) + 1;
        hist = (char *) allocate(histlen);
        strcpy(hist,progname);
        for (i = 1; argv[i] != NULL; i++) {     /* copy name and arguments  */
            strcat(hist, " ");
            strcat(hist, argv[i]);
        }
        strcat(hist," VERSION=");            /* append version identifcation */
	if (version_i) strcat(hist, version_i);
    }
    app_history(hist);
    free(hist);
#endif /* HISTORY */
}

/*
 * PRINTHELP: print out help message in various ways.
 */

local void printhelp(string help)
{
    string *bp;
    int i, numl;

    dprintf(1,"printhelp: help_string=%s\n",help);

    if (strchr(help,'?')) {
        printf("Help options are any combination of:\n\n");
        printf("  a       key=val   (all)\n");
        printf("  p,k     key       (parameter,key)\n");
        printf("  d,v     val       (default,value)\n");
        printf("  n       >> add's a newline between above\n");
        printf("  q       >> quit's after above work done\n");
        printf("  u       >> show the usage string\n");
	printf("  h	  >> show key and help strings\n");
        printf("  t       >> show as MIRIAD doc file\n");
        printf("  z       >> show as KHOROS pane file\n");
        printf("  i       >> show some internal variables\n");
	printf("  o       >> show the output key names\n");
	printf("  c       >> show cpu usage at the end of the run\n");
	printf("  m       >> show memory usage at the end of the run\n");
	printf("  I       >> cvs id\n");
        printf("  ?       >> this help (always quits)\n\n");
        printf("Numeric helplevels determine degree and type of assistence:\n");
        printf("They can be added to give combined functionality\n");
        printf("  1       keyword files used (prognam.def)\n");
        printf("  2       interactive prompting\n");
        printf("  4       menu interface using EDITOR environment variable\n");
	printf("  8       reserved\n");
	printf(" 16       reserved\n");
        printf(" VERSION_ID = %s\n",GETPARAM_VERSION_ID);
        printf(" NEMO VERSION = %s\n",NEMO_VERSION);
        showconfig();
	showsystem();
        local_exit(0);
        /*NOTREACHED*/
    }

    if (strchr(help,'i')) {
    	printf("NEMO version: %s\n",NEMO_VERSION);
        printf("help: %s yapp: %s error: %s\n",
                help_string  ? help_string  : "",
		yapp_string  ? yapp_string  : "",
		error_string ? error_string : "");
        printf("debug_level=%d error_level=%d\n",
                debug_level,   error_level);
        printf("argv: %s\n",
                argv_string  ? argv_string  : "");
    }
    if (strchr(help,'V')) {
        for (i=1; i<nkeys; i++)
            if (streq(keys[i].key,"VERSION")) 
                printf("%s  %s (%s)\n",
                       keys[0].val,keys[i].val,keys[i].help);

        local_exit(0);
        /*NOTREACHED*/
    }
    if (strchr(help,'I')) {
      printf("%s\n",cvsid);
        local_exit(0);
        /*NOTREACHED*/
    }

    if (strchr(help,'h')) {
        for (i=1; i<nkeys; i++)
            printf("%-16s : %s [%s]\n",keys[i].key, keys[i].help, keys[i].val);
        local_exit(0);
        /*NOTREACHED*/
    }

    numl = ((strchr(help,'n')) ? 1 : 0);    /* add newlines between key=val ? */

    if (strchr(help,'a') || strpbrk(help,"oapdqntvkzucm")==NULL) { /* arguments */
        printf("%s", progname);
        for (i=1; i<nkeys; i++) {
            newline(numl);
            printf(" %s=%s",keys[i].key, keys[i].val);
        }
        newline(1);
        if (strpbrk(help,"oapdqntvkzu")==NULL) {
	  local_exit(0);
            /*NOTREACHED*/
        }
    }

    if (strchr(help,'p') || strchr(help,'k')) {      /* parameters keywords */
        printf("%s", progname);
        for (i=1; i<nkeys; i++) {
            newline(numl);
            printf(" %s", keys[i].key);
        }
        newline(1);
    } /* 'p', 'k' */

    if (strchr(help,'d') || strchr(help,'v')) {              /* arguments */
        printf("%s", progname);
        for (i=0; i<nkeys; i++) {
            newline(numl);
            printf(" %s", keys[i].val);
        }
        newline(1);
    } /* 'd' 'v' */
    
    if (strchr(help,'u')) {
        printf("%s\n",usage);
    } /* 'u' */

    if (strchr(help,'o')) {
      string *sp = outdefv;

      warning("New option help=o in the user interface\nUnformatted output");
      if (sp) {
	while (*sp) {
	  printf("%s\n",*sp++);
	}
      } else {
	warning("No output keys defined for this program");
      }
      local_exit(0);
      /*NOTREACHED*/
    }

    if (strchr(help,'t')) {     /* nemotool - special doc file output */
        printf("%%N %s\n",progname);
        printf("%%D %s\n",usage);
        printf("%%B\n");
        printf("  This doc file has been produced with NEMO help=t option\n");
        printf("  Try 'man %s' for more extensive online help\n",progname);
	printf("  Defaults of keywords are given between square brackets\n");
        for (i=1; i<nkeys; i++) 
            printf("%%A %s\n\t%s [%s]\n",
                    keys[i].key,
                    keys[i].help ? keys[i].help : "No help",
                    keys[i].val);
        local_exit(0);                    /* always exit in this case */
        /*NOTREACHED*/
    } /* 't' */

    if (strchr(help,'z')) {     /* khoros - special pane file output */
        int vcount = 2, opt;
        char *parval;

        printf("-F 4.2 1 0 170x7+10+20 +35+1 'CANTATA for KHOROS' cantata\n");
        printf("-M 1 0 100x40+10+20 +23+1 'A NEMO program' nemo\n");
        printf("-P 1 0 80x38+22+2 +0+0 '%s ' %s\n",usage,progname);
        for (i=1; i<nkeys; i++) {           /* help */
            parval = keys[i].val;
            opt = !streq(parval,"???");
            if (strncmp(keys[i].key,"in",2)==0)
                printf("-I 1 0 %d 1 0 1 50x1+2+%d +0+0 '%s' '%s ' '%s' %s\n",
                           opt,          
                           vcount,
                           (*parval == 0 ? " " : (opt?parval:" ")),
                           keys[i].key,keys[i].help,keys[i].key);

            else if (strncmp(keys[i].key,"out",3)==0)
                printf("-O 1 0 %d 1 0 1 50x1+2+%d +0+0 '%s' '%s ' '%s' %s\n",
                           opt,          
                           vcount,
                           (*parval == 0 ? " " : (opt?parval:" ")),
                           keys[i].key,keys[i].help,keys[i].key);
            else
                printf("-s 1 0 %d 1 0 50x1+2+%d +0+0 '%s' '%s ' '%s' %s\n",
                           opt,          
                           vcount,
                           (*parval == 0 ? " " : (opt?parval:" ")),
                           keys[i].key,keys[i].help,keys[i].key);
            vcount += 2;
        }
        vcount++;
        printf("-H 1 13x2+1+%d 'Help' 'Help for %s' nemo.help\n",
                            vcount, progname);
        printf("-R 1 0 1 13x2+39+%d 'Run' 'RunMe' khoros2nemo %s\n",
                            vcount, progname);
        printf("-E\n-E\n-E\n");
        local_exit(0);                    /* always exit in this case */
        /*NOTREACHED*/
    } /* 'z' */
    if (strchr(help,'q')) {
        local_exit(0);                        /* quit - don't run program */
        /*NOTREACHED*/
    }
    if (strchr(help,'c')) {
        report_cpu = 1;
    }
    if (strchr(help,'m')) {
        report_mem = 1;
    }
}

local void newline(int n)        /* print 'n' newlines */
{
    while (n-- > 0) printf("\n");
}

/*
 * SHOWCONFIG: show what non-standards have been compiled into the
 *             the interface
 */

local void showconfig()
{
    printf("INTERACT   ");
#if defined(INTERACT)
    printf("on\n");
#else
    printf("off\n");
#endif
    printf("NEMOINP    ");
#if defined(NEMOINP)
    printf("on\n");
#else
    printf("off\n");
#endif
    printf("PUTPARAM   ");
#if defined(PUTPARAM)
    printf("on\n");
#else
    printf("off\n");
#endif
    printf("MACROREAD  ");
#if defined(MACROREAD)
    printf("on\n");
#else
    printf("off\n");
#endif
    printf("HISTORY    ");
#if defined(HISTORY)
    printf("on\n");
#else
    printf("off\n");
#endif
    printf("INTERRUPT  ");
#if defined(INTERRUPT)
    printf("on\n");
#else
    printf("off\n");
#endif
    printf("MINMATCH   ");
#if defined(MINMATCH)
    printf("on\n");
#else
    printf("off\n");
#endif
    printf("TCL        ");
#if defined(TCL)
    printf("on\n");
#else
    printf("off\n");
#endif
    printf("INDEXED    ");
#if defined(INDEXED)
    printf("on\n");
#else
    printf("off\n");
#endif
    printf("OUTKEYS    ");
#if defined(OUTKEYS)
    printf("on (testing)\n");
#else
    printf("off\n");
#endif

    printf("Precision: %s\n",Precision);
}
/*
 * SHOWSYSTEM: show the system keywords
 */

local void showsystem()
{
  printf("System keywords:\n");
  printf("  help=     various inline help options (see above) for program keywords\n");
  printf("  debug=    change the debug level [0]\n");
  printf("  error=    change the number of allowed fatal errors [0]\n");
  printf("  yapp=     change the default plotting device [$YAPP]\n");
  printf("  outkeys=  **under development**  outkeys to be written\n");
  printf("  review=   interrupt mode to review keywords before start\n");
  printf("  argv=     addition cmdline arguments not parsed by NEMO\n");
  printf("  tcl=      go into tcl (deprecated)\n");
}
/*
 * PRINTUSAGE: print out helpful usage info.
 */

local void printusage(string *defv)
{
    int i;
    bool otherargs;

    if(mpi_proc) fprintf(stderr,"@%d: ",mpi_rank);
    fprintf(stderr,"Insufficient parameters, try 'help=', 'help=?' or 'help=h',\n");
    if(mpi_proc) fprintf(stderr,"@%d: ",mpi_rank);
    fprintf(stderr,"Usage: %s", progname);
    otherargs = FALSE;
    for (i=1; i<nkeys; i++)
        if (streq(keys[i].val,"???"))
	  fprintf(stderr," %s=???", keys[i].key);
        else
	  otherargs = TRUE;
    fprintf(stderr,otherargs ? " ...\n" : "\n");
    if (usage) fprintf(stderr,"%s\n",usage);
}

/*
 * GET_MACRO:  see if string is of a '@file', if so, read it into
 *             newly allocated memory and return pointer to that
 *             string;
 *             if not, return string itself.
 *   Note: not hierarchical/recursive
 *         @ must be first character of the keyword 'value' part
 *         does not expand $ variables
 */
local string get_macro(char *mname)
{
#if defined(MACROREAD)
    char *cp, *mp;
    int n;
    stream fstr;    

    if (*mname != '@')  return mname;
    mname++;
    dprintf(1,"getparam[get_macro]: Opening macro file %s\n",mname);
    n = nemo_file_size(mname);
    if (n<0)
        error("(getparam) macro include file \"%s\" does not exist",mname);
    mp = (char *) allocate(n+1);        /* extra 1 for 0 at the end */
    if (n==0) {              /* make sure n=0 also works properly */
        *mp = '\0';
        return mp;
    }

    fstr = stropen(mname,"r");      /* open macro file (guaranteed ok)     */
    fread(mp,sizeof(char),(unsigned)n,fstr);    /* oops, unix only here ?? */
    strclose(fstr);                 /* file is read in 'as is' in binary   */
    mp[n] = 0;  /* terminate - since mp could point anywhere */

    for (cp=mp; *cp; cp++) {
      if (*cp=='\n') {
	if (cp[1]==0)
	  *cp = 0;            /* last line: newline trimmed */
	else
	  *cp=' ';            /* replace newline with blanks */
      }
    }
    return mp;
#else
    if (*mname == '@') 
      warning("getparam.c: macro not compiled (use: -DMACROREAD)");
    return mname;              /* no processing, return string */
#endif
}

/*
 * EVAL_KEYS: scan key=val pairs for cross references
 *            should do $ as well as @ parsing - only does $ now
 */
 
local void eval_keys(void)
{
    int i, j;
    char *c1,*c2,*c3,*c4, *cd, check, newval[MAXBUF], key[MAXKEYLEN];
    bool more = TRUE;

    while (more) {                   /* loop until no more $ references found */
        more = FALSE;
        for (i=1; i<nkeys; i++) {                  /* loop over all key=val's */
            c1 = keys[i].val;   /* from */
            c2 = newval;     /* to */
            cd = strchr(c1,'$');                   /* see if cross ref used */
            if (cd == 0) continue;
            dprintf(1,"eval_keys: parsing %s\n",c1);
            more = TRUE;
            while (cd) {                       /* loop, while more '$' found */
                if (cd[1] == '{')       check = '}';
                else if (cd[1] == '(')  check = ')';
                else                    check = 0;
                while (*c1 && *c1 != '$') *c2++ = *c1++;    /* copy before $ */
                if (*c1 == 0)
                    error("improper $-reference in %s",keys[i].val);

                c3 = key;
                c1++;
                if (check) {
                    c1++;
                    while (*c1 && *c1 != check) *c3++  = *c1++;
                    c1++;
                } else
                    while (*c1 && (isalpha(*c1) || *c1=='_')) *c3++ = *c1++;
                *c3 = 0; /* terminate key */
                j = findkey(key);
                if (j<0) {
                    c4 = getenv(key);
                    if (c4==NULL) 
                        error("eval_keys: bad $%s in %s=%s",
                            key,keys[i].key,keys[i].val);
                } else if (i==j) {
                    error("eval_keys: recursive %s=%s",keys[i].key,keys[i].val);
                    c4 = "";
                } else
                    c4 = keys[j].val;
                dprintf(3,"Patching %s with %s=%s\n",keys[i].val,key,c4);
                while (*c4) *c2++ = *c4++;
                cd = strchr(c1,'$');        /* search for next */
            }
            while (*c1) *c2++ = *c1++;
            *c2 = 0;
            dprintf(3,"eval_keys: Modifying %s\n",keys[i].val);
            free(keys[i].val);
            keys[i].val = scopy(newval);
            /* keys[i].count++; */
        }
    }
}


/* GETARGV:  find out if any parameters were unparsed (only if DASHDASH enabled)
 *
 */

string *getargv(int *argc)
{
  *argc = dd_argc;
  return dd_argv;
}
/*
 * CNTPARAM: return number of parameters entered on commandline 
 */
int cntparam(void)
{
    return getparam_argc;
}
/*
 * ISAPARAM: check if this is truely a parameter for this program
 *           The only way to prevent program from crashing through error(3NEMO)
 */

bool isaparam(string name)
{
    if (nkeys==0) local_error("isaparam: called before initparam");
    return findakey(name) != NULL;
}

/*
 * HASVALUE:  check if a parameter has been given a value
 */

bool hasvalue(string name)     
{
  char *cp, key[MAXKEYLEN];
  int n = strlen(name);
  int idx;
  keyword *kw;

  strcpy(key,name);
  dprintf(2,"hasvalue: checking indexing on %s\n",key);

  kw = findakey(name);
  if (kw == NULL) error("keyword %s does not exist",name);

  if (kw->val == NULL) return FALSE;
  return  strlen(kw->val) > 0;
}

/*
 * UPDPARAM:  check if a parameters value has been updated since last read
 */

bool updparam(string name)
{
    keyword *kw = findakey(name);
    if (kw == NULL) error("(updparam) \"%s\" unknown keyword",name);
    return kw->upd == 1;
}

/*
 *  GETPARAMSTAT
 */

int getparamstat(string name) {
  error("getparamstat is a ZENO feature, not implemented in NEMO yet");
}



/*
 * GETPARAM: look up parameter value by name
 */

string getparam(string name)
{
    char *cp, *cp1;
    keyword *kw;
    int i;

    if (nkeys == 0) local_error("(getparam) called before initparam");

    kw = findakey(name);
    if (kw == NULL) error("(getparam) \"%s\" unknown keyword", name);
    kw->upd = 0;
#if defined(MACROREAD)
    cp = kw->val;
    if (*cp == '@') {
        cp1 = kw->val;
        kw->val = get_macro(cp);
        free(cp1);
    }
#endif
    return kw->val;                     /* return value part */
}

/* 
 * cost:  10M in 5" on a P600 (w/ gcc_2.95.3 -g) 
 */

string getparam_idx(string name, int idx)
{
    int i;
    char *cp, *cp1, key[MAXKEYLEN+1];
    keyword *kw;

    if (nkeys == 0)  local_error("(getparam_idx) called before initparam");
    strcpy(key,name);
    strcat(key,"#");
    i = findkey(key);
    if (i<0) error("(getparam_idx) \"%s\" unknown keyword", name);
    kw = &keys[i];
    if (kw->indexed < -1) error("%s is not an indexed keyword",name);
    while (kw->next) {
      dprintf(1,"Checking linked list w/ %s for %d, %d\n",kw->key,kw->indexed,idx);
      kw = kw->next;
      if (kw->indexed == idx+1) {
#if defined(MACROREAD)
	cp = kw->val;
	if (*cp == '@') {
	  cp1 = kw->val;
	  kw->val = get_macro(cp);
	  free(cp1);
	}
#endif
	return kw->val;
      }
    }
    return NULL;
}

/*
 * indexparam: 
 *       if IDX < 0, it will return 
 *           >=0  the largest index available in an 0-based indexed parameter
 *            -1  if no indexed one for this keyword has been used yet
 *            -2  if the parameter is not indexed at all
 *
 *       if IDX >= 0 is given, it will check for that index and returns:
 *             1  if that index exists
 *             0  if that index does not exist
 */

int indexparam(string basename, int idx)
{
  int i, idxmax = -1;
  char *cp, key[MAXKEYLEN+1];
  keyword *kw;

  if (nkeys == 0)  local_error("(indexparam) called before initparam");
  strcpy(key,basename);
  strcat(key,"#");
  i = findkey(key);        /* first find if it's indexed */
  if (i<0) {
    i = findkey(basename);     /* if not, see if it is indexed at all */
    if (i < 0) return -2;      /* no, not indexed, return -2 */
  }
  kw = &keys[i];
  if (kw->indexed < -1) error("%s is not an indexed keyword",basename);
  if (idx < 0) {        /* check the max idx and return it */
    while (kw->next) {
      kw = kw->next;
      idxmax = MAX(idxmax, kw->indexed);
    }
    return idxmax;
  } else {              /* else check if a specific index exists */
    while (kw->next) {
      kw = kw->next;
      if (kw->indexed == idx+1) return 1;
    }
    return 0;
  }
  error("indexparam: reached an impossible part of code");
  return -3;   /* never reached */
}

/*
 * GETIPARAM, ..., GETDPARAM: get int, long, bool, or double parameters.
 *          and allow some optional parsing if -DNEMOINP
 */

int getiparam(string par)
{
    string val;
    int    ipar, nret;

    val = getparam(par);                        /* obtain value of param */
#if !defined(NEMOINP)
    return (atoi(val));                         /* convert to an integer */
#else
    if (strncmp("0x",val,2)==0)                 /* parse as single hex # */
      return strtol(val,0,16);
    nret = nemoinpi(val,&ipar,1);
    if (nret < 0)
        error("getiparam(%s=%s) parsing error %d, assumed %d\n",
                    par,val,nret,ipar);
    return (nret==0) ? 0 : ipar;
#endif /* !NEMOINP */
}

int getiparam_idx(string par, int idx)
{
    string val;
    int    ipar, nret;

    val = getparam_idx(par,idx);                /* obtain value of param */
#if !defined(NEMOINP)
    return atoi(val);                           /* convert to an integer */
#else
    nret = nemoinpi(val,&ipar,1);
    if (nret < 0)
        error("getiparam_idx(%s=%s,%d) parsing error %d, assumed %d\n",
                    par,val,idx,nret,ipar);
    return (nret==0) ? 0 : ipar;
#endif /* !NEMOINP */
}

long getlparam(string par)
{
    string val;
    long lpar;
    int nret;

    val = getparam(par);                        /* obtain value of param */
#if !defined(NEMOINP)
    return atol(val);                           /* convert to an integer */
#else
    if (strncmp("0x",val,2)==0)                 /* parse as single hex # */
      return strtol(val,0,16);
    nret = nemoinpl(val,&lpar,1);
    if (nret < 0)
        error("getlparam(%s=%s) parsing error %d assumed %l\n",
                    par,val,nret,lpar);
    return (nret==0) ? 0 : lpar;
#endif /* !NEMOINP */
}

long getlparam_idx(string par, int idx)
{
    string val;
    long lpar;
    int nret;

    val = getparam_idx(par,idx);                        /* obtain value of param */
#if !defined(NEMOINP)
    return atol(val);                           /* convert to an integer */
#else
    nret = nemoinpl(val,&lpar,1);
    if (nret < 0)
        error("getlparam(%s=%s) parsing error %d assumed %l\n",
                    par,val,nret,lpar);
    return (nret==0) ? 0 : lpar;
#endif /* !NEMOINP */
}

bool getbparam(string par)
{
    char *val;
    bool bpar;
    int  nret;

    val = getparam(par);                        /* obtain value of param */
#if !defined(NEMOINP)
    if (*val=='.') val++;                       /* catch .TRUE. .FALSE. */
    if (strchr("1tTyYjJ", *val) != NULL)        /* is value true? */
        return TRUE;
    if (strchr("0fFnN", *val) != NULL)          /* is value false? */
        return FALSE;

    error("getbparam: %s=%s not bool", par, val);
    return 0;   /*turboc*/
#else
    nret = nemoinpb(val,&bpar,1);
    if (nret < 0)
        error("getbparam(%s=%s) parsing error %d, assumed %d (FALSE)",
                        par,val,nret,bpar);
    return (nret<=0) ? FALSE : bpar;
#endif /* !NEMOINP */
}

bool getbparam_idx(string par, int idx)
{
    char *val;
    bool bpar;
    int  nret;

    val = getparam_idx(par,idx);                /* obtain value of param */
#if !defined(NEMOINP)
    if (*val=='.') val++;                       /* catch .TRUE. .FALSE. */
    if (strchr("tTyY1jJ", *val) != NULL)        /* is value true? */
        return TRUE;
    if (strchr("fFnN0", *val) != NULL)          /* is value false? */
        return FALSE;

    error("getbparam: %s=%s not bool", par, val);
    return 0;   /*turboc*/
#else
    nret = nemoinpb(val,&bpar,1);
    if (nret < 0)
        error("getbparam(%s=%s) parsing error %d, assumed %d (FALSE)",
                        par,val,nret,bpar);
    return (nret<=0) ? FALSE : bpar;
#endif /* !NEMOINP */
}

double getdparam(string par)
{
    string val;
    double dpar;
    int    nret;

    val = getparam(par);                        /* obtain value of param */
#if !defined(NEMOINP)
    return atof(val);                         /* convert to a double */
#else
    nret = nemoinpd(val,&dpar,1);
    if (nret < 0)
        error("getdparam(%s=%s) parsing error %d, assumed %g",
                            par,val,nret,dpar);
    return (nret==0) ? 0.0 : dpar;
#endif /* !NEMOINP */
}

double getdparam_idx(string par, int idx)
{
    string val;
    double dpar;
    int    nret;

    val = getparam_idx(par,idx);                        /* obtain value of param */
#if !defined(NEMOINP)
    return atof(val);                         /* convert to a double */
#else
    nret = nemoinpd(val,&dpar,1);
    if (nret < 0)
        error("getdparam(%s=%s) parsing error %d, assumed %g",
                            par,val,nret,dpar);
    return (nret==0) ? 0.0 : dpar;
#endif /* !NEMOINP */
}

float getfparam(string par)
{
    return (float) getdparam(par);
}

float getfparam_idx(string par, int idx)
{
    return (float) getdparam_idx(par,idx);
}

/*
 *
 *  PUTPARAM:       resetting the value of a keyword.
 *  PROMPTPARAM:    interactive setting the value of a keyword
 *
 *  both function call a local work horse:
 * 
 *  SETPARAM:
 *      
 *      prompt==NULL:  both 'par' and 'val' are assumed a valid string
 *                      and the keyword 'par' is updated with the value 'val'
 *      prompt!=NULL:  The string 'val' is ignored, the string 'prompt' is
 *                      a user prompt for interactive updating of the value
 *                      of the keyword 'par': the user must enter a new 'val'
 *                      If need be, the bindvec[] pointer to the 'par=val'
 *                      string will be reallocated.
 *
 *      setparam does not check if standard input is a redirected unit,
 *      in which case it's use is questionable or dangerous.
 *
 */
void putparam(string par,string val)
{

    setparam(par,val,(char *)NULL);
}

void promptparam(string par,string prompt)
{
    /* BUG: check if OK to do IO, i.e. isatty(stdin) && isatty(stdout).. */
    setparam(par,(char *)NULL,prompt);
}

local void setparam (string par, string val, string prompt)
{
#if defined(PUTPARAM) && defined(INTERACT)
    char line[80];      
    keyword *kw;
    int  i, idx;

    if (par==NULL || *par==0)
        error("setparam: no parameter supplied?");
    if (nkeys==0)
        local_error("setparam: called before initparam");

    kw = findakey(par);
    if (kw == NULL) {
      // check if it's an indexed keyword, we'll allow a new one to be entered
      char *keyval;
      i=set_indexed(par,&idx);
      if (i==0)
	error("setparam: parameter \"%s\" unknown",par);
      keyval = (char *) allocate(strlen(par)+strlen(val)+2);
      strcpy(keyval,par);
      strcat(keyval,"=");
      strcat(keyval,val);
      addindexed(i, keyval, idx);    /* create this new indexed key */
      free(keyval);
      return;
    }

    if (prompt && *prompt!=0) {        /* interactive entry */
        beep();
        fprintf (stderr,"%s: %s=",par, prompt);     /* DPRINTF ??? */
        fflush(stderr);
	clearerr(stdin);
#if 0
        if (gets(line) == NULL) error("Null input");
#else
        // if (fgets(line,80,stdin) == NULL) error("Null input");
	error("Can't do prompting anymore until fgets() is fixed");
#endif
        val = line;
    }
    kw->val = scopy(val);
    kw->upd = 2;              /* mark it as being updated */
#else
    error("setparam: not compiled into getparam.c");
#endif /* PUTPARAM && INTERACT */
}

/*
 * OUTPARAM:  a set of routines that store keyword values
 *            in the registered database of valid keywords meant
 *            for output
 *            For efficiency one local routine, writeparam_out
 *            does all the work talking to the database, is has
 *            a number of derived functions that deal with each
 *            type, to ensure type-safety.
 *
 *             
 */
local void writparam_out(string key, char type, int ndat, void *dat)
{
  /* find the entry that matches key
   * then check ndat and type
   * if so, copy dat into
   * for strings, need to free(), malloc() and strdup()
   */
  warning("writeparam_out(%s,%c,%d,0x%x) not implemented yet",key,type,ndat,dat);
}

void outparam(string key, string val)
{
  int n = strlen(val)+1;    /* count the terminating 0 as part of the data */
  writparam_out(key, 'c', n, val);
}

void outdparam(string key, double val)
{
  writparam_out(key, 'd', 1, &val);
}

void outfparam(string key, float val)
{
  writparam_out(key, 'f', 1, &val);
}

void outlparam(string key, long val)
{
  writparam_out(key, 'l', 1, &val);
}

void outiparam(string key, int val)
{
  writparam_out(key, 'i', 1, &val);
}

void outbparam(string key, bool val)
{
  writparam_out(key, 'b', 1, &val);
}



void outdparams(string key, int n, double *val)
{
  writparam_out(key, 'd', n, val);
}

void outfparams(string key, int n, float *val)
{
  writparam_out(key, 'f', n, val);
}

void outlparams(string key, int n, long *val)
{
  writparam_out(key, 'l', n, val);
}

void outiparams(string key, int n, int *val)
{
  writparam_out(key, 'i', n, val);
}

void outbparams(string key, int n, bool *val)
{
  writparam_out(key, 'b', n, val);
}


#if defined(INTERACT)
#define ESC 0x1b
/* 
 * READPARAM:  interactive input of a keyword
 *             normally returns 0 , but 1 if user wants global exit
 *             through some kind of HOT key (not really implemented)
 *
 *      BUG: not all can deal with ioctl() 
 */

local int readparam(string buffer,string text)
{
#if defined(TIOCSTI)
    char *p;
    int  i,n, go = 0;

    if (!isatty(fileno(stdin)))
        error("readparam: Cannot use redirected input in interactive mode");
    fflush(stdin);                  /* some programs seem to need this */
    beep();
    if (text && *text)
        printf ("%s\n",text);
    n = strlen(buffer);
    for (i=0; i<n; i++)             /* simulate as if this was already input */
        ioctl(fileno(stdin),TIOCSTI,&buffer[i]);
    p = buffer;
    p--;  
    do {                            /* read it */
        p++;
        *p = getchar();
        n++;
    } while (*p != '\n');
    *p = 0;                      /* terminate string */
    dprintf(1,"readparam: buffer=%d p=%d\n",buffer,p);
    return go;
#else
    error("readparam: feature not available (missing TIOCSTI)");
#endif
}
#endif /* INTERACT */


#define BELL 0x07

local void beep()
{   
    if (bell_level)
        putchar((char) BELL);
}

/*
 * PARNAME: extract name from name=value string.
 * W A R N I N G :  returns ptr to static storage, that may also be empty
 */

local string parname(string arg)
{
    permanent char namebuf[64];
    char *ap, *np;

    ap = (char *) arg;
    while (*ap == ' ')          /* skip initial whitespace */
        ap++;
    np = &namebuf[0];               /* point to spot where to copy to */
    while ((*np = *ap) != 0) {       /* copy until '=' or ' ' */
#if 0
        if (*np == '=' || *np == ' ') {
#else
        if (*np == '=') {
#endif
            *np = 0;
            return ((string) namebuf);
        }
        np++;
        ap++;
    }
#if 1
    namebuf[0] = 0;
    return ((string) namebuf);  /* return pointer to null string */
#else
    return NULL;                /* return NULL (nothing found) */
#endif    
}

/*
 * PARVALUE: extract value from name=value string, skipping initial whitespace
 *           if no value, pointer to 0 (last element) is returned!
 *           if no key, pointer to (local copy of ) first element is returned
 *  ???      BUG:  when HELPVEC is not defined, helpvec also returned
 *           Note: returns unsafe pointer into the input string
 *           
 */

local string parvalue(string arg)
{
#if 0
    permanent char *valbuf = NULL;      /* dynamic array, with (re)allocate */
    permanent int size = 0;
#else
    permanent char valbuf[MAXBUF];         /* static array: dangerous */
#endif

    char *ap;

    ap = (char *) arg;
    while (*ap) {
        if (*ap++ == '=') {
            while (*ap && *ap == ' ')      /* skip whitespace after '=' */
                ap++;
            strncpy(valbuf,ap,255);
            valbuf[255] = 0;
            ap = valbuf;
            while (*ap) {
                if (*ap == '\n') {
                    *ap = 0;
                    return valbuf;          /* return patched "value\nhelp" */
                }
                ap++;
            }
            return valbuf;                  /* return unmodified "value" */
        }
    }
#if 1        
    return ((string) ap);	    /* return pointer to null string */
#else
    return NULL;                    /* return nothing found */
#endif    
}
/* 
 * PARHELP: extract help from a defv[] "keyword=value\n help" string
 *          If no help part (or \n), returns zero string.
 *          Note: returns pointer into original string, assumed
 *          to be R/O (might be in text space of code)
 */
local string parhelp(string arg)
{
    char *cp = (char *) arg;

    while (*cp  && *cp != '\n')      /* scan until EOS or BOH */
        cp++;
    if(*cp == '\n')
        cp++;
    while (*cp && (*cp == ' ' || *cp == '\t'))   /* skip white space > BOH */
        cp++;
    return cp;
}

/*
 * FINDKEY:  scan valid keywords and return index (>=0) for a keyname
 *           return -1 if none found
 *	     Optionally match can be done in minimum match mode
 *
 *  to be deprecated in favor of 'findakey' (see below)
 */

local int findkey(string name)
{
    int i, j, l, count, last;

    if (nkeys<=0) return -1;                /* quick: no keywords at all */
    for (i = 0; i < nkeys; i++)             /* try exact match */
        if (streq(keys[i].key,name)) return i;

#if defined(MINMATCH)
    l = strlen(name);                       /* try minimum match */
    count = 0;                              /* count # matches on */
    for (i=1; i<nkeys; i++) {               /* any of the program keys */
        if (strncmp(keys[i].key, name,(unsigned)l)==0 &&
	    keys[i].indexed == -2) {
            last = i;
            count++;
        }
    }
    if (count==1) {
        warning("Resolving partially matched keyword %s= into %s=",
                name, keys[last].key);
        return last;
    } else if (count > 1) {
        dprintf(0,"Minimum match failed for \"%s\", found: ",name);
        for (j=0; j<nkeys; j++)
            if (strncmp(keys[j].key,name,(unsigned)l)==0)
                dprintf(0,"%s ",keys[j].key);
        dprintf(0,"\n");
        error("Ambiguous keyword %s=",name);
    }
#endif

    return -1;          /* if all else fails: return -1 */
}

/* findakey:  like findkey, but now returns pointer to the keyword
 *            will also find indexed keywords if so desired
 */

local keyword *findakey(string name)
{
  int i, j, l, count, last;

  if (nkeys<=0) return NULL;              /* quick: no keywords at all yet */

  for (i = 0; i < nkeys; i++)             /* first try an exact match */
    if (streq(keys[i].key,name)) return &keys[i];

#ifdef INDEXED
  {
  char *cp, *cp1, key[MAXKEYLEN];
  int n = strlen(name);
  int idx;
  keyword *kw;

  strcpy(key,name);
  dprintf(2,"findakey: checking indexing on %s\n",key);

  /* split basename from index number by working from the back */
  /* should go in a private function , a.k.a. refactoring :-)  */

  cp = &key[n-1];
  while (isdigit(*cp))
    cp--;
  if (*(cp+1)) {        /* ok, this keyword name must be indexed now */
    idx = atoi(cp+1);
    cp++;
    *cp = 0;
    strcat(key,"#");
    n = findkey(key);
    if (n < 0)  return NULL;
    // error("findakey:  #=%d findkey(%s) -> %d NOT INDEXED",idx,key,n);
    kw = &keys[n];
    if (kw->indexed < -1) 
      error("findakey(%s): not an indexed keyword, %s: %d n=%d", 
	    name,kw->key,kw->indexed,n);
    while (kw->next) {
      kw = kw->next;
      if (kw->indexed == idx+1) {
#if defined(MACROREAD)
	cp = kw->val;
	if (*cp == '@') {
	  cp1 = kw->val;
	  kw->val = get_macro(cp);
	  free(cp1);
	}
#endif
	return kw;
      }
    }
    return NULL;
  } else {              /* keyword was definitely not indexed */
    return NULL;
  }

  n = indexparam(key,idx);
  dprintf(2,"Re-Checking indexparam(%s,%d) -> %d\n",key,idx,n);
  //return n>0;
  return NULL;
  }
#endif

#if defined(MINMATCH)
  l = strlen(name);                       /* try minimum match */
  count = 0;                              /* count # matches on */
  for (i=1; i<nkeys; i++) {               /* any of the program keys */
    if (strncmp(keys[i].key, name,(unsigned)l)==0 &&
	keys[i].indexed == -2) {
      last = i;
      count++;
    }
  }
  if (count==1) {
    warning("Resolving partially matched keyword %s= into %s=",
	    name, keys[last].key);
    return &keys[last];
  } else if (count > 1) {
    dprintf(0,"Minimum match failed for \"%s\", found: ",name);
    for (j=0; j<nkeys; j++)
      if (strncmp(keys[j].key,name,(unsigned)l)==0)
	dprintf(0,"%s ",keys[j].key);
    dprintf(0,"\n");
    error("Ambiguous keyword %s=",name);
  }
#endif

  return NULL;    // will it ever get here....  it should not !
}

/* 
 * this routine is only used to enter a new indexed keyword, if
 * it turns out it is a indexed keyword. Returns 0 if the keyword
 * is not indexed.
 *
 * If the keyword is indexed, it returns the index into defv[] array
 * (1 being the first program keyword, argv0 was already reserved for 0)
 */

local int set_indexed(string name, int *idx)
{
  char *cp, key[MAXKEYLEN+1], keyidx[MAXKEYLEN];
  int j;

  dprintf(1,"set_indexed(%s)\n",name);
  *idx = -1;
  strcpy(key,name);                 /* local copy of the keyword to mess with */
  cp = &key[strlen(key)-1];
  if (!isdigit(*cp)) return 0;      /* definitely not indexed */
  while (isdigit(*cp))              /* work backwords through all digits */
    *cp--;
  cp++;                           
  strcpy(keyidx,cp);
  *idx = atoi(keyidx);
  *cp = 0;
  strcat(key,"#");
  j = findkey(key);
  dprintf(1,"set_indexed: now at %c, base=%s j=%d idx=%s -> %d\n",*cp,key,j,keyidx,*idx);
  if (j>0) return j;
  return 0;
}

/*
 * addindexed: 
 */
local int addindexed(int i, string keyval, int idx)
{
  keyword *kw;

  kw = &keys[i];       /* start at the indexed (base) keyword */
  while (kw->next) {        /* link through all indexed ones */
    dprintf(1,"Link List Skipping %s\n",kw->key);
    kw = kw->next;
    if (kw->indexed == idx+1) {
      if (kw->count) {
        kw->count++;
	dprintf(1,"Skipping entering duplicated indexed keyword %s",keyval);
	return 0;
      }
      error("Duplicated indexed keyword %s",keyval);
    }
  }
  kw->next = (keyword *) allocate(sizeof(keyword));
  kw = kw->next;
  kw->keyval = scopy(keyval);
  kw->key = scopy(parname(keyval));
  kw->val = scopy(parvalue(keyval));
  kw->help = 0;
  kw->count = 1;
  kw->upd = 0;
  kw->indexed = idx+1;
  kw->next = NULL;
  dprintf(1,"Link List new keyword %s, idx=%d\n",keyval,idx);
  return 0;
}

#if defined(INTERACT)

/* READKEYS, WRITEKEYS:
 * 
 *  These local routines read and write the keyword files.
 *
 *  Some peculiar notion to the readkeys routine is that
 *  the first time around the keyword file does not have to exist
 *
 *  The writekeys routine always writes out the internal version
 *  ID, not the external. This may result in excessive warnings
 *  when lot's of keyroutine I/O is done with an old key file
 */

local void writekeys(string mesg)
{
    FILE *keyfile;
    int i;
    keyword *kw;

    keyfile = fopen(key_filename,"w");           /* open to (over) write */
    if (keyfile==NULL)
        error("%s: Cannot write to keyfile \"%s\", - no write permission?",
                                                 mesg, key_filename);
    dprintf(5,"Writing to keyfile %s\n",key_filename);
    fprintf(keyfile,
        "# Program: %s\n",progname);
    fprintf(keyfile,
        "# keyword file written by nemo (help level=%d)\n",help_level);
    for (i=1; i<nkeys; i++) {
        if (streq(keys[i].key,"VERSION"))
	  fprintf(keyfile,"VERSION=%s\n",version_i ? version_i : "*");
        else {
	  if (keys[i].next) {
	    dprintf(1,"writing indexed keys");
	    kw = &keys[i];
	    while (kw->next) {
	      kw = kw->next;
	      fprintf(keyfile,"%s=%s\n",kw->key,kw->val);
	    }
	  } else if (keys[i].indexed == -2) {
            fprintf(keyfile,"%s=%s\n",keys[i].key,keys[i].val);
	  } 
	}
    }
    fprintf(keyfile,
        "#### end of keywords - Save file and exit editor to execute program\n");
    fprintf(keyfile,
        "#### To get more help on keywords use help=h command line option\n");
    fclose(keyfile);
} /* writekeys */


/*
 * READKEYS: read defaults from a keyword file. The first time around the
 *           keyfile doesn't have to exist, if it does, already assigned
 *           defaults (from the command line) are not overwritten
 *           The second time around (assumed after keyfile written by
 *           writekeys) all values are re-assigned, if different.
 *
 */

local void readkeys(string mesg, bool first)
{
    FILE *keyfile;
    keyword *kw;
    int i, idx;

    keyfile = fopen(key_filename,"r");
    if (keyfile==NULL && !first)
        error("%s Cannot read keyfile \"%s\" ",mesg,key_filename);
    if (keyfile==NULL) return;
    dprintf(5,"Reading from keyfile %s\n",key_filename);
    while (fgets(keybuf,MAXBUF,keyfile) != NULL) {
        if (keybuf[0]=='#' || keybuf[0]=='\n' || keybuf[0]==' ')  /* skip */
            continue;
        if (keybuf[strlen(keybuf)-1] != '\n')
            warning("readkeys: reading incomplete lines from %s",key_filename);
        keybuf[strlen(keybuf)-1] = '\0';          /* get rid of newline */
        if(streq(parname(keybuf),"VERSION")) 
	  if (version_i && !streq(version_i,parvalue(keybuf)))
                warning("readkeys: internal[%s] and external[%s] VERSION differ",
                    version_i,parvalue(keybuf));
        /* should process system keywords here... */

        i = findkey(parname(keybuf));             /* get index */
        if (i<=0) {                                  /* illegal keyword ? */
	  i = set_indexed(parname(keybuf),&idx);
	  if (i)
	      addindexed(i,keybuf,idx);
	  continue;
        }
        if (keys[i].count && first) {
            continue;    /* ignore: it was already in argv */
        } else {
            if (!streq(keys[i].val,parvalue(keybuf))) {     /* if different */
                free(keys[i].val);                          /* free old */
                keys[i].val = scopy(parvalue(keybuf));      /* and patch in */
            }
        }
    } 
    fclose(keyfile);
    return;
} /* readkeys */
#else
/*          dummy routines */
local void writekeys(string mesg) {}
local void readkeys(string mesg, bool first) {}
#endif

#if defined(INTERRUPT)
        /* should not be called if I/O redirection is performed */
local void review()
{
        char cmd[80], *key, *cp, prompt;
        int  i, pid;

        if (!initparam_done) {
            fprintf(stderr,"Cannot interrupt yet, let initparam() finish\n");
            return;
        }
        if (review_flag == 0) return;
                
        fflush(stdin);
        fflush(stdout);
        fprintf(stderr,"REVIEW called, enter your command: \n");
        fflush(stderr);

        prompt = '>';

        for (;;) {                              /* requesting loop */
           sprintf(cmd,"(%d)", getpid());       /* get process id */
           if (!debug_level)
              cmd[0] = '\0';
           fprintf(stderr,"REVIEW|%s%s%c ",getargv0(),cmd,prompt);
           fflush(stderr);
#if 0
           gets(cmd);
#else
           fgets(stdin,80,cmd);
#endif
           if (strcmp(cmd,"end")==0 || strcmp(cmd,"quit")==0) { /* all done */
	     local_exit(0);                                    /* plain exit */
	     /*NOTREACHED*/
           }
           if (strcmp(cmd,"stop")==0) {                 
              stop(0);                                   /* graceful exit */
           }else if (strcmp(cmd,"continue")==0 || strcmp(cmd,"go")==0) {
              fflush(stdin);
              fflush(stdout);
              break;
           }else if (strncmp(cmd,"task ",5)==0) {
              cp = &cmd[5];
              pid = fork();
              if (pid<0) {
                 dprintf(0,"Unable to fork: %s\n",cmd);
              } else if (pid==0) {      /* child - will run the job */
                 /* execvp(path,args); */    /* do it */
                 dprintf(0,"[STILL TO DO] Unable to execute: %s\n",cmd);
                 local_exit(1);   /* if child didn't work out, kill this thread */
                 /*NOTREACHED*/
              } 
              /* now in parent : can just live on...and on...and on...*/
           }else if (strncmp(cmd,"set ",4)==0) {        /* set new value for key */
              cp = &cmd[4];
              key = parname(cp);
              i = findkey(key);
              if (i==0)         
                  fprintf(stderr,"Cannot change name of program....\n");
              else if (i>0) {   /* regular keyword */
                  if (strlen(bindvec[i]) < strlen(cp))
                       bindvec[i] = (char *) allocate(strlen(cp));
                  strcpy (bindvec[i],cp);
              } else {          /* unknown keyword, check if system */
                  i = check_syskey(cp,"debug",&debug_level) |
                      check_syskey(cp,"help",&help_level)   |
                      check_syskey(cp,"yapp",&yapp_dev)     |
                      check_syskey(cp,"error",&error_level); 
                  if (!i)
                      fprintf(stderr,"Syntax error (%s)\n",cmd);
              }
           } else if (strncmp(cmd,"show ",5)==0) {      /* show a keyword */
              key = &cmd[5];
              i = findkey(key);
              if (i>0)
                  fprintf(stderr,"%s\n",bindvec[i]);
              else
                  fprintf(stderr,"keyword %s unknown\n",key);
           } else if (strcmp(cmd,"keys")==0) {          /* show all keywords */
              for (i=0; bindvec[i] != NULL; i++)
                 fprintf(stderr,"%s\n",bindvec[i]);
           } else if (strcmp(cmd,"syskeys")==0) {   /* show system keywords */
              fprintf(stderr,"debug=%d\n",debug_level);
              fprintf(stderr,"error=%d \n",error_level);
              fprintf(stderr,"help=%d\n",help_level);
              fprintf(stderr,"yapp=%d (%s)\n",yapp_dev,yapp_string);
           } else if (strcmp(cmd,"version")==0) {
              fprintf(stderr,"version id = %s\n",GETPARAM_VERSION_ID);
           } else if (strcmp(cmd,"time")==0) {
              fprintf(stderr,"cputime = %g\n",cputime());
           } else if (cmd[0]=='!' && cmd[1] != '\0') {  /* run a command */
              if (system(&cmd[1])<0)
                  fprintf(stderr,"Couldn't spawn another task\n");
           } else if (cmd[0]=='?' || strncmp(cmd,"help",4)==0) {
              if (strlen(cmd) < 6) {
                fprintf(stderr,"help,?             help (this list)\n");
                fprintf(stderr,"continue,go        back to program\n");
                fprintf(stderr,"end,quit           exit() to shell\n");
                fprintf(stderr,"stop               graceful finiparam() exit\n");
                fprintf(stderr,"task name par(s)   start up new task\n");
                fprintf(stderr,"set <key>=<val>    enter new val for key\n");
                fprintf(stderr,"show <key>         show val of key\n");
                fprintf(stderr,"keys,set           show values of all keys\n");
                fprintf(stderr,"syskeys,sets       show values of all syskeys\n");
                fprintf(stderr,"set <syskey>       show value of system keyword\n");
                fprintf(stderr,"set <syskey>=<val> set new val to system keyword\n");
                fprintf(stderr,"time               show cputime used\n");
                fprintf(stderr,"!cmd               run a shell command\n");
                fprintf(stderr,"<HOTKEY>           cause coredump anyhow\n");
                fprintf(stderr,"version            display version of INITPARAM\n");
              } else {
                review_help(&cmd[5]);
              }
           } else if (cmd[0]=='\0')
              continue;
           else
              fprintf(stderr,"Syntax error, try 'help'\n");
        } /* for(;;) */
}

local int check_syskey(char *cmd,char *syskey,int *value)
{
    int len;

    len=strlen(syskey);
    if (strncmp(cmd,syskey,len))
        return(0);
    if (cmd[len]=='\0') {         /* display system keyword */
        fprintf(stderr,"%s=%d\n",cmd,*value);
        return(1);
    } else if (cmd[len]=='=') {     /* set system keyword */
        if (cmd[len+1]=='\0') {
            fprintf(stderr,"System keyword must be given a value\n");
            return(1);
        }
        *value = atoi(&cmd[len+1]);
        return(1);
    } else                          /* error or unknow (system) keyword */
        return(0);
}

local void review_help(char *cmd)
{
    char command[MAXBUF], *pager;
    
    pager = getenv("PAGER");

    sprintf(command,"%s %s/%s.doc",
            (pager == NULL ? "more" : pager),
            getenv("NEMODOC"),cmd);
    printf("%s\n",command);
    system(command);
}
#endif  /* INTERRUPT */

#if defined(sun) && defined(USE_IEEE)

local void ieee_nemo(
    int sig,                    /* sig == SIGFPE always */
    int code,
    struct sigcontext *scp,	/* solaris warning: dubious tag declaration */
    char *addr
) {
 /* warning("ieee exception code 0x%x occurred at pc 0X%X", code, scp->sc_pc); */
    if (debug_level>0) abort();
}
                                                   
#endif

#if defined(TCL)
/*
 * TCL: if TCL is installed, it builds an extra comand processing layer
 *      into NEMO after the commandline and keyword files have been read,
 *      but before the program goes off and does it's thing (nemo_main).
 *      All program keywords have now become TCL variables, and may be
 *      modified in the usual TCL way (with the set command), as well as
 *      some system variables (but on a lower level)
 *
 * Commands:
 *      inp         - shows all keywords and their values
 *      go          - done, return to NEMO program
 *	keys        - return all keywords as a list
 *      set key val - change
 *      help        - rudimentary help
 */

int cmdInp(ClientData cd, Tcl_Interp *ip, int argc, char *argv[])
{
    int i;

    for (i=1; i<nkeys; i++)
        printf("%-16s= %s\n",keys[i].key,keys[i].val);
    return TCL_OK;
}

int cmdKeys(ClientData cd, Tcl_Interp *ip, int argc, char *argv[])
{
    int i;

    for (i=1; i<nkeys; i++) {
	Tcl_AppendElement(ip,keys[i].key);
    }
    return TCL_OK;
}

int cmdCd(ClientData cd, Tcl_Interp *ip, int argc, char *argv[])
{
    int i;

    if (argc != 2) {
	sprintf(ip->result, "wrong # args: should be \"%s directory\"",argv[0]);
        return TCL_ERROR;
    }
    i = chdir(argv[1]);
    if (i==0)
      return TCL_OK;
    else {
      sprintf(ip->result, "%s: No such directory",argv[1]);
      return TCL_ERROR;
    }
}

/* 
 *  load all commands from a directory 'dir' as new tcl commands.
 *  this gives some possibility to retrieve command information
 *  back from command. Of course they have to know about it..
 */

inc cmdLoad(ClientData cd, Tcl_Interp *ip, int argc, char *argv[])
{
    FILE *fp;       /* pipe */
    char line[64];
    int k;

    if (argc != 2) {
	sprintf(ip->result, "wrong # args: should be \"%s directory\"",argv[0]);
        return TCL_ERROR;
    }
    
    sprintf(line,"/bin/ls %s",argv[1]);
    fp = popen(line,"r");
    if (fp==NULL) {
        sprintf(ip->result,"Could not find directory %s",argv[1]);
        return TCL_ERROR;
    }
    k=0;
    while (fgets(line,64,fp) != NULL) {
        line[strlen(line)-1] = '\0';
        dprintf(1,"Loading command %s\n",line);
        Tcl_CreateCommand(tcl_interp, line, cmdTest, (ClientData) line, dP);
        k++;
    }
    sprintf(ip->result,"Loaded %d commands from %s.\n",k,argv[1]);
    return TCL_OK;
}

int cmdTest(ClientData cd, Tcl_Interp *ip, int argc, char *argv[])
{
    char *name = (char *) cd;
    int i, status, pid;
    
    printf("CmdTest: %s\n",argv[0]);
    printf("  Command executed as a special TCL command:");
    if ((pid=fork()) < 0) {				/* parent in error */
    	sprintf(ip->result,"Failed to fork %s",argv[0]);
    	return TCL_ERROR;
    } else if (pid==0) {				/* child */
    	execvp(argv[0],argv);
    	sprintf(ip->result,"Failed to exec process %s",argv[0]);
    	local_exit(1);
    } else {					/* parent */
    	while ( i != wait(&status))		/* wait for child to finish */
    		;
    }

    return TCL_OK;
}

int cmdGo(ClientData cd, Tcl_Interp *ip, int argc, char *argv[])
{
    tcl_flag = 0;       /* flag command parser to exit */
    return TCL_OK;
}

int cmdHelp(ClientData cd, Tcl_Interp *ip, int argc, char *argv[])
{
    printf("Try the TCL command 'info' to get more help\n");
    printf("Many commands accept the ? as argument to get help\n");
    printf("Some import commands though:\n");
    printf("  inp           overview of the current keywords\n");
    printf("  set KEY VAL   reset KEY keyword (this is TCL syntax) with VAL\n");
    printf("  go            run the program\n");
    printf("  exit          exit from the program, don't run it\n");
    printf("  load DIR      load all files from DIR as commands\n");
    printf("  cd DIR        change directory to DIR\n");
    printf("  echo ...      echo argument list\n");
    return TCL_OK;
}

int cmdEcho(ClientData cd, Tcl_Interp *ip, int argc, char *argv[])
{
    int i;

    for (i = 1; ; i++) {
	if (argv[i] == NULL) {
	    if (i != argc) {
		echoError:
		sprintf(ip->result,
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

run_tcl()
{
  int i, result;
  char buffer[MAXBUF], *cmd;

  dprintf(0,"[Running TCL]\n");

  tcl_interp = Tcl_CreateInterp();

  for (i=1; i<nkeys; i++) {             /* program keywords */
    Tcl_SetVar(tcl_interp,keys[i].key,keys[i].val,0);
    Tcl_LinkVar(tcl_interp,keys[i].key,(char *)&keys[i].val,TCL_LINK_STRING);
  }
  Tcl_LinkVar(tcl_interp,"error",(char *) &error_level, TCL_LINK_INT);
  Tcl_LinkVar(tcl_interp,"debug",(char *) &debug_level, TCL_LINK_INT);

  Tcl_CreateCommand(tcl_interp,"inp",cmdInp, (ClientData)NULL, dP);
  Tcl_CreateCommand(tcl_interp,"keys",cmdKeys, (ClientData)NULL, dP);
  Tcl_CreateCommand(tcl_interp,"cd",cmdCd, (ClientData)NULL, dP);
  Tcl_CreateCommand(tcl_interp,"load",cmdLoad, (ClientData)NULL, dP);
  Tcl_CreateCommand(tcl_interp,"go",cmdGo, (ClientData)NULL, dP);
  Tcl_CreateCommand(tcl_interp,"help",cmdHelp, (ClientData)NULL, dP);
  Tcl_CreateCommand(tcl_interp,"echo",cmdEcho, (ClientData)NULL, dP);
  
  result = Tcl_Eval(tcl_interp,"source ~/.tclrc");
  if (*tcl_interp->result != 0)
    printf("TCL_STARTUP: %s\n",tcl_interp->result);
  load_environ();    

  for (;;) {                  /* loop reading interactive commands */
    sprintf(buffer,"TCL: %s>",progname);    /* make prompt */
    myreadline(buffer,MAXBUF);              /* get command from user */
    if (buffer[0]=='!') {
    	system(&buffer[1]);
    	continue;
    }
    result = Tcl_RecordAndEval(tcl_interp, buffer, 0);
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
    if (!tcl_flag) break;
  }
  Tcl_DeleteInterp(tcl_interp);     /* free up all that memory */
}

/* a dummy routine for TCL to stub procedure deletions */
dP(char *clientData)
{
}

local void myreadline(char *buffer, int buflen)
{
#if defined(READLINE)
    char *s;

    for (;;) {                     /* Read until something is present. */
      if ((s = readline(buffer)) != NULL) {
        stripwhite(s);      /* Get rid of leading and trailing blanks. */ 
        if (*s != 0) {
           strcpy(buffer, s);       /* stuff it back into buffer */
           free(s);
           break;
        }
      }
    }
    if (buffer[strlen(buffer)-1] == '\n')     /* Need to overwrite '\n'? */
      buffer[strlen(buffer)-1] = 0;               /* ...then do it. */
    add_history(buffer);
#else    
    printf("%s",buffer);
    fflush(stdout);
    fgets(stdin,buflen,buffer);
#endif
}

local void ini_readline()
{
#if defined(READLINE)
    read_history(".nemohis");
#endif
}

local void end_readline()
{
#if defined(READLINE)
    write_history(".nemohis");
#endif
}

local void stripwhite(char *s)
{
    register int i, j;

    for (i = 0; isspace(s[i]); i++)
      /* NULL */ ;

    if (i > 0) /* Do not use strcpy() here! */
      for (j = 0; s[j] = s[i]; j++)
        i++;

    for (i = strlen(s) - 1; ((i > 0) && (isspace(s[i]))); i--)
      s[i] = 0;

}



load_environ()
{
    char *cp, *cp1, line[512];
    int i, k, kbad, len;
    
    k = kbad = 0;
    for (i=0, cp=environ[0]; cp != NULL; cp=environ[++i]) {
    	len = strlen(cp);
        if (len > 512-1) {
            cp1 = strchr(cp,'=');
            if (cp1==NULL) continue;	/* should never happen */
            *cp1 = 0;
            dprintf(0,"%s not saved -- too long (%d)\n",cp,len);
            kbad++;
            continue;
        }
        strcpy(line,cp);
        cp = strchr(line,'=');
        if (cp==NULL) {
            printf("### %s: not a proper environment variable\n");
            continue;
        }
        *cp = '\0';  /* patch it */
        cp++;
        Tcl_SetVar(tcl_interp,line,cp,1);
        k++;
    }
    printf("Loaded %d environment variables, discarded %d (too long).\n",
		k,kbad);
}

#endif

/*
 * SET_XXX: 	set system variables
 */
 
local void set_help(string arg)
{
    char *cp;
    
    help_string = scopy(arg);
    if ((cp = strpbrk(help_string,"0123456789"))!=NULL)  /* isnum? */
        help_level = atoi(cp);  /* if so, change help_level */
}

local void set_outkeys(string arg)
{
    char *cp;
#ifdef OUTKEYS    
    warning("No output key processing done yet (%s)",arg);
#else
    error("outkeys= not supported, or recompile with -DOUTKEYS");
#endif

}

local void set_debug(string arg)
{
    debug_level = atoi(arg);
}

local void set_yapp(string arg)
{
    yapp_string = scopy(arg);
    yapp_dev = atoi(yapp_string);
}

local void set_review(string arg)
{
    review_flag = atoi(arg);
}

local void set_tcl(string arg)
{
    tcl_flag = atoi(arg);
}

local void set_error(string arg)
{
    error_level = atoi(arg);
}

local void local_error(string msg)
{
    fprintf(stderr,"### Fatal error in getparam.c: %s\n", msg);
    exit(-1);
}
local void local_exit(int level)
{
    finiparam();
    exit(level);
}

#if defined(TESTBED)

string defv[] = {
    "flag1=???\n        This needs a value. Perhaps even more",
    "flag2=false\n      A false flag2",
    "infile=r000.snap\n just a name, but what's is a name",
    "tstop=10.0\n       just a number (you may try parsing an expression)",
    "foobar=12345\n     just some digits",
    "donot=read\n       check if complains about never read params",
    "prompt=t\n		Checking interactive prompting routines?",
#ifdef INDEXED
    "naxis#=\n          testing an indexed keyword",
#endif
    "VERSION=1.4\n      16-sep-01 PJT",
    NULL,
};

string usage = "A testbed for getparam.c";

void check(string);

void nemo_main(void)
{
    bool flag2, prompt;
    double tstop;
    int i, n, my_argc;
    string *my_argv;

    my_argv = getargv(&my_argc);
    if (my_argc) {
      printf("===> unparsed arguments <===\n");
      for (i=0; i<my_argc; i++)
	printf("argv[%d] = %s\n",i,my_argv[i]);
    }

    printf("===> program %s: <===\n", getargv0());
    check("argv0");
    check("flag1");
    check("flag2");
    check("infile");
    check("tstop");
    check("foobar");
    check("prompt");

    dprintf(0,"debug level 0 printout\n");
    dprintf(1,"debug level 1 printout\n");
    dprintf(2,"debug level 2 printout\n");
    dprintf(9,"debug level 9 printout\n");
    dprintf(11,"debug level 11 printout\n");

#ifdef INDEXED
    n = indexparam("naxis",0);
    if (n<1) warning("No indexed keyword naxis### used");
    printf("Indexed keyword:\n");
    for (i=0; i<=n; i++)
      if (getparam_idx("naxis",i))
            printf("naxis%d=%s\n",i,getparam_idx("naxis",i));
     
#endif


    flag2 = getbparam("flag2");    
    tstop = getdparam("tstop");
    prompt = getbparam("prompt");
    printf("flag2=%d tstop=%f prompt=%d\n",flag2,tstop,prompt);
    if (!prompt) stop(0);

#if defined(PUTPARAM)
    putparam("flag1","false");      /* resetting flag1 to false */
    if(getbparam("prompt")) {
        for (;;) {              /* infinite loop */
#if defined(INTERRUPT)
            promptparam("foobar",
                "(REVIEW interruptable) Enter a new value for foobar");
            check("foobar");
#else
            promptparam("flag1","flag1 determines end of infinite loop");
            check("flag1");
#endif 
            if (getbparam("flag1"))
                break;
        }
    }
    warning("Next you should see a fatal error");
    putparam("nosuchkey","nosuchprompt");
    warning("If you read this, you bypassed the previous fatal error");
#endif

} /* nemo_main */

void check(string par)
{
    string val;

    val = getparam(par);
    printf("keyword %s has value \"%s\"\n", par, val);
}

#endif  /* TESTBED */


