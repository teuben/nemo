#define VERSION_ID "2.4b 11-nov-90"
/*= mirtool */
/*& pjt */
/*: tools */
/*+*/
/*  Suntools/view based executioner for miriad programs			*/
/*  For more information see appropriate tex file/documents		*/
/*--*/
/***********************************************************************\
*									*
*	ExecTool Version (version # see below)		                *
*									*
*	Written by B. Sutin - University of Maryland                    *
*	as part of the Mirth/Miriad package				*
*	Now being updated by P. Teuben                                  *
*                                                                       *
*	Please send comments suggestions, bug fixes, gripes,            *
*	improvements, and changes to teuben@astro.umd.edu               *
*                                                                       *
* Updates:                                                              *
*      xx-jun-89:  V2.2 compile switch -DNEMO creates a NemoTool        *
*      13-jul-89:  help button, host button, 'rsh' bugs removed         *
*      20-jul-89:  V2.3 nemo also uses .def for keyword filename        *
*      29-jul-89:  3 lines for top lines; prepare for load status and   *
*                   dynamic search directory                            *
*      29-sep-89:  more stuff						*
*  rjs 28-nov-89:  Fixed bug when creating the executable menu on Sun-4.*
*  pjt 27-may-90:  all versions include mirth.icon now                  *
*  mjs 30-jun-90:  patch ... thoroughly de-activated the "hostname"     *
*                  option because of incompatibility with NIS on OS4.1; *
*                  routed select-action to a no-op subroutine and       *
*                  stuck a "" as the menu's only entry ... search for   *
*                  string "spool" to see the pertinent subroutine.      *
*	<< All files bundled into one .c file >>			*
*  pjt 15-jul-90   Renamed file to exectool.c - made it less verbose    *
*                  rearranged sorting (bin -> sbin)                     *
       18-jul-90   too many NEMO def's -> MIR
       12-nov-90   removed stuparisms

Known bugs/deficience
  - putenv/getenv don't seem to work on Sun, hence the -b,-p,-l etc.
	don't work. Perhaps use the miriad.c solution. Aha, got
	to share files here... Perhaps put it all in one dir eh?
        Now solved by getting environment before processing cmdln flags

Work To Do:
  - see if pipe's can have different buffer size (ioctl?) to improve
    interactive programs writing their output while user is still using
    program. 
  - remote execution: fix it for Mark
  - dynamic loading of all doc's and bin's as through NEXEC and NDOC
    which also should fix the problem that the tool doesn't have to be
    restarted when new programs have been added to the BIN area
  - make log screen as wide as first window
  - fix help (middle button for keywords, right button for mirtool
  - ecount returned by cmpbin_doc() is wrong.
  + resort list of binaries for column wise alphabetise [**done**]
\***********************************************************************/

#include <suntool/sunview.h>
#include <suntool/panel.h>
#include <suntool/textsw.h>
#include <suntool/tty.h>
#include <suntool/scrollbar.h>
#include <stdio.h>
#include <sys/file.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/dir.h>
#include <string.h>
#include <sys/ioctl.h>
/* #include <sys/tty.h> */
#include <protocols/rwhod.h>

/* 	Pick a sensible name for the tool from supplied compiler flag */
#if defined(NEMO)
#define TOOL_NAME "nemotool"
#else
#define TOOL_NAME "mirtool"
#endif




/**************** begin install section ****************/
/***************** end install section *****************/

/************************************************************************
 *									*
 *	Global Variables						*
 *									*
 ************************************************************************/

Frame frame ;					/* main menu frame */
  Panel control ;				/* main control panel */
    static void quit() ;			/* quit button */
    static void run() ;				/* run button */
      Notify_value reap() ;			/* reap dead jobs */
    static void save() ;			/* save button */
    static void help() ;                        /* help button */
    static void hosttoggle();                   /* local/host toggle button */
    Panel_item dir_text ;			/* current directory */
      Panel_setting dir_notify_proc() ;		/* directory change */
    Panel_item message ;			/* messages */
    Panel_item name_text ;			/* routine name */
      Menu name_menu ;				/* menu of routine names */
      name_event_proc() ;			/* name events */
      Panel_setting name_notify_proc() ;	/* name events */
    Panel_item host_text ;			/* remote host name */
      Menu host_menu ;                          /* menu of host names */
      host_event_proc() ;                       /* host events */
      Panel_setting host_notify_proc() ;        /* hosts events */
    Panel_item log_text ;			/* logfile name */
      Menu log_menu[4] ;			/* menu of log functions */
      log_event_proc() ;			/* log events */
      Panel_setting log_notify_proc() ;		/* log events */
  Panel inputs ;				/* routine input list */
    Scrollbar input_scroll ;			/* input panel scrollbar */
    input_help_proc() ;				/* input help */
  Panel job_win ;				/* job list panel */
    int job_event_proc() ;			/* job events */
    Menu job_menu ;				/* menu of job signals */
    draw_jobs() ;				/* update job list */
Frame help_frame ;				/* help popup window */
  Textsw help_text ;				/* help text */
Frame log_frame ;				/* logfile window */
  Tty log_tty ;					/* logfile tty */

Panel_item citem5;				/* toggle local/remote */

#define NEXEC	200				/* max executable files */
#define NDOC    200                             /* max doc files */
#define NCHAR	128				/* max characters in a line */
#define NJOBS	 10				/* max concurrent jobs */
#define NARGS	100				/* max arguments on a cmd */
#define NHOST    50                             /* max remote hosts */

#define PWIDTH	20				/* window depth */

char doc_path[NCHAR] ;				/* current document path */
char bin_path[NCHAR] ;				/* current binary path */
char def_path[NCHAR] ;				/* current defaults path */

char *LOG_FILE ;				/* current LOG_DESC path */
char LOG_SAVE[NCHAR] ;				/* saved log file name */

static char *arglist[NARGS+10] ;		/* exec argument list */

char *getenv(), *malloc(), *strcpy(), *strcat() ;

struct {					/* executable kludge */
    char *files[NEXEC] ;
    } bin ;

struct {					/* vert. sorted bins */
    char *files[NEXEC] ;
    } sbin ;

struct {					/* doc kludge */
    char *files[NDOC] ;
    } doc ;

struct {
    char *names[NHOST] ;
    } hosts;

struct {					/* argument structure */
    Panel_item	item ;
    char 	name[20] ;
    char	*value ;
    } arg[NARGS] ;

struct JOB_STRUCT {				/* current job structure */
    int		njob ;
    int		pid ;
    char	name[20] ;
    char        hostname[20] ;
    char        status;             /* 'Running', 'Paused, 'Killed' */
    } jobs[NJOBS] ;

Panel_item job_message[NJOBS] ;			/* job messages */
int njob = 0 ;					/* number of last job */
char *cwd;                           /* name of current working directory */
int remote_exec = 0;            /* toggle to execute remote/local */
int debug_level = 0;		/* see -d flag */


main(ac,av)
int ac;
char *av[];
{
    int arg ;

    set_defs(ac,av);	/* install command line options and check environment */
    chk_env() ;                             /* now check if all is set */

    for( arg = 0 ; arg < NARGS + 10 ; arg++ )
	arglist[arg] = malloc(NCHAR) ;

    mk_frame() ;				/* make main window */
    mk_control() ;				/* make control panel */
    mk_inputs() ;				/* make argument panel */
    mk_jobs() ;					/* make job panel */
    mk_help() ;					/* make help window */

    window_fit(frame) ;				/* let 'er rip */
    window_main_loop(frame) ;
    }

/************************************************************************
 *									*
 *	chk_env -- check the environment variables			*
 *									*
 ************************************************************************/

int NUL_DESC ;					/* desc of /dev/null */
int LOG_DESC ;					/* desc of logfile */

char *DOC_DIR ;					/* document directory */
char *BIN_DIR ;					/* binary directory */
char *DEF_DIR ;					/* defaults directory */

#if defined(NEMO)
# define BIN_DIR_NAME  "NEMOBIN"
# define DOC_DIR_NAME  "NEMODOC"
# define DEF_DIR_NAME  "NEMODEF"
# define LOG_FILE_NAME "NEMOLOG"
#else
# define BIN_DIR_NAME  "MIRBIN"
# define DOC_DIR_NAME  "MIRPDOC"
# define DEF_DIR_NAME  "MIRDEF"
# define LOG_FILE_NAME "MIRLOG"
#endif

chk_env()
{
    NUL_DESC = open( "/dev/null", 0 ) ;
    LOG_DESC = open( "/dev/null", 0 ) ;

    if( !BIN_DIR ) {
	(void) fprintf( stderr, "%s not set\n", BIN_DIR_NAME );
	exit(0) ;
    } else if (debug_level)
        printf("%s = %s\n",BIN_DIR_NAME,BIN_DIR);

    if( !DOC_DIR ) {
	(void) fprintf( stderr, "%s not set\n", DOC_DIR_NAME ) ;
	exit(0) ;
    } else if (debug_level)
        printf("%s = %s\n",DOC_DIR_NAME,DOC_DIR);

    if( !DEF_DIR ) {
	(void) fprintf( stderr,
	    "%s not set -- using current directory\n", DEF_DIR_NAME ) ;
	DEF_DIR = "." ;
    } else if (debug_level)
        printf("%s = %s\n",DEF_DIR_NAME,DEF_DIR);

    if (debug_level)
        printf("%s = %s\n",LOG_FILE_NAME,LOG_FILE);

}

set_defs(ac,av)
int ac;
char *av[];
{
/*
    This gets any command line parameters
    -g      turn on debug mode
    -p      doc area
    -b      bin area
    -d      def area
    -l      log file
    -
------------------------------------------------------------------------*/
    int i;
    char *cp, newenv[NCHAR];

    BIN_DIR = getenv( BIN_DIR_NAME );
    DOC_DIR = getenv( DOC_DIR_NAME );
    DEF_DIR = getenv( DEF_DIR_NAME ) ;
    LOG_FILE = getenv( "LOG_FILE_NAME" ) ;
    i=1;
    while(i < ac) {         /* avoid using getopt here */
        cp = av[i];
        if (*cp++=='-') {
            switch(*cp) {
              case 'g':
                    debug_level = 1;
                    break;
              case 'd':                    
                    DEF_DIR = av[++i];
                    break;
              case 'b':
                    BIN_DIR = av[++i];
                    break;
              case 'p':
                    DOC_DIR = av[++i];
                    break;
              case 'l':
                    LOG_FILE = av[++i];
                    break;
              case '?':
              case 'h':
                    usage(av[0]);    /* will also exit(0) */
              default:
                    fprintf(stderr,"%s: Illegal option %s\n",av[0],cp);
            }
            i++;
        } else
            fprintf(stderr,"%s: Illegal option %s\n",av[0],--cp);
    } /* while */
}

usage(name)
char *name;
{
    fprintf(stderr,"Usage: %s [-d] \n",name);
    fprintf(stderr,"Options:\n");
    fprintf(stderr,"  -g          Turn on debugging output mode\n");
    fprintf(stderr,"  -b bindir   Actual directory for binaries\n");
    fprintf(stderr,"  -p docdir   Actual directory for doc\n");
    fprintf(stderr,"  -d defdir   Actual directory for def files\n");
    fprintf(stderr,"  -l logfile  Logfile used\n");
    fprintf(stderr,"  -?/h        This help\n");
    exit(0);
}

mk_frame()
{
    char buff[NCHAR] ;

    static short icon_image[] = {
/*
    Format_version=1, Width=64, Height=64, Depth=1, Valid_bits_per_item=16
*/
        0xFFFF,0xFFFF,0xE0FE,0x07C0,0xFFFF,0xFFE0,0x1CFE,0x07C0,
        0xFC03,0xFD81,0xDF3E,0x07C0,0xFE7D,0xF183,0x9F06,0x07E0,
        0xFE7C,0xC1C3,0x9F01,0x03F0,0xFE7C,0x01C3,0x9F80,0x83F8,
        0xFE7D,0x01C3,0x8F80,0xC1FF,0xFE7A,0x7F65,0x8980,0x60FF,
        0xFE0E,0x1965,0x8980,0x787F,0xFE63,0x1965,0x84C0,0x3C3F,
        0xFE41,0x9935,0x8CC0,0x1E0F,0xFE01,0x9939,0x8EE0,0x1F00,
        0xFE81,0x9939,0x9260,0x0F80,0xFE81,0x9911,0x9270,0x0380,
        0xFF83,0x1911,0x9178,0x01C0,0xFFFE,0x1B93,0xFE3C,0x0060,
        0xF800,0x1800,0x203C,0x0000,0xF800,0x1800,0x203E,0x000F,
        0xF000,0x1800,0x405F,0x8007,0xF000,0x1800,0x405F,0xC007,
        0xE000,0x7E00,0xE0EF,0xE003,0xE000,0x0000,0x000F,0xF803,
        0xE000,0x0000,0x0003,0xFE03,0xC000,0x0000,0x0001,0xFFE1,
        0xC000,0x0000,0x0000,0xFFFE,0xC000,0x0000,0x0000,0x3FFE,
        0xC000,0x0000,0x0000,0x0FFE,0x8000,0x0001,0x8000,0x03FF,
        0x8000,0x0001,0x8000,0x003F,0x8000,0x0000,0x0000,0x0000,
        0x8000,0x0000,0x0000,0x0000,0x8000,0x0001,0x8000,0x0000,
        0x8000,0x0001,0x8000,0x0000,0x8000,0x0000,0x0000,0x0000,
        0x8000,0x0000,0x0000,0x0000,0x8000,0x1999,0x9998,0x0000,
        0x8000,0x1999,0x9998,0x0000,0x8000,0x0000,0x0000,0x0000,
        0xC000,0x0000,0x0000,0x0001,0xC000,0x0000,0x0000,0x0001,
        0xC000,0x0000,0x0000,0x0001,0xC000,0x0000,0x0000,0x0001,
        0xE000,0x0000,0x0000,0x0003,0xE000,0x0000,0x0000,0x0003,
        0xE000,0x0000,0x0000,0x0003,0xF000,0x0000,0x0000,0x0007,
        0xF000,0x0000,0x0000,0x0007,0xF800,0x0000,0x0000,0x000F,
        0xF800,0x0000,0x0000,0x000F,0xFC00,0x0000,0x0000,0x001F,
        0xFE00,0x0000,0x0000,0x003F,0xFF00,0x0000,0x0000,0x007F,
        0xFF00,0x0400,0x8000,0x307F,0xFF86,0x8051,0xE30C,0x10FF,
        0xFFC5,0x4C68,0x8492,0x11FF,0xFFE5,0x4440,0x8492,0x13FF,
        0xFFF5,0x4440,0x8492,0x17FF,0xFFF9,0x4440,0x630C,0x1FFF,
        0xFFFE,0x0000,0x0000,0x1FFF,0xFFFF,0x8000,0x0000,0x7FFF,
        0xFFFF,0xC000,0x0001,0xFFFF,0xFFFF,0xF000,0x0007,0xFFFF,
        0xFFFF,0xFE00,0x003F,0xFFFF,0xFFFF,0xFFE0,0x03FF,0xFFFF
	} ;
    DEFINE_ICON_FROM_IMAGE( icon, icon_image ) ;

    strcpy (buff,TOOL_NAME);
    strcat (buff,VERSION_ID);
    frame = window_create(NULL, FRAME,
	FRAME_LABEL,			buff,
	WIN_ERROR_MSG,			"Cannot create window",
	FRAME_SUBWINDOWS_ADJUSTABLE,	FALSE,
	FRAME_ICON,			&icon,
	0) ;
    }

mk_control()
{
    int  ecount , dcount, hcount, i, j, k, l, nrows, ncols;
    char buff[NCHAR], **s;
    int nop();

    control = window_create( frame, PANEL,
	PANEL_SHOW_MENU,	FALSE,
	PANEL_LABEL_BOLD,	TRUE,
	PANEL_LAYOUT,		PANEL_HORIZONTAL,
	WIN_COLUMNS,		100,
	WIN_ROWS,		3,
	0 ) ;

    panel_create_item( control, PANEL_BUTTON,
	PANEL_LABEL_IMAGE,	panel_button_image( control, "Quit", 0, 0 ),
	PANEL_NOTIFY_PROC,	quit,
	0 ) ;
    panel_create_item( control, PANEL_BUTTON,
	PANEL_LABEL_IMAGE,	panel_button_image( control, "Run", 0, 0 ),
	PANEL_NOTIFY_PROC,	run,
	0 ) ;
    panel_create_item( control, PANEL_BUTTON,
	PANEL_LABEL_IMAGE,	panel_button_image( control, "Save", 0, 0 ),
	PANEL_NOTIFY_PROC,	save,
	0 ) ;
    panel_create_item( control, PANEL_BUTTON,
	PANEL_LABEL_IMAGE,	panel_button_image( control, "Help", 0, 0 ),
	PANEL_NOTIFY_PROC,	help,
	0 ) ;
    citem5 = panel_create_item( control, PANEL_BUTTON,
	PANEL_LABEL_IMAGE,	panel_button_image( control, "Local ", 0, 0 ),
	PANEL_NOTIFY_PROC,	hosttoggle,
	0 ) ;

    dir_text = panel_create_item( control, PANEL_TEXT,
	PANEL_LABEL_STRING,		"directory: ",
	PANEL_VALUE,			getwd(buff),
	PANEL_NOTIFY_STRING,		"\n\r",
	PANEL_NOTIFY_PROC,		dir_notify_proc,
	PANEL_VALUE_DISPLAY_LENGTH,	35,
	0 ) ;

    message = panel_create_item( control, PANEL_MESSAGE,
	PANEL_SHOW_ITEM,	FALSE,
	PANEL_ITEM_Y,		ATTR_ROW(0),
	0 ) ;

    name_text = panel_create_item( control, PANEL_TEXT,
	PANEL_LABEL_STRING,		"name: ",
	PANEL_NOTIFY_STRING,		"\n\r",
	PANEL_EVENT_PROC,		name_event_proc,
	PANEL_NOTIFY_PROC,		name_notify_proc,
	PANEL_VALUE_DISPLAY_LENGTH,	15,
	PANEL_ITEM_X,			ATTR_COL(0),
	PANEL_ITEM_Y,			ATTR_ROW(1),
	0 ) ;

    ecount = search( BIN_DIR, &bin, NEXEC ) ;    /* get exe's */
    dcount = search (DOC_DIR, &doc, NDOC) ;     /* get doc's */
    ecount = cmpbin_doc (&bin, ecount, &doc, dcount); /* delete non-doc exe's */
    fprintf(stderr,"A total of %d executables found in %s\n",ecount,BIN_DIR);


    ncols = 1 + ecount/20;
    nrows = (ecount-1)/ncols + 1;
    name_menu = menu_create(MENU_NCOLS, ncols, 0);

    for (l=0; l<ecount; l++) {      /* rearrange bin -> sbin for */
        i = l/nrows;                /* vertical alphabet */
        j = l - i*nrows;
        k = j*ncols+i;
        sbin.files[k] = bin.files[l];
    }

    s = sbin.files;     /* THIS IS VERT ALPHA - set to bin.files for HORI */
    for(i=0; i < ecount; i++) {
      menu_set(name_menu, MENU_STRING_ITEM, *s++, i+1, 0);
    }

    host_text = panel_create_item( control, PANEL_TEXT,
	PANEL_LABEL_STRING,		"hostname: ",
	PANEL_NOTIFY_STRING,		"\n\r",
/*
	PANEL_EVENT_PROC,		host_event_proc,
	PANEL_NOTIFY_PROC,		host_notify_proc,
*/
	PANEL_EVENT_PROC,		nop,
	PANEL_NOTIFY_PROC,		nop,
	PANEL_VALUE_DISPLAY_LENGTH,	15,
	PANEL_ITEM_Y,			ATTR_ROW(1),
	0 ) ;

    hcount = hsearch( "/usr/spool/rwho", &hosts, NHOST);
    host_menu = menu_create(
        MENU_NCOLS,             1 + hcount / 20,
/*
        MENU_STRINGS,           hosts,
*/
	MENU_STRINGS,		"",0,
        0 ) ;


    log_text = panel_create_item( control, PANEL_TEXT,
	PANEL_LABEL_STRING,		"logfile: ",
	PANEL_NOTIFY_STRING,		"\n\r",
	PANEL_EVENT_PROC,		log_event_proc,
	PANEL_NOTIFY_PROC,		log_notify_proc,
	PANEL_VALUE_DISPLAY_LENGTH,	35,
	PANEL_ITEM_Y,			ATTR_ROW(1),
	PANEL_VALUE,			LOG_FILE,
	0 ) ;

    log_menu[0] = menu_create( MENU_STRINGS,
	"show",   "page", "erase", "print", 0, 0 ) ;
    log_menu[1] = menu_create( MENU_STRINGS,
	"hide",   "page", "erase", "print", 0, 0 ) ;
    log_menu[2] = menu_create( MENU_STRINGS,
	"show", "nopage", "erase", "print", 0, 0 ) ;
    log_menu[3] = menu_create( MENU_STRINGS,
	"hide", "nopage", "erase", "print", 0, 0 ) ;
    }

static void hosttoggle()
{
    char buff[10];

    remote_exec = !remote_exec;         /* toggle */
    strcpy(buff, (remote_exec ? "remote" : "local") );
    panel_set (citem5,
        PANEL_LABEL_IMAGE, panel_button_image(control,buff, 0, 0 ),
        0 ) ;
    
}
/*
 * This function needs to be fixed - first of all NEMO vs. MIR vs.
 *			anything else -
 *			Second the RIGHT help button doesn't work
 *			reliably : There. It's not called anymore
 */
static void help()
{
    int spot;
    char buff[NCHAR];
    Textsw_index first, last;
    Panel_item   item;
    Event *event;
    
    if (debug_level) {
        printf("Help pushed\n");  fflush(stdout);
    }
    
    strcpy (doc_path, DOC_DIR);
    strcat (doc_path, "/mirtool.doc");

    if (debug_level) {
    printf("helpfile = %s\n",doc_path); fflush(stdout);
    }
    window_set( help_text, TEXTSW_FILE, doc_path, 0) ;
    first = (Textsw_index) 0;

    strcpy (buff,"suntools");
    spot = textsw_find_bytes(
         help_text, &first, &last, buff, strlen(buff), 0);

    window_set( help_text, TEXTSW_FIRST, &first, 0);
    window_set( help_text, WIN_SHOW, TRUE, 0) ;
    scrollbar_paint_bubble(window_get( help_text, WIN_VERTICAL_SCROLLBAR ) ) ;
#if 0
    panel_default_handle_event (item, event );    
#endif
}

static void quit()
{
    int len ;

    window_destroy( frame ) ;

    len = lseek( LOG_DESC, 0, 2 ) ;
    if( !len )
	unlink( LOG_SAVE ) ;
    }

/************************************************************************
 *									*
 *	Run -- the big routine						*
 *									*
 ************************************************************************/

static void run()
{
    int next, narg, hoobie ;
    char *value, *host, *name, *argsave, *p, *path, *hostname;
    char buff[1000] ;

    name = panel_get_value( name_text ) ;
    if (debug_level) {
        printf ("name_panel = (%s)\n",name);
    }
    if( *name == '\0' )         
	return ;                /* no name: nothing to run */

    save() ;			/* do we really want this? */

/*
 *  if remote, add "rsh" <host> to the argument list
 *  of the form:
 *          rsh <host> -n (cd <pwd> ;<oldcmd> &)
 *
 */

#define WHITE(x)	( (x == ' ') || (x == '\t') )

    host = panel_get_value( host_text ) ;   /* get hostname */
    hostname = host;            /* save name for later */

    if (debug_level) {    
        printf("remote execution: hostpanel =(%s)\n",host);
    }
    while( WHITE( *host ) )             /* skip white space */
	host++ ;

    narg = 0 ;
    if (remote_exec ) {
        if (*host == NULL)
            return;             /* no hostname given */

        if (debug_level) {        
            printf("remote execution: hostname = (%s)\n",host);
        }
	path = "rsh" ;

	strcpy( arglist[narg++], "rsh" ) ;

	while( *host ) {             /* loop (?) to copy name (and options) */
	    p = arglist[narg] ;              /* start copying here */
	    while( !WHITE(*host) && *host )     /* copy until first white */
		*p++ = *host++ ;
	    *p = '\0' ;
	    narg++ ;                    /* next argument */
	    while( WHITE( *host ) )     /* skip white stuff first */
		host++ ;
	}
	cwd =  panel_get_value( dir_text );
        if (debug_level) {        
            printf ("Remote executing from directory %s\n",cwd);
        }
	strcpy( arglist[narg++], "-n" ) ;
        strcpy( arglist[narg++], "(cd" ) ;
        strcpy( arglist[narg++], cwd  ) ;
        strcpy( arglist[narg++], ";") ;
    } else {
        *hostname = '\0';
	path = bin_path ;
    }

/*
 * set up an argv[] list of the form:
 *
 *	argv[0]:	name
 *	argv[1]:	arg1=value1
 *	argv[2]:	arg2=value2
 *	  ...
 *	argv[n]:	NULL
 *
 */

    hoobie = time(0) ;
    (void) strcpy( buff, asctime( localtime( &hoobie ) ) ) ;
    (void) strcat( buff, "Executing " ) ;
    (void) strcat( buff, name ) ;
    (void) strcat( buff, "\n" ) ;
    arglist[narg++] = name ;
    for( next = 0 ; next < NARGS ; next++ ) {
	value = panel_get_value( arg[next].item ) ;
	if( value && strlen(value) ) {
	    (void) strcpy( arglist[narg], arg[next].name ) ;
	    (void) strcat( arglist[narg], "=" ) ;
	    (void) strcat( arglist[narg], value ) ;
	    (void) strcat( buff, "	" ) ;
	    (void) strcat( buff, arglist[narg] ) ;
	    (void) strcat( buff, "\n" ) ;
	    narg++ ;
	}
    }
    if (remote_exec)
        strcpy( arglist[narg++], "&)" ) ;       /* tail end */

    write( LOG_DESC, buff, strlen(buff) ) ;

    argsave = arglist[narg] ;
    arglist[narg] = (char *) 0 ;

    local_exec( name, path , hostname) ;

    arglist[narg] = argsave ;
}

local_exec( name, path, hostname )
char *name, *path, *hostname;
{
    int job ;

    njob++ ;
    for( job = 0 ; jobs[job].njob ; job++ ) ;
    jobs[job].pid = lexec( path, arglist, NUL_DESC, LOG_DESC ) ;
    jobs[job].njob = njob ;
    (void) strcpy( jobs[job].name, name ) ;
    (void) strcpy( jobs[job].hostname, hostname);
    jobs[job].status = 'R';
    draw_jobs() ;
    notify_set_wait3_func(
	(char *) jobs + njob, reap, jobs[job].pid ) ;
}

Notify_value reap( handle, pid, status, rusage )
char *handle ;
int pid ;
union wait *status ;
struct rusage *rusage ;
{
    int job ;

    if( kill(pid, 0) != -1 )
	return NOTIFY_IGNORED ;

    for( job = 0 ; job < NJOBS ; job++ ) {
	if( pid == jobs[job].pid ) {
	    jobs[job].njob = 0 ;
	    jobs[job].status = 'K' ;
	    draw_jobs() ;
	    return NOTIFY_DONE ;
	    }
	}
    if (debug_level) {            
        (void) printf( "What the hell?\n" ) ;
    }
    return NOTIFY_IGNORED ;
}

mk_inputs()
{
    int next ;

    inputs = window_create( frame, PANEL,
	PANEL_SHOW_MENU,	FALSE,
	PANEL_LABEL_BOLD,	TRUE,
	PANEL_LAYOUT,		PANEL_HORIZONTAL,
	WIN_BELOW,		control,
	WIN_X,			0,
	WIN_COLUMNS,		60,
	WIN_ROWS,		PWIDTH,
	WIN_ROW_GAP,		1,
	PANEL_EVENT_PROC,	input_help_proc,
	0 ) ;

    input_scroll = scrollbar_create(0) ;

    for( next = 0 ; next < NARGS ; next++ ) {
	arg[next].item = panel_create_item( inputs, PANEL_TEXT,
	    PANEL_NOTIFY_STRING,	"\n\r",
	    PANEL_SHOW_ITEM,		FALSE,
	    PANEL_ITEM_X,		ATTR_COL(1),
	    0 ) ;
	}
}

input_help_proc( item, event )
Panel_item item ;
Event *event ;
{
    if( (event_id( event ) == MS_MIDDLE) && event_is_down( event ) )
	find_help_proc(item) ;

#if 0
    if( (event_id( event ) == MS_RIGHT) && event_is_down( event ) ) {
        help();     /* A TERRIBLE KLUDGE */
    }
#endif
    panel_default_handle_event( item, event ) ;

}

find_help_proc(item)
Panel_item item ;
{
    int i, spot ;
    char buff[NCHAR] ;
    Textsw_index first, last ;

    for( i = 0 ; i < NARGS ; i++ ) {
	if( item == arg[i].item ) {
	    window_set( help_text, TEXTSW_FILE, doc_path, 0 ) ;
	    (void) strcpy( buff, "%A " ) ;
	    (void) strcat( buff, arg[i].name ) ;
	    first = (Textsw_index) 0 ;
	    spot = textsw_find_bytes( help_text,
		&first, &last, buff, strlen(buff), 0 ) ;
	    if( spot == -1 ) {
		(void) fprintf( stderr, "cannot find %s in %s\n",
		    arg[i].name, doc_path ) ;
		return ;
		}
	    window_set( help_text, TEXTSW_FIRST, first, 0 ) ;
	    window_set( help_frame, WIN_SHOW, TRUE, 0 ) ;
#if 1
	    scrollbar_paint_bubble(
		window_get( help_text, WIN_VERTICAL_SCROLLBAR ) ) ;
#endif
	    return ;
	    }
	}
    }

mk_jobs()
{
    int next ;

    job_win = window_create( frame, PANEL,
	PANEL_SHOW_MENU,	FALSE,
	PANEL_LABEL_BOLD,	TRUE,
	PANEL_LAYOUT,		PANEL_HORIZONTAL,
	WIN_BELOW,		control,
	WIN_RIGHT_OF,		inputs,
	WIN_ROWS,		PWIDTH,
	WIN_ROW_GAP,		1,
	0 ) ;

    job_menu = menu_create( MENU_STRINGS, "stop", "cont", "kill", 0, 0 ) ;

/*  input_scroll = scrollbar_create(0) ; */

    for( next = 0 ; next < NJOBS ; next++ ) {
	job_message[next] = panel_create_item( job_win, PANEL_TEXT,
	    PANEL_ITEM_X,	ATTR_COL(1),
	    PANEL_SHOW_ITEM,	FALSE,
	    PANEL_EVENT_PROC,	job_event_proc,
	    0 ) ;
	}
    }

job_event_proc(item, event)
Panel_item item ;
Event *event ;
{
    int next, value, njob ;
    char *label ;

    if( event_id(event) != MS_RIGHT )
	return panel_default_handle_event(item, event) ;

    if( !panel_get( item, PANEL_SHOW_ITEM ) )
	return panel_default_handle_event(item, event) ;

    label = panel_get( item, PANEL_LABEL_STRING ) ;
    njob = -1 ;
    sscanf( label, "%d", &njob ) ;
    for( next = 0 ; jobs[next].njob ; next++ )
	if( jobs[next].njob == njob )
	    break ;
    if( !jobs[next].njob )
	return ;

    value = (int) menu_show( job_menu, job_win, event, 0 ) ;

    if( jobs[next].pid <= 0 )
	return ;

    if( value == 1 ) {                                 /*  pausing */
	kill( jobs[next].pid, 17 ) ;                   
	jobs[next].status = 'P';        
    } else if( value == 2 ) {                          /* continuing */
	kill( jobs[next].pid, 19 ) ;
	jobs[next].status = 'R';        	
    } else if( value == 3 ) {                          /* killing */
	kill( jobs[next].pid, 9 ) ;
	jobs[next].status = 'K';        	
    }
}

draw_jobs()
{
    int next, buff[NCHAR], job_cmp() ;

    qsort( jobs, NJOBS, sizeof(struct JOB_STRUCT), job_cmp ) ;
    for( next = 0 ; jobs[next].njob ; next++ ) {
	(void) sprintf( buff, "%4d:%c: %s (%s)",
	                jobs[next].njob, jobs[next].status, 
	                jobs[next].name, jobs[next].hostname) ;
	panel_set( job_message[next],
	    PANEL_LABEL_STRING,	buff,
	    PANEL_SHOW_ITEM,	TRUE,
	    0 ) ;
	}
    while( next < NJOBS ) {
	panel_set( job_message[next++], PANEL_SHOW_ITEM, FALSE, 0 ) ;
	}
    }

int job_cmp(A, B)
struct JOB_STRUCT *A, *B ;
{
    if( !A->njob && !B->njob ) return  0 ;
    if( !A->njob &&  B->njob ) return  1 ;
    if(             !B->njob ) return -1 ;
    if(  A->njob  >  B->njob ) return  1 ;
    if(  A->njob  <  B->njob ) return -1 ;
    return 0 ;
    }

mk_help()
{
    help_frame = window_create( frame, FRAME, WIN_ROWS, 24, 0 ) ;
    help_text = window_create( help_frame, TEXTSW,
	TEXTSW_BROWSING,	TRUE,
	TEXTSW_DISABLE_CD,	TRUE,
	TEXTSW_DISABLE_LOAD,	TRUE,
	TEXTSW_MENU,		NULL,
	0 ) ;
    }

/*		handle name field change */
name_event_proc(item, event)
Panel_item item ;
Event *event ;
{
    int value ;

    if( event_id(event) == MS_MIDDLE ) {
	window_set( help_text, TEXTSW_FILE, doc_path, 0 ) ;
	window_set( help_frame, WIN_SHOW, TRUE, 0 ) ;
	return panel_default_handle_event(item, event) ;
	}

    if( event_id(event) != MS_RIGHT )
	return panel_default_handle_event(item, event) ;

    value = (int) menu_show( name_menu, frame, event, 0 ) ;
    if( value ) {           /* Also Used vertical ALPHABET here */
	panel_set( name_text, PANEL_VALUE, sbin.files[value-1], 0 ) ;
	event_set_id( event, '\n' ) ;
	name_notify_proc( item, event ) ;
	}
    }

Panel_setting name_notify_proc(item, event)
Panel_item item ;
Event *event ;
{
    char buff[20], *name ;
    int next, count, length ;

    name = panel_get_value( name_text ) ;

    (void) strcpy( bin_path, BIN_DIR ) ;
    (void) strcat( bin_path, "/" ) ;
    (void) strcat( bin_path, name ) ;

    (void) strcpy( doc_path, DOC_DIR ) ;
    (void) strcat( doc_path, "/" ) ;
    (void) strcat( doc_path, name ) ;
    (void) strcat( doc_path, ".doc" ) ;

    (void) strcpy( def_path, DEF_DIR ) ;
    (void) strcat( def_path, "/" ) ;
    (void) strcat( def_path, name ) ;
    (void) strcat( def_path, ".def" ) ;

    count = readdoc() ;
    if( count > PWIDTH - 3 ) {
	panel_set( inputs, WIN_VERTICAL_SCROLLBAR, input_scroll, 0 ) ;
	length = (int) panel_get( arg[count].item, PANEL_ITEM_Y ) ;
	scrollbar_set( input_scroll, SCROLL_OBJECT_LENGTH, length, 0 ) ;
	}
    else
	panel_set( inputs, WIN_VERTICAL_SCROLLBAR, 0, 0 ) ;

    for( next = 0 ; next < NARGS ; next++ ) {
	panel_set( arg[next].item,
	    PANEL_SHOW_ITEM,		FALSE,
	    PANEL_VALUE,		"",
	    0 ) ;
	}
    panel_paint( inputs, PANEL_CLEAR ) ;

    for( next = 0 ; next < count ; next++ ) {
	(void) sprintf( buff, "%10s: ", arg[next].name ) ;
	panel_set( arg[next].item,
	    PANEL_LABEL_STRING,		buff,
	    PANEL_SHOW_ITEM,		TRUE,
	    0 ) ;
	}
    restore() ;
    return panel_text_notify(item, event) ;
    }

/* 	handle host name field change */
host_event_proc(item, event)
Panel_item item ;
Event *event ;
{
    int value , la;
    char buffer[NCHAR];

#if 0
    if( event_id(event) == MS_MIDDLE ) {
	window_set( help_text, TEXTSW_FILE, doc_path, 0 ) ;
	window_set( help_frame, WIN_SHOW, TRUE, 0 ) ;
	return panel_default_handle_event(item, event) ;
    }
#endif

    if( event_id(event) != MS_RIGHT )
	return panel_default_handle_event(item, event) ;

    value = (int) menu_show( host_menu, frame, event, 0 ) ;
    if( value ) {
        if (debug_level) {            
	    printf("Set a new host: %s\n",hosts.names[value-1], 0 ) ;
        }
        la = loadav(hosts.names[value-1]);
        if (debug_level) {            
            printf("The load on that host is %d\n",la);
        }
        sprintf(buffer,"%s load %5.2f",hosts.names[value-1], 0.01*((float)la));
        my_message(buffer); 
	panel_set( host_text, PANEL_VALUE, hosts.names[value-1], 0 ) ;
	event_set_id( event, '\n' ) ;
	name_notify_proc( item, event ) ;
    }
}

Panel_setting host_notify_proc(item, event)
Panel_item item ;
Event *event ;
{
    char buff[20], *host ;
    int next, count, length ;

    host = panel_get_value( host_text ) ;


/*
 *	deleted some junk ? here
 */

    restore() ;
    return panel_text_notify(item, event) ;
    }

/*		process action within the logfile name/menu */
log_event_proc(item, event)
Panel_item item ;
Event *event ;
{
    int value, show, page ;
    char buff[NCHAR] ;

    if( event_id(event) != MS_RIGHT )
	return panel_default_handle_event(item, event) ;

    if( !log_frame )
	mk_log() ;
    if( !log_frame )
	return ;
	
    show = (int) window_get( log_frame, WIN_SHOW ) ;
    page = (int) window_get( log_tty, TTY_PAGE_MODE ) ;

    value = (int) menu_show( log_menu[show + 2 * page], frame, event, 0 ) ;

    if( value == 1 )					/* show/hide */
	window_set( log_frame, WIN_SHOW, TRUE - show, 0 ) ;
    else
    if( value == 2 )					/* page/nopage */
	window_set( log_tty, TTY_PAGE_MODE, TRUE - page, 0 ) ;
    else
    if( value == 3 ) {					/* truncate file */
	close( open( LOG_SAVE, O_WRONLY | O_TRUNC, 0644 ) ) ;
	mk_log() ;
	}
    else
    if( value == 4 ) {					/* print */
	strcpy( buff, "pr -f " ) ;
	strcat( buff, LOG_SAVE ) ;
	strcat( buff, " | lpr" ) ;
	system( buff ) ;
	}
    }

Panel_setting dir_notify_proc(item, event)
Panel_item item ;
Event *event ;
{
    cwd = panel_get_value( dir_text );
    chdir( cwd );

    return panel_text_notify(item, event) ;
    }

Panel_setting log_notify_proc(item, event)
Panel_item item ;
Event *event ;
{
    LOG_FILE = panel_get_value( log_text ) ;
    if( !LOG_FILE || !LOG_FILE[0] || !strcmp( LOG_SAVE, LOG_FILE ) )
	return panel_text_notify(item, event) ;

    mk_log() ;
    return panel_text_notify(item, event) ;
    }

int readdoc()
{
    char buff[NCHAR] ;
    FILE *f ;
    int next ;

    f = fopen( doc_path, "r" ) ;
    if( f == NULL )  {
	(void) fprintf( stderr, "readdoc: cannot open: %s\n", doc_path ) ;
	return 0 ;
	}

    next = 0 ;
    while( fgets(buff, NCHAR-1, f) ) {
	if( sscanf( buff, "%%A %s", arg[next].name ) == 1 ) {
	    next++ ;
	    }
	}
    arg[next].name[0] = (char) NULL ;
    (void) fclose(f) ;
    return next ;
    }

static void save()
{
    FILE *f ;
    char *value ;
    int next ;

    if( !def_path )
	return ;

    f = fopen( def_path, "w" ) ;
    if( f == NULL )  {
	(void) fprintf( stderr, "save: cannot open %s\n", def_path ) ;
	return ;
	}
    for( next = 0 ; arg[next].name[0] ; next++ ) {
	value = panel_get_value( arg[next].item ) ;
	if( strlen( value ) )
	    (void) fprintf( f, "%10s = %s\n", arg[next].name, value ) ;
	}

    (void) fclose(f) ;
    }

restore()
{
    FILE *f ;
    char buff[NCHAR], keyword[32] ;     /* let's have keywords max 32 chars */
    int i, next , nval;
    extern int errno ;

    f = fopen( def_path, "r" ) ;
    if( f == NULL ) {
	if( errno != 2 )
	    (void) fprintf( stderr, "restore: cannot open: %s\n", def_path ) ;
	return ;
	}

    while( fgets(buff, NCHAR-1, f) ) {
	(void) sscanf( buff, "%s", keyword ) ;      /* get keyword */
        for (i=0;keyword[i]!=NULL;i++)     /* make sure there's no '=' left */
            if (keyword[i] == '=') {        /* in keyword */
                keyword[i] = '\0';          /* this is to catch 'key=val' */
                break;                      /* as opposed to ' key = val' */
            }
        for (nval=0;buff[nval]!=NULL;nval++)/* find where value begins: */
	    if (buff[nval] == '=') {        /* locate the '=' */
                while (buff[++nval] == ' ') /* skip until first non-blank */
                    ;
                break;                      /* &buff[nval] is now value */
            }

/*	fprintf(stderr,"scanning, keyword=%s\n",keyword);	/* DEBUG */
/*        fprintf(stderr,"value @ %d in:%s\n",nval,buff);         /* DEBUG */
	for( next = 0 ; arg[next].name[0] ; next++ ) {
	    if( strcmp( keyword, arg[next].name ) == 0 ) {
		for( i = 0 ; i < NCHAR ; i++ )
		    if( buff[i] == '\n' )
			buff[i] = '\0' ;
		panel_set( arg[next].item, PANEL_VALUE,	&buff[nval], 0 ) ;
		break ;
		}
	    }
	}
    (void) fclose(f) ;
    }

int open_log()
{
    int newfile ;

    LOG_FILE = panel_get_value( log_text ) ;
    if( !LOG_FILE || !LOG_FILE[0] )
	return -1 ;

    newfile = open( LOG_FILE, O_RDWR | O_CREAT | O_APPEND, 0644 ) ;
    if( newfile < 0 ) {
	(void) fprintf( stderr, "unable to open %s\n", LOG_FILE ) ;
	return -1 ;
	}

    (void) close( LOG_DESC ) ;
    LOG_DESC = newfile ;
    strcpy( LOG_SAVE, LOG_FILE ) ;
    return 1 ;
    }

#define XOFFSET	0
#define YOFFSET	420

mk_log()
{
    static char *args[] = { "tail", "-f", (char *) 0, (char *) 0 } ;
    int x, y, show, rows, cols ;

    if( open_log() < 0 )
	return ;
    show = TRUE ;
    args[2] = LOG_FILE ;
    if( log_frame ) {
	window_set( log_frame, FRAME_NO_CONFIRM, TRUE, 0 ) ;
	show = (int) window_get( log_frame, WIN_SHOW ) ;
	x = (int) window_get( log_frame, WIN_X ) ;
	y = (int) window_get( log_frame, WIN_Y ) ;
	rows = (int) window_get( log_frame, WIN_ROWS ) ;
	cols = (int) window_get( log_frame, WIN_COLUMNS ) ;
    if( x < 0 ) x = XOFFSET ;
    if( y < 0 ) y = YOFFSET ;
	window_destroy( log_frame ) ;
	}
    else {
	x = XOFFSET ;
	y = YOFFSET + 5 ;
	rows = 16 ;
	cols = 100 ;
	}

    log_frame = window_create( frame, FRAME,
	WIN_ROWS,	rows,
	WIN_COLUMNS,	cols,
	WIN_X,		x,
	WIN_Y,		y,
	0 ) ;

    log_tty = window_create( log_frame, TTY, TTY_ARGV, args, 0 ) ;
    window_set( log_frame, WIN_SHOW, show, 0 ) ;
    }

my_message( string )
char *string ;
{
    if( string ) {
	panel_set( message,
	    PANEL_LABEL_STRING,	string,
	    PANEL_SHOW_ITEM,	TRUE,
	    0 ) ;
	}
    else {
	panel_set( message, PANEL_SHOW_ITEM, FALSE, 0 ) ;
	}
    }

/************************************************************************
 *									*
 *  hsearch -- get a list of all the hostnamed from the rwho directory  *
 *									*
 *	char *dname		the directory to be searched		*
 *	char *flist[]		a NULL terminated list of file names	*
 *									*
 ************************************************************************/
char *malloc() ;

int hsearch( dname, flist, nmax)
char *dname, *flist[] ;
int   nmax;
{
    DIR *dir ;
    struct direct *fentry ;
    int next, cmp() ;

    dir = opendir( dname ) ;
    if( !dir ) {
	perror( dname ) ;
	return 0 ;
	}

    next = 0 ;
    while( fentry = readdir( dir ) ) {
	if( fentry->d_name[0] != '.' ) {
	    if ((next+2) >= nmax) {
		printf("Warning: only first %d entries in directory %s read\n",
			nmax, dname);
		break;
            }
	    flist[next] = malloc( fentry->d_namlen + 1 ) ;
                                                    /* skip whod. */
                                                    /*     v      - a kludge */
	    (void) strcpy( flist[next++], (fentry->d_name)+5  ) ; 
	}
    }
    qsort( (char *) flist, next, sizeof(char *), cmp ) ;
    flist[next++] = (char *) NULL ;
    flist[next++] = (char *) NULL ;
    closedir( dir ) ;

    if (next == 1)
        printf("Warning: no files in remote host directory %s\n",dname);
    return next ;
    }
/*
 *  lexec -- local exec call
 *
 */
int lexec(path, alist, inp_desc, out_desc)
char *path, *alist[] ;
int inp_desc, out_desc ;
{
    int pid, i ;
    static char error[] = "unable to execute " ;

    for (i=0; alist[i]!=NULL; i++)               /* debug */
        if (debug_level) {
            printf("exec-debug>%s\n",alist[i]);
        }
    if( pid = vfork() )
	return pid ;

    (void) dup2( inp_desc, 0 ) ;		/* reassign i/o */
    (void) dup2( out_desc, 1 ) ;
    (void) dup2( out_desc, 2 ) ;

    for( i = 3 ; i < 32 ; i++ )			/* close everything else */
	(void) close(i) ;

    i = open( "/dev/tty", 0 ) ;			/* detach from terminal */
    (void) ioctl( i, TIOCNOTTY, 0 ) ;
    (void) close(i) ;

	    /* execvp so that search paths are used for "rsh" */
    execvp( path, alist ) ;			/* execute it */

    (void) write( 2, error, strlen(error) ) ;	/* write error message */
    (void) write( 2, path, strlen(path) ) ;
    (void) write( 2, "\n", 1 ) ;
    _exit(127) ;				/* exit if we bombed */
    }
/*
 *	gets the load average of a remote host under investigation
 *			29-jul-89	Peter Teuben
 */
loadav(name)
char *name;
{
    char        fname[128];
    FILE        *fp;
    struct whod w;
    int         n;

    strcpy(fname,"/usr/spool/rwho/whod.");
    strcat(fname,name);

    if ((fp=fopen(fname,"r"))==NULL) {
        fprintf(stderr,"File %s cannot be opened\n", fname);
        return(-1);
    }
    
                /* read first part of structure until the 'whoent' struct */
    n  = read (fileno(fp), &w, 4*sizeof(int) + 36*sizeof(char));
    fclose(fp);
    if (n < 0) {
        fprintf(stderr,"Error reading %s\n",fname);
        return(-1);
    }
    if (debug_level) {
        printf("Load averages on %s are: %d %d %d boottime=%d\n",
                w.wd_hostname, 
                w.wd_loadav[0], w.wd_loadav[1], w.wd_loadav[2],
                w.wd_boottime);
    }
    return(w.wd_loadav[0]);
}
/************************************************************************
 *									*
 *  search -- get a list of all the files in a directory		*
 *									*
 *	char *dname		the directory to be searched		*
 *	char *flist[]		a NULL terminated list of file names	*
 *									*
 ************************************************************************/
char *malloc() ;

int search( dname, flist, nmax)
char *dname, *flist[] ;
int   nmax;
{
    DIR *dir ;
    struct direct *fentry ;
    int next, cmp() ;

    dir = opendir( dname ) ;
    if( !dir ) {
	perror( dname ) ;
	return 0 ;
	}

    next = 0 ;
    while( fentry = readdir( dir ) ) {
	if( fentry->d_name[0] != '.' ) {
            if ( (next+2) >= nmax) {
                printf("Warning: could only read %d entries from %s\n",
                        nmax, dname);
                break;
            }
	    flist[next] = malloc( fentry->d_namlen + 1 ) ;
	    (void) strcpy( flist[next++], fentry->d_name ) ;
	}
    }
    qsort( (char *) flist, next, sizeof(char *), cmp ) ;
    flist[next++] = (char *) NULL ;
    flist[next++] = (char *) NULL ;
    closedir( dir ) ;

    if (next == 1)
        printf("Warning: no files found in directory %s\n",dname);
    return next ;
}

static int cmp(a, b)
char **a, **b ;
{
    return strcmp( *a, *b ) ;
}


/*
 *  This is a nice feature to avoid seeing executables which do not
 *  have documentation files, perhaps in debug mode this should
 *  warn user of programs which do not have doc's and doc's which
 *  do not have exe's 
 *
 *  It returns the new number of executables, i.e. <= ne.
 */

int cmpbin_doc( elist, ne, dlist, nd)
char *elist[], *dlist[];
int   ne, nd;
{
    int ie, id, i, n;
    char exedoc[128], *ep, *dp;

    for (ie=0, id=0; ; ) {
        strcpy(exedoc,elist[ie]);
        strcat(exedoc,".doc");
        n=strcmp(exedoc,dlist[id]);
        if (n==0) {
            ie++;
            id++;
        } else if (n<0) {
            if (debug_level)
                printf("Warning: executable %s has no doc file\n",elist[ie]);
#if 1
            i = ie;
            do {
                elist[i] = elist[i+1];
            } while (elist[++i] != NULL);
#else
	    ie++
#endif
        } else {
            if (debug_level)
                printf("Warning: no executable for doc file %s\n",dlist[id]);
            id++;
        }
        if (elist[ie]==NULL)
            break;
        if (dlist[id]==NULL)
            break;
    } 
    return ne;
}

int nop()
{
	return;
}
