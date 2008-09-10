/*
 *	NEMO:  simple intro for poor users
 *
 *	23-oct-87	1.0: first version written	PJT
 *	25-Mar-88	1.1: lot's of more output for beginning user
 *	  -Apr-88	1.2: split up into several options
 *	21-Apr-88	1.3: help-level implemented
 *	21-May-88	1.4: define option implemented
 *	 7-Jun-88	1.5: added HISTORY env var
 *	 2-dec-88	1.6: online documentation upgraded
 *	19-jun-89	1.7: added yapp & help strings 
 *	 7-jul-89	1.8: added NEMOTOOL stuff
 *	21-nov-89	1.8a: added _trace_ #define
 *	14-dec-89	1.8b: help vector a bit expanded	PJT
 *	11-sep-90	2.0: changed name to 'nemoshow'		PJT
 *	10-dec-90       2.1: few new things for new release	PJT
 *	19-feb-92	2.1a: usage
 *	24-mar-97	3.0: NEMOPATH is deprecated now		pjt
 *      10-sep-08       3.1: fixed for g++ but deprecated       pjt
 */

#include <stdinc.h>		/* also gets <stdio.h>	*/
#include <getparam.h>

extern char     **environ;	/* environment variables */
				
extern int	debug_level;	/* see: getparam.c ?? or dprintf.c */
extern int	help_level; 	/* see: getparam.c */
extern int      yapp_dev;  	/* see: getparam.c */
extern int	history;  	/* see: getparam.c ?? or history.c */

extern string   yapp_string;	/* see: getparam.c */
extern string   help_string;    /* see: getparam.c */

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "env=f\n		List some used environment variables",
    "def=f\n		List various defaults",
    "VERSION=3.1\n	10-sep-08 PJT",
    NULL,
};

string usage = "show NEMO intrinsics";

bool  bin, man;		/* see directory listing */

char *help[]={
 "NEMO is a collection of programs designed to help you in doing stellar",
 "dynamics.  All programs have a common user-interface, based on a",
 "commandline with keyword=value parameters.  Those keywords and their",
 "default values (if they have one) are listed if you try 'nemo_command",
 "help='.  In this case 'help' is a so called hidden system keyword.  NEMO",
 "uses a few of such system keywords (help,yapp,debug,host).",
 "The 'help=h' option may list a short explanation for keywords.\n",
 "The directory $NEMOBIN should contain all possible nemo_command's.\n",
 "The userinterface is not always very friendly with respect to illegal",
 "input, so usually looking up the help file using (UNIX) man command can be",
 "helpfull.  i.e.  'man nemo_command' to look up documentation for",
 "nemo_command.",
 "\n\n",
 "Type 'nemo help=' for additional options to this command\n\n",
   NULL,
};

void show_env(void);
void show_def(void);
void show (char *e, char *s);

extern int nemo_history;                        /* see history.c */

int main(int argc, char *argv[])		/* This is an example where we don't want nemo_main */
{
    int i;
    char *np;
    
    i=0;
    if (argc==1) {			/* give help when no arguments given */
	while (help[i]!=NULL)
	   printf ("%s\n",help[i++]);
        np = getenv("NEMO");
        if (np==NULL)
            printf ("No NEMO environment set (NEMO=(null))\n");
        else
            printf ("NEMO = %s\n",np);
    }
        	
    initparam(argv, defv);              /* initialize command line pars */

    if (getbparam("env"))
	show_env();

    if (getbparam("def"))
	show_def();
    
    return 0;
}

void show_env(void)
{
    printf ("environment variable which may be used in NEMO:\n");
    show("NEMO","where NEMO is");
    show("NEMOHOST","useful in multi-CPU environment");
    show("NEMOBIN","directory for executables");
    show("NEMOLIB","directory for libraries");
    show("NEMOOBJ","root directory for 'loadobj' object files");
    show("MONGOPATH","where MONGO is (must be defined for yapp_mongo");
    show("NUMRECPATH","where NUMREC is (some routines need it at compile time");
    show("YAPPLIB","default link library yapp");
    show("NEMONEWS","addnews/readnews/scannews: see news.1");
    show("BTRPATH","body transformations");
    show("POTPATH","potential functions");
    show("FLOAT_OPTION","sun specific auto float_option compilation");
    show("INCLUDE","extra include directory for #include <...> files");
    show("LIBRARY","extra lib directory for resolving -l... cc-switches");
    show("YAPP","default yapp_mongo plot device");
    show("DEBUG","default level of helpful debug output");
    show("HELP","default help level");
    show("HISTORY","need to write history?");    
    printf ("\n Hidden system keywords:\n\n");
    printf ("(YAPP)       yapp= %d   - string = %s\n",yapp_dev,yapp_string);
    printf ("(DEBUG)     debug= %d\n",debug_level);
    printf ("(HELP)       help= %d   - string = %s\n",help_level,help_string);
    printf ("(HISTORY) history= %d\n",nemo_history);
}

void show (char *e, char *s)
{
   int i, n;

   for (i=0; environ[i]!=NULL; i++) {
      n = strlen(e);
      if (strncmp(e,environ[i],n)==0 && *(environ[i]+n)=='=') {
         printf ("%15s set to %s\n",e,environ[i]+n+1);
	 dprintf (1," %s\n",s);
	 return;
      }
   }
   printf ("%15s not set\n",e);
   dprintf (1," %s\n",s);
}

void show_def(void)
{
   printf ("\nNow follow a few check on #define's:\n\n");
   printf ("\nThey can be used in code to isolate system dependancies\n");
   printf ("Either explicitely through a cc-compilation option -Dsystem\n");
   printf ("or by adding an appropriate #if defined() in options.h\n\n");

#if defined(sun)
   printf ("sun     is defined\n");
#else
   printf ("sun     is not defined\n");
#endif
#if defined(sun3)
   printf ("sun3    is defined\n");
#else
   printf ("sun3    is not defined\n");
#endif
#if defined(sun4)
   printf ("sun4    is defined\n");
#else
   printf ("sun4    is not defined\n");
#endif
#if defined(gould)
   printf ("gould   is defined\n");
#else
   printf ("gould   is not defined\n");
#endif
#if defined(convex)
   printf ("convex  is defined\n");
#else
   printf ("convex  is not defined\n");
#endif
#if defined(alliant)
   printf ("alliant is defined\n");
#else
   printf ("alliant is not defined\n");
#endif
#if defined(turboc)
   printf ("turboc  is defined\n");
#else
   printf ("turboc  is not defined\n");
#endif
#if defined(unicos)
   printf ("unicos  is defined\n");
#else
   printf ("unicos  is not defined\n");
#endif
#if defined(unixpc)
   printf ("unixpc  is defined\n");
#else
   printf ("unixpc  is not defined\n");
#endif
#if defined(_trace_)
   printf ("_trace_ is defined\n");
#else
   printf ("_trace_ is not defined\n");
#endif
#if defined(bsd)
   printf ("bsd     is defined\n");
#else
   printf ("bsd     is not defined\n");
#endif
#if defined(sys5)
   printf ("sys5    is defined\n");
#else
   printf ("sys5    is not defined\n");
#endif
#if defined(mc68k)
   printf ("mc68k    is defined\n");
#else
   printf ("mc68k    is not defined\n");
#endif
#if defined(SAVE_OBJ)
   printf ("SAVE_OBJ    is defined\n");
#else
   printf ("SAVE_OBJ    is not defined\n");
#endif
#if defined(NEMO_INC)
   printf ("NEMO_INC    is defined\n");
#else
   printf ("NEMO_INC    is not defined\n");
#endif
#if defined(DLOPEN)
   printf ("DLOPEN    is defined\n");
#else
   printf ("DLOPEN    is not defined\n");
#endif
#if defined(NUMREC)
   printf ("NUMREC    is defined\n");
#else
   printf ("NUMREC    is not defined\n");
#endif
}

