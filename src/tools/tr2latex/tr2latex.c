/*
** tr2latex - troff to LaTeX converter
** $Id$
** COPYRIGHT (C) 1987 Kamal Al-Yahya, 1991,1992 Christian Engel
** 
** Module: tr2latex.c
**
** This module contains the main function inititating the translation
** and supporting the Usage.
*/

char *usage_doc[] = {
"tr2latex - troff to LaTeX converter, $Revision$",
"SYNTAX:  tr2latex [-m] [-t] [-<n>] [-s <style>] [-o <outfile>] [<file>...]",
"options: -m            for manual",
"         -t            twoside page style",
"         -<n>          a number n gives the font size (default is 12pt",
"                       for man, 11pt otherwise)",
"         -s <style>    use documentstyle <style> (default is article)",
"         -o <outfile>  send output to <outfile> (default is stdout)",
#ifdef DEBUG
"         -do           debug output",
"         -dv           verbose informations",
#endif
0
};

#define MAIN

#include	"setups.h"
#include	"protos.h"

#ifdef __TURBOC__
/*--- The Turbo C stack needs to be increased in size,
      6k seems to be enough; 4k is not                 ---*/
extern unsigned _stklen = 0x1800 ;
#endif


extern char *mktemp();
char scratch_file[MAXWORD];

bool	man,		/* option -m: manual page */
		fontsize,	/* option -9/-10/-11/-12: font size */
		twoside,	/* option -t: twoside */
		piped_in;

char *document = "article";	/* document type, see also -s option */
extern char version [];

FILE *out_file;		/* in case they can't use *NIX redirecting or piping */

char *prgname;
char inbuf[MAXLEN],
     outbuf[MAXLEN];



int main (int argc, char *argv[])
{
	char *pin = inbuf,
	     *pout = outbuf;
	FILE *in_file;
	long timeval;		/* clock value from time() for ctime()	*/
	int argi;

	prgname = argv [0];
	out_file = stdout;		/* default output */

	/* process option flags */
	getopts (&argc, argv);

#ifdef DEBUG
	if (debug_v) {
		printf ("Arguments:");
		for (argi = 1; argi < argc; argi++)
			printf (" %s", argv [argi]);
		putchar ('\n');
	}
#endif

#ifdef	VMS
	if ((out_file == stdout) &&
		((out_file = freopen ("Sys$Output:", "w", stdout,
							  "rat=cr", "rfm=var")) == NULL)) {
		fprintf (stderr, "Can't reopen stdout\n");
	    errexit (GOOD);
	}
#endif

	/* initialize spacing and indentation parameters */
	strcpy(linespacing.def_units,"\\normalbaselineskip");
	strcpy(linespacing.old_units,"\\normalbaselineskip");
	strcpy(indent.def_units,"em");
	strcpy(indent.old_units,"em");
	strcpy(tmpind.def_units,"em");
	strcpy(tmpind.old_units,"em");
	strcpy(space.def_units,"\\baselineskip");
	strcpy(space.old_units,"\\baselineskip");
	strcpy(vspace.def_units,"pt");
	strcpy(vspace.old_units,"pt");
	linespacing.value = 1.;
	linespacing.old_value = 1.;
	indent.value = 0.;
	indent.old_value = 0.;
	tmpind.value = 0.;
	tmpind.old_value = 0.;
	space.value = 1.;
	space.old_value = 1.;
	vspace.value = 1.;
	vspace.old_value = 1.;
	linespacing.def_value = 0;
	indent.def_value = 0;
	tmpind.def_value = 0;
	space.def_value = 1;
	vspace.def_value = 1;
	
	math_mode = 0;					/* start with non-math mode */
	de_arg = 0;                     /* not a .de argument */
	
	/* start of translated document */
	
	timeval = time(0);
	fprintf (out_file,
"%% -*-LaTeX-*-\n\
%% Converted automatically from troff to LaTeX\n\
%% by %s\n\
%% on %s\
%% tr2latex was written by Kamal Al-Yahya at Stanford University\n\
%% (Kamal%%Hanauma@SU-SCORE.ARPA)\n\
%% and substantially enhanced by Christian Engel at RWTH Aachen\n\
%% (krischan@informatik.rwth-aachen.de).\n\
%%\n\
%% troff input file%s:%s",
			 version, ctime(&timeval),
			 argc>2?"s":"",
			 argc==1?" <stdin>":"");
	for (argi = 1; argi < argc; argi++) {
		if (strcmp (argv [argi], "-") == 0)
			fprintf (out_file, " <stdin>");
		else
			fprintf (out_file, " %s", argv[argi]);
	}

	/* document style and options */
	fprintf (out_file,"\n\n\\documentstyle[%s", man? "troffman": "troffms");
	if (fontsize == 0 && !man)
		fontsize = 11;
	if (fontsize != 0)
	fprintf (out_file,",%dpt", fontsize);
	if (twoside)
		fputs (",twoside", out_file);
	fprintf (out_file,"]{%s}\n\\begin{document}\n", document);

	if (argc == 1)
		process (stdin, "<stdin>", pin, pout);
	else {
		for (argv++; --argc; argv++) {
			if (strcmp (*argv, "-") == 0)
				process (stdin, "<stdin>", pin, pout);
			else if ((in_file = fopen(*argv,"r")) == NULL)
				fprintf(stderr,"%s: Cannot open input file `%s'\n",
						prgname,*argv);
			else {
				process (in_file, *argv, pin, pout);
				fclose(in_file);
			}
		}
	}
	/* close translated document */
	fputs("\\end{document}\n",out_file);

	exit(GOOD);
}


void process (FILE *in_file, char *f_name, char *pin, char *pout)
{
	static char sep[] = "--------------------------------------------------";

	tmpbuf (in_file, pin);
	fprintf (out_file, "%%%s\n%% start of input file: %s\n%%\n", sep, f_name);
	troff_tex (pin, pout, 0, 0);
	fputs (pout, out_file);
	fprintf (out_file, "%%\n%% end of input file: %s\n%%%s\n", f_name, sep);
}


# define shift_arg(cnt)         for (i = argind; i < *p_argc - cnt; i++)    \
                                    argv [i] = argv [i + cnt];              \
                                argind -= cnt;                              \
                                *p_argc -= cnt;

void getopts (int *p_argc, char *argv[])
{
	char *p_opt;
	int argind, i;

    for (argind = 1; argind < *p_argc; argind ++)
        if (*argv [argind] == '-') {
			if (strcmp (argv+1, "help") == 0)
				usage (0);
            for (p_opt = argv [argind] + 1; ; p_opt++) {
                switch (*p_opt) {
                case '\0':
                    if (p_opt != argv [argind] + 1) {
						shift_arg (1);
					}
                    break;
				case 'm':
					man = 1;
					continue;
				case 't':
					twoside = 1;
					continue;
				case 's':
					document = argv [argind+1];
					shift_arg(2);
					continue;
				case '9':
					fontsize = 9;
					continue;
				case 'o':
					if (argind + 1 >= *p_argc) {
						fprintf (stderr, "%s: Missing output file name\n",
								 prgname);
						usage (1);
					}
#ifdef	VMS
					if ((out_file = fopen(argv[argind+1],"w",
										  "rat=cr","rfm=var")) == NULL)
#else
					if ((out_file = fopen(argv[argind+1],"w")) == NULL)
#endif
					{
						fprintf(stderr, "%s: can't open %s\n",
						  prgname, argv [argind+1]);
						usage (errno);
					}
					shift_arg(2);
					break;
#ifdef DEBUG
				case 'd':
					if (*++p_opt == '\0') {
						fprintf (stderr, "%s: missing debug option\n",
								 prgname);
						usage (1);
					}
					switch (*p_opt) {
					case 'o':
						debug_o = 1;
						break;
					case 'v':
						debug_v = 1;
						break;
					default:
						fprintf (stderr, "%s: unknown debug option %c\n",
								 prgname, *p_opt);
						usage (1);
					}
					continue;
#endif
                case '?':
                    usage (0);

				case '1':
					if (isdigit (p_opt [1])) {
						fontsize = 10 + *++p_opt - '0';
						continue;
					}
					/* NO BREAK */
                default:
					fprintf (stderr,"%s: Unknown option %c\n",prgname,*p_opt);
                    usage (1);
                }
                break;
            }
		}
}


void usage (int exitcode)
{
	int i;

	printf ("tr2latex (c) 1986/1987 Kamal Al-Yahya, 1991 C. Engel\nVersion %s\n",
			version);
	for (i=0; usage_doc [i]; i++)
		printf ("%s\n", usage_doc [i]);
	exit (exitcode);
}
	

void errexit (int exitcode)
{
	fprintf (stderr, "%s: Error #%03d ", prgname, exitcode);
	exit (exitcode);
}


