/*
 *    crck.c  --  Unix version of the CP/M utility
 *
 *    Usage:
 *        crck [-t|c|u|i] [filename ... ]
 *
 *    Calculates a 16-bit cyclic redundancy check on its input files and
 *    prints the resulting number in hexadecimal.  CRCK is manually invoked
 *    on two versions of the same file to ensure that they are the same.
 *    Prints the the filename, size in bytes or sectors depending on the
 *    option flag, and the CRC.
 *
 *    Options:
 *        -t    Perform CRCK on a CP/M text file that is stored as a Unix
 *              file.  CP/M record delimters are substituted for Unix 
 *		newlines and the file is padded to a 128-byte sector
 *  		boundary with ^Z's.
 *
 *	  -u    Perform CRCK on a regular Unix file.  Takes whatever is there.
 *
 *	  -c    Perform CRCK on a CP/M ".COM" file.  Handled the same as
 *		the -u flag except sectors are reported rather than bytes.
 *
 *	  -i    Perform CRCK on a CP/M ".COM" file stored in ITS binary 
 *		format.  Ignores the first four bytes of the file.  Implies
 *  		the -c option.
 * 
 *	  If no option flag is chosen, the default option is -u.
 *    
 *    Wildcards may, of course, be used and the program may be used in 
 *    pipelines.
 *
 *
 *----------------------------------------------------------------------------
 *    Version 2.0  mods by Ben Goldfarb (decvax!ucf-cs!goldfarb or
 *                                            Goldfarb.ucf-cs @ Rand-Relay)
 *
 *		   Added CP/M text and binary file modes and changed to 
 *		   use standard stream I/O.
 *
 *    Version 2.1  mods by Peter Teuben :
 *		   at least don't crash when no cmdln arguments are given!!
 */

#include <stdio.h>

#define VERSION 20		/*  Version 2.0				*/


#define SECSIZ	128		/*  CP/M sector size for fill  		*/
#define CTLZ 	'Z' - 0x40	/*  Fill character for text files	*/ 

int cflag = 0;			/*  CP/M binary file			*/
int tflag = 0;			/*  CP/M text file			*/
int uflag = 0;			/*  Unix file				*/
int iflag = 0;			/*  ITS Binary file			*/

unsigned short crck();		/*  CRC calculation function		*/


main(argc, argv)
int argc;
char **argv;
{
	int err = 0;

	if (argc == 1 )
		usage();

	if ((++argv)[0][0] == '-') {
		switch (argv[0][1]) {
		  case 't':
			++tflag;
			break;
		  case 'c': 
			++cflag;
			break;
		  case 'u':
			++uflag;
			break;
		  case 'i':
			++iflag; ++cflag;
			break;
		  default: {
			printf ("not a valid switch\n");
			usage();
		  }
		}
		--argc; ++argv;
	}
	else  /* for right now, "Unix mode" is the default */
		++uflag;

	if (isatty(2)) {  /*  print header if not in a pipe  */
		printf("CRCK program for Unix\n");
		printf("Version %1d.%1d\n\n", VERSION/10, VERSION%10);
		if (iflag) printf("ITS mode selected\n");
		if (tflag) printf("CP/M text mode selected\n");
		if (cflag) printf("CP/M \".COM\" mode selected\n");
		if (uflag) printf("Unix file mode selected\n");
		putchar('\n');
	}

	/*  we're doing a named file or a list of same  */
	for (; --argc; ++argv) {
		err = 0;
		if (freopen(*argv, "r", stdin) == NULL) {
			perror(*argv);
			++err;
		}
		if (!err)
			docrck(*argv);
	}
}

usage()
{
	printf("Usage: crck [-t|c|u|i] [filename ... ]\n");
	exit(0);
}
	
docrck(name) /* Calculate and print a "CRCK" for stdin stream. */
char *name;
{
	int icount;
	long nsec;
	unsigned short oldcrck;
	unsigned char crbuf[SECSIZ];
	unsigned char *bufend = crbuf + SECSIZ;

	register c;
	register unsigned char *cp;

	oldcrck = 0;
	nsec = 0L;
	cp = crbuf;
	if (tflag) {
	    /*
   	     *  To calculate CRC for CP/M text files that are
	     *  presently stored in the form of Unix text files, 
	     *  an ASCII CR must be inserted in the buffer before 
             *  each '\n'.  Then the file must  be padded out to 
	     *  an even SECSIZ boundary with CTLZ's.
	     */
            while ((c = getchar()) != EOF) {
		if (((*cp++ = c) == '\n') && tflag) {  /*  need to add \r */
			*(cp - 1) = '\r';
			if (cp >= bufend) {  /* sector boundary */
				oldcrck = crck(crbuf, SECSIZ, oldcrck);
				++nsec;
				cp = crbuf;
			}
			*cp++ = '\n';
		}
		if (cp >= bufend) {  /* sector boundary */
			oldcrck = crck(crbuf, SECSIZ, oldcrck);
			++nsec;
			cp = crbuf;
		}
	    }

	    if (cp != crbuf) { /* EOF during partial sector  */
		register fillcnt;

		/*  Fill from current pointer pos to end of buffer  */
		if (tflag)
			for (fillcnt = bufend - cp; fillcnt--; *cp++ = CTLZ)
					;

		/*  Need to update CRC for final sector. */
		oldcrck = crck(crbuf, cp-crbuf, oldcrck);
		++nsec;
	    }
	}
	else {  /* uflag or cflag  -- do a quick crck  */
	    unsigned char tbuf[BUFSIZ];
	    int count = 0;
	    
	    if (iflag) /* skip first four bytes of ITS binary file */
		for (icount = 0; 
		    (c = getchar()) != EOF && icount < 4;
		    icount++)
			;
	    while ((count = fread(tbuf, 1, BUFSIZ, stdin)) > 0) {
		oldcrck = crck(tbuf, count, oldcrck);
		if (!cflag)
			nsec += count;  /* here nsec is for BYTES  */
		else
			nsec += ((count + SECSIZ - 1) / SECSIZ);
	    }
	}

	printf("%-18s", name);
	printf(" (%7ld %-7s)   ==>   %04x\n", 
		nsec, uflag ? "bytes" : "sectors", oldcrck);
}

unsigned short
crck(crbuf, count, ldcrc)  /* uses algorithm of CRCK.COM prior to V5.0 */
register unsigned char *crbuf;
register count;
register unsigned short ldcrc;
{
/*
 *  This routine copyright (c) 1983, Computer Development, Inc., 
 *  Beaverton, OR, USA.
 *  All rights reserved.
 */
	register unsigned short n;

	while (--count >= 0) {
		n = (ldcrc << 1);
		n = (n & 0xff00) | (( n + *crbuf++ ) & 0x00ff);
		if(ldcrc & 0x8000)
			n ^= 0xa097;
		ldcrc = n;
	}
	return(ldcrc);		
}
