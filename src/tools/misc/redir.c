/* redir
 * a program to execute a command, redirecting stdout and stderr separately,
 * a feature available in sh but not csh.
 *
 * THIS SOFTWARE EXISTS IN THE PUBLIC DOMAIN.  ABSOLUTELY NO
 * RESTRICTIONS ON ITS USAGE, MODIFICATION OR DISTRIBUTION
 * EXIST.  NO GUARANTEES ARE MADE CONCERNING THE PROPER
 * FUNCTIONING OF THIS SOFTWARE.
 *
 * Usage:
 *
 * redir [options] command-words...
 *
 * where options are
 *  -o file
 *   to redirect stdout to file "file," without clobbering "file" if it exists
 *  -e file
 *   to redirect stderr to file "file," without clobbering "file" if it exists
 *  -O file
 *   to redirect stdout to file "file," clobbering "file" if it exists
 *  -E file
 *   to redirect stderr to file "file," clobbering "file" if it exists
 *  -a
 *   the following redirection option is to append to the named file, as in:
 *    redir -O foo.log -ae foo.err command command-arg1 command-arg2 ...
 *   (stdout (destructively) to foo.log, stderr appended to foo.err.
 *   The -a option doesn't care if the named file exists or not.  This option
 *   must be used twice if you want stdout and stderr both to be appended
 *   to their respective files;
 *    redir -ao foo.log -ae foo.err command ...
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/file.h>

#define USAGESTRING ("Usage: %s [-[a][oeOE] file] ... command-words\n")

main(argc, argv)
int             argc;
char          **argv;
{
    extern int      optind;
    extern char    *optarg;
    char           *stderrfile = NULL, *stdoutfile = NULL;
    int             rstderr = 0, rstdout = 0, cstderr = 0, cstdout = 0, astderr
= 0, astdout = 0, appendopt = 0, c, fd;

	if (argc==1) {
		printf("redir [-[a][oeOE] [file] ... command arg1 arg2 ...\n");
		exit(1);
	}

    while ((c = getopt(argc, argv, "ao:e:O:E:")) != EOF) {
       switch (c) {
           case 'a':
               if (appendopt) {
                   fprintf(stderr, USAGESTRING, argv[0]);
                   exit(1);
               }
               ++appendopt;
               break;
           case 'e':
               if (rstderr) {
                   fprintf(stderr, USAGESTRING, argv[0]);
                   exit(1);
               }
               ++rstderr;
               stderrfile = optarg;
               if (appendopt) {
                   appendopt = 0;
                   ++astderr;
               }
               break;
           case 'o':
               if (rstdout) {
                   fprintf(stderr, USAGESTRING, argv[0]);
                   exit(1);
               }
               ++rstdout;
               stdoutfile = optarg;
               if (appendopt) {
                   appendopt = 0;
                   ++astdout;
               }
               break;
           case 'E':
               if (rstderr) {
                   fprintf(stderr, USAGESTRING, argv[0]);
                   exit(1);
               }
               ++rstderr;
               ++cstderr;
               stderrfile = optarg;
               if (appendopt) {
                   appendopt = 0;
                   ++astderr;
               }
               break;
           case 'O':
               if (rstdout) {
                   fprintf(stderr, USAGESTRING, argv[0]);
                   exit(1);
               }
               ++rstdout;
               ++cstdout;
               stdoutfile = optarg;
               if (appendopt) {
                   appendopt = 0;
                   ++astdout;
               }
               break;
           default:
               fprintf(stderr, USAGESTRING, argv[0]);
               exit(1);
               break;
       }
    }
    if (rstdout) {
       if ((fd = open(stdoutfile, (astdout ?
                                   (O_WRONLY | O_APPEND | O_CREAT) :
                                   (cstdout ?
                                    (O_WRONLY | O_CREAT | O_TRUNC) :
                                    (O_WRONLY | O_CREAT | O_EXCL))),
                      0644)) < 0) {
           fprintf(stderr,
                   "%s: couldn't open destination %s for standard output\n",
                   argv[0], stdoutfile);
           exit(2);
       }
       if (dup2(fd, 1) < 0) {
           fprintf(stderr,
                   "%s: couldn't open destination %s for standard output\n",
                   argv[0], stdoutfile);
           exit(2);
       }
       close(fd);
    }
    if (rstderr) {
       if ((fd = open(stderrfile, (astderr ?
                                   (O_WRONLY | O_APPEND | O_CREAT) :
                                   (cstderr ?
                                    (O_WRONLY | O_CREAT | O_TRUNC) :
                                    (O_WRONLY | O_CREAT | O_EXCL))),
                      0644)) < 0) {
           fprintf(stderr,
                   "%s: couldn't open destination %s for standard error\n",
                   argv[0], stderrfile);
           exit(2);
       }
       if (dup2(fd, 2) < 0) {
           fprintf(stderr,
                   "%s: couldn't open destination %s for standard error\n",
                   argv[0], stderrfile);
           exit(2);
       }
       close(fd);
    }

    if (argv[optind] == NULL) {
	fprintf(stderr, USAGESTRING, argv[0]);
        exit(1);
    }
    execvp(argv[optind], argv + optind);
    fprintf(stderr, "%s: ", argv[0]);
    perror(argv[optind]);
}
