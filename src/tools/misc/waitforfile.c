/*
 * waitforfile.c
 *
 * Wait until a file exists, then continue
 */

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>

extern int getopt();
extern int optind;

#define SECOND (1)
#define MINUTE (60)
#define HOUR (60*60)
#define DAY (24*60*60)


int main(argc, argv)
int argc;
char **argv;
{
	time_t clock;
	struct tm *now;
	struct stat st;
	int c, mode = SECOND;
	char *fname;

	if (argc==2)
   	    fname = argv[1];
	else
	    usage();

        for (;;) {
     	    if (stat(fname,&st) == 0) break;
     	    sleep(5);
     	}
	exit(0);
}

usage()
{
  fprintf(stderr, "usage: %s file\n","age");
  fprintf(stderr, "return when file exists, otherwise wait until it does\n");
  exit (1);
}			

