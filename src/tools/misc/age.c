/*
 * age.c
 *
 * Age of a file
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

	while ((c = getopt(argc, argv, "smhd")) != -1)
		switch (c) {
		case 'd':
			mode = DAY;  break;
		case 'h':
			mode = HOUR; break;
		case 'm':
			mode = MINUTE; break;
		case 's':
			mode = SECOND; break;
		default:
			usage();
		}

	if (optind < argc)
		fname = argv[optind];
	else
		usage();

	time(& clock);
	/* now = localtime(& clock); */
	stat(fname,&st);
	clock -= st.st_mtime;
	now = localtime(&clock);
	
	printf("%ld\n",clock/mode);
	
	exit(0);
}

usage()
{
fprintf(stderr, "usage: %s [-dhms] file\n","age");
fprintf(stderr, "return age of a file in days, hour, minutes or seconds\n");
exit (1);
}			

