/*
 * This program expects args mmddyy dir1 dir2 ....
 *
 * For each directory named, all files under that directory recursively
 * that are newer than the specified date will be listed on stdout, but
 * it does not check subdirectories that are not on the same device as
 * their parent directory.
 *
 * It should be useful if you want to find out what files on your
 * root and usr devices have been changed since your last upgrade without
 * listing every user file, unless you keep users on /usr.
 *
 * Warren Burstein (pluto!warren) 7/22/87
 * Peter Teuben			  8/13/87  small mods for SUN
 */

#include <stdio.h>
#include <errno.h>
#ifdef SUN
#   include <sys/time.h>
#else
#   include <time.h>
#endif
#include <sys/types.h>
#include <sys/dir.h>
#include <sys/stat.h>
#include <sys/param.h>

char *progname;

main(argc, argv)
char **argv;
{
	int arg;
	long date, parse_date();
	char cwd[MAXPATHLEN];

	progname = argv[0];

	if (argc < 3) {
		fprintf(stderr, "Lists changed files since specified date in directories and below\n");
		fprintf(stderr, "usage: %s mmddyy dir...\n", progname);
		exit(1);
	}

	date = parse_date(argv[1]);

	getwd(cwd);
	for (arg = 2; arg < argc; arg++) {
		/*
		 * look_at may change working directory, so come back
		 * to where we started if this directory arg doesn't start
		 * at /.  Don't bother if this is the first time thru here.
		 */
		if (arg != 2 && *argv[arg] != '/' && chdir(cwd) < 0) {
			fprintf(stderr, "%s: cannot return to %s: ", progname, cwd);
			perror("");
			exit(1);
		}
		look_at(argv[arg], date);
	}

	exit(0);
}

/*
 * Turn a string that looks like mmddyy into a date in the format
 * returned by time(2), seconds since 00:00:00 GMT Jan 1, 1970.
 */

static int months[13] = {
	0,				/* let January be 1, not zero */
    31, 28, 31, 30, 31, 30,
	31, 31, 30, 31, 30, 31
};

long parse_date(mmddyy)
char *mmddyy;
{
	int month, day, year;
	int m;
	long date;
	struct timeval tv;
	struct timezone tz;

	(void) gettimeofday(&tv, &tz);

	if (sscanf(mmddyy, "%2d%2d%2d", &month, &day, &year) != 3) {
		fprintf(stderr, "%s: date must be in mmddyy format\n", progname);
		exit(1);
	}

	if (year < 70) {
		fprintf(stderr, "%s: the world was created in 1970!\n", progname);
		exit(1);
	}

	/*
	 * Adjust February
	 */
	months[2] = is_leap(year + 1900) ? 29 : 28;

	if (month < 1 || month > 12) {
		fprintf(stderr, "%s: month must be between 1 and 12\n", progname);
		exit(1);
	}

	if (day < 1 || day > months[month]) {
		fprintf(stderr, "%s: day must be between 1 and %d\n",
			   progname, months[month]);
		exit(1);
	}

	/*
	 * Find days between Jan 1, 1970 and Jan 1 of year (ignore Julian
	 * corrections here), add days in months and day in month.
	 */
	date = (year - 70) * 365 + year / 4 - 70 / 4;
	for (m = 1; m < month; m++)
	  date += months[m];
	date += day - 1;

	/*
	 * Turn into seconds, correct for longitude and daylight savings time.
	 */
	date *= 24 * 60 *60;
	date += tz.tz_minuteswest * 60;
	if (localtime(&date)->tm_isdst)
	  date -= 60 * 60;

	return date;
}

/*
 * Return true if year is a leap year.
 */
int is_leap(year) {
  return year % 4 == 0 && year % 100 != 0 || year % 400 == 0;
}

/*
 * List all entries in dir that are newer than date.  Recurse on
 * subdirectories that are on the same device as dir.  Treat symbolic
 * links like regular files - list them if the link is new, not the file
 * referenced.
 *
 * Exit on any error since it's a pain to figure out what the current
 * directory is afterwards.
 */

look_at(dir, date)
char *dir;
long date;
{
	extern int errno;
	DIR *dirp;
	struct direct *dp;
	struct stat st;
	dev_t dev;

	if (stat(dir, &st) < 0) {
		fprintf(stderr, "%s: cannot stat directory %s\n", progname, dir);
		fprintf(stderr, "Errno=%d\n",errno);
		exit(1);
	}
	dev = st.st_dev;

	if ( (dirp = opendir(dir)) == NULL) {
		fprintf(stderr, "%s: cannot open directory %s\n", progname, dir);
		exit(1);
	}

	if (chdir(dir) < 0) {
		fprintf(stderr, "%s: cannot chdir to %s\n", progname, dir);
		exit(1);
	}

	while ( (dp = readdir(dirp)) != NULL) {
		/*
		 * Skip self and parent
		 */
		if (strcmp(dp->d_name, ".") == 0 || strcmp(dp->d_name, "..") == 0)
		  continue;

		/*
		 * Use lstat, so as not to follow symbolic links.
		 */
		if (lstat(dp->d_name, &st) < 0) {
			fprintf(stderr, "cannot stat file %s/%s\n", dir, dp->d_name);
			continue;
		}

		if ((st.st_mode & S_IFMT) == S_IFDIR) {
			if (st.st_dev == dev) {
				char dirname[MAXPATHLEN];

				/*
				 * look_at chdir's to a subdirectory, go back up
				 * when done.
				 */
				(void) sprintf(dirname, "%s/%s", dir, dp->d_name);
				look_at(dirname, date);
				(void) chdir("..");
			}
		} else if (st.st_mtime > date)
		  printf("%s/%s\n", dir, dp->d_name);
	}
	(void) closedir(dirp);
}
