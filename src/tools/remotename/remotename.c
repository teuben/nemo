/*
 * remotehost		W.Sebok Nov.3, 1988
 * 	modified May 10, 1990 to handle possible truncation of local domain
 *	qualified name to 16 characters
 *
 * try to return name of host the user is sitting in front of:
 *  if remotely logged in from some host print name of that host
 *  if not remotely logged in print name of this host
 *
 * can be fooled if one remotely logs into host then from there remotely logs
 * into a second host
 *
 * may fail if name of remote host is larger than sizeof u.ut_host == 16 char
 */
#include <sys/types.h>
#include <sysexits.h>
#include <utmp.h>
#include <stdio.h>

#define UTMP "/etc/utmp"

typedef enum {
	OKAY,
	TRUNCATED,
	LOCAL,
	BAD
} Machret;

char *index(), *strcpy(), *strncpy();
void exit();
off_t lseek();

static Machret getmach();
static char *getdomain();

main(ac,av)
int ac;
char *av[];
{
	char *p, *d, *dom;
	int n, nh, debug;
	char host[BUFSIZ];
	char thishost[BUFSIZ];

	debug = (ac>1);
	if(debug)fprintf(stderr,"%s in debug mode\n",av[0]);
	
	thishost[sizeof thishost - 1] = '\0';
	(void)gethostname(thishost, sizeof thishost - 1);

	switch (getmach(host)) {
	case OKAY:
		fprintf(stderr,"OKAY\n");
		(void)puts(host);
		exit(0);
		/*NOTREACHED*/
	case LOCAL:
		fprintf(stderr,"LOCAL\n");
		(void)puts(thishost);
		exit(0);
		/*NOTREACHED*/
	case BAD:
		fprintf(stderr,"BAD\n");
		(void)puts(thishost);
		exit(EX_OSERR);
		/*NOTREACHED*/
	}

	fprintf(stderr,"Allright, let's go\n");
	/* truncated */
	p = host;
	dom = getdomain(thishost);
	nh = strlen(host);

	/* see if it is somewhere in the local domain */
	while ((p = index(p, '.')) != NULL) {
		p++;
		if ((n = nh - (p - host)) <= 0)
			break;
		for (d=dom-1; d != NULL; d = index(d, '.')) {
			d++;
			if (strncmp(p, d, n) == 0) {
				(void)strcpy(p, d);
				(void)puts(host);
				exit(0);
				/*NOTREACHED*/
			}
		}
	}
	(void)puts(host);
	exit(EX_NOHOST);
	/*NOTREACHED*/
}

static Machret
getmach(buf)
	char *buf;
{
	int f, n;
	struct utmp u;
	if ((n = ttyslot()) == 0)
		return(BAD);

	if ((f = open(UTMP, 0)) < 0)
		return(BAD);

	if (lseek(f, (off_t)(n*sizeof(struct utmp)), 0) < 0 ||
	    read(f, (char *)&u, sizeof(struct utmp)) != sizeof (struct utmp)) {
		(void)close(f);
		return(BAD);
	}

	(void)close(f);

	if (u.ut_host[0] == '\0')
		return(LOCAL);

	(void)strncpy(buf, u.ut_host, sizeof u.ut_host);
	buf[sizeof u.ut_host] = '\0';

	return((strlen(buf) == sizeof u.ut_host)? TRUNCATED: OKAY);
}

static char *
getdomain(thishost)
	char *thishost;
{
	register char *d;
#ifndef DOMAIN
	static char domain[66];
#endif
	/* see if domain is contained within host name */
	if ((d = index(thishost, '.')) != NULL) {
		d++;
		return(d);
	}

#ifdef DOMAIN
	return(DOMAIN);
#else
	/* use sun's getdomainname() function */
	if (getdomainname(domain, sizeof domain) < 0)
		return("");

	if ((d = index(domain,'.')) == NULL)
		return(domain);
	d++;
	return(d);
#endif
}
