#include <stdio.h>

int cdplayer_pid = -1;	/* PID of process that plays the CD */

void cdrom ()
{
#if defined CDPLAYER	/* The CDPLAYER is the CD-ROM player (Sun demo) */
	char buf[200];

	cdplayer_pid = fork();
	if (cdplayer_pid == NULL) {
		sprintf (buf, "%s", CDPLAYER);
		execl (buf, buf, "-n", (char *)0);	/* Exec the CDPLAYER */
	} else if (cdplayer_pid == -1)
		fprintf (stderr, "Couldn't fork CDPLAYER process\n");
#endif
}
