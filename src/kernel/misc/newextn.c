/*
 * NEWEXTN: append new extension to file name.  Note: this function
 * returns a pointer to a static area -- be careful!
 *	dark-ages	written					JEB
 *	19-oct-90	rindex() -> strrchr()			PJT
 *	25-feb-92	happy gcc2.0				PJT
 */

#include <stdinc.h>

string newextn(string name, string extn)
{
    permanent char buf[64];
    char *tmp;

    strncpy(buf, name, sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = 0;
    tmp = strrchr(buf, '.');
    if (tmp == NULL)
	tmp = buf + strlen(buf);
    strncpy(tmp, extn, sizeof(buf) - 1 - (tmp - buf));
    return buf;
}
