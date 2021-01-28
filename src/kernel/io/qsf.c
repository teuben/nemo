/*
 * QSF: check if a file is a structured file
 *
 *	Note: this program returns a status to the parent shell
 *		0 if OK, 1 if no structured file
 *
 *	V1.0  13-feb-92	Created				PJT
 *       1.0a 15-may-92 fixed bug: isatty(fileno(str))  PJT
 *	    b  5-mar-94 ansi
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>


string defv[] = {               /* DEFAULT INPUT PARAMETERS */
    "in=???\n                     input file name to test",
    "VERSION=1.0c\n		  6-jun-08 PJT ",
    NULL,
};

string usage = "check if a file is a structured file";

string cvsid="$Id$";

void nemo_main(void)
{
    string in = getparam("in");
    stream str= stropen(in,"r");

    dprintf(2,"Fileno(\"%s\") = %d\n",in,fileno(str));
    if (isatty(fileno(str))) {
        dprintf(1,"File %s is not a proper file (isatty)\n",in);
        stop(1);
    } else if (qsf(str)) {
        dprintf(1,"File %s is a proper binary structured file\n",in);
        stop(0);
    } else {
        dprintf(1,"File %s is not a proper binary structured file\n",in);
        stop(1);
    }
}
