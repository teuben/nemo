/*
 * PLTEXT: simple routine to call pltext with specified parameters.
 * Mostly to be used as pltext_ps to generate figure labels, etc.
 */

#include <stdinc.h>
#include <getparam.h>

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "msg=Hello World",		/* message to plot */
    "x=10.0",			/* x coordinate of message */
    "y=10.0",			/* y coordinate of message */
    "hgt=0.4",			/* height of characters */
    "ang=0.0",			/* angle, anti-clockwise in deg */
    "just=-1",			/* justification: -1, 0, 1 */
    "VERSION=1.0",		/* JEB  17 October 1987  IAS */
    NULL,
};

main(argc, argv)
int argc;
string argv[];
{
    string msg;
    real x, y, hgt, ang;
    int just;

    initparam(argv, defv);
    msg = getparam("msg");
    x = getdparam("x");
    y = getdparam("y");
    hgt = getdparam("hgt");
    ang = getdparam("ang");
    just = getiparam("just");
    plinit(streq(getargv0(), "pltext_ps") ? "ps-" : "***",
	   0.0, 20.0, 0.0, 20.0);
    pljust(just);
    pltext(msg, x, y, hgt, ang);
    plstop();
}
