/*
 * PLARROW.C: draw a arrow.
 */

#include <stdinc.h>
#include <getparam.h>

string defv[] = {
    "x=5.0",
    "y=5.0",
    "len=0.5",
    "wid=0.2",
    "theta=45.0",
    "lwidth=1",
    "VERSION=1.0",
    NULL,
};

main(argc, argv)
int argc;
string argv[];
{
    initparam(argv, defv);
    plinit("", 0.0, 20.0, 0.0, 20.0);
    plltype(getiparam("lwidth"), 0);
    plarrow(getdparam("x"), getdparam("y"),
	    getdparam("len"), getdparam("wid"),
	    getdparam("theta"));
    plstop();
}

plarrow(x, y, l, w, theta)
real x, y;
real l, w;
real theta;
{
    real sin(), cos(), s, c, x1, y1;

    s = sin(PI * theta / 180.0);
    c = cos(PI * theta / 180.0);
    x1 = x - 0.5 * l * c;
    y1 = y - 0.5 * l * s;
    plmove(x1 - 0.5 * w * s, y1 + 0.5 * w * c);
    plline(x, y);
    plline(x1 + 0.5 * w * s, y1 - 0.5 * w * c);
    plmove(x, y);
    plline(x - l * c, y - l * s);
}
