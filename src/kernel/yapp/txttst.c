#include <stdinc.h>
#include <getparam.h>

string defv[] = {
    "msg=FOO",
    "sizes=0.24,0.32,0.4,0.5,0.64",
    "angle=0.0",
    NULL,
};

main(argc, argv)
int argc;
string argv[];
{
    string msg, *burststring(), *szptr;
    real angle, ybase, atof();

    initparam(argv, defv);
    msg = getparam("msg");
    szptr = burststring(getparam("sizes"), ", ");
    angle = getdparam("angle");
    plinit("ps-", 0.0, 20.0, 0.0, 20.0);
    ybase = 17.5;
    while (*szptr != NULL) {
	pljust(-1);
	pltext(*szptr, 0.5, ybase, 0.24, 0.0);
	txttst(msg, ybase, atof(*szptr++), angle);
	ybase = ybase - 2.5;
    }
    plstop();
}

txttst(msg, ybase, size, angle)
string msg;
real ybase, size, angle;
{
    plmove( 2.0, ybase);
    plline(18.0, ybase);
    plmove( 5.0, ybase-size/2);
    plline( 5.0, ybase+size/2);
    pljust(-1);
    pltext(msg, 5.0, ybase, size, angle);
    plmove(10.0, ybase-size/2);
    plline(10.0, ybase+size/2);
    pljust(0);
    pltext(msg, 10.0, ybase, size, angle);
    plmove(15.0, ybase-size/2);
    plline(15.0, ybase+size/2);
    pljust(1);
    pltext(msg, 15.0, ybase, size, angle);
}
