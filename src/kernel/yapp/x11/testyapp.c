/*
 *  TESTYAPP:  TESTBED code for any of the yappXXX.c files
 */

#include <stdinc.h>
#include <getparam.h>
#include <yapp.h>

string defv[] = {
    "name=***\n     Override for yapp= system keyword (*** = default)",
    "screendump=\n  Filename if screendump requested",
    "headline=yapp\nHeadline for plot",
    "pages=1\n      Number of pages to test",
    "VERSION=1.1\n  26-oct-90 PJT",
    NULL,
};

nemo_main()
{
    int i, j, np;
    string name, dumpfile;

    name = getparam("name");
    dumpfile = getparam("screendump");
    np = getiparam("pages");
    printf("Testing wth pages=%d\n",np);
    plinit(name, 0.0, 20.0, 0.0, 20.0);     /* open device */
    x_init_plobj();
    plmove(0.0, 0.0);
    plcolor(0);
    plline(20.0, 0.0);
    plcolor(1);
    plline(20.0, 20.0);
    plcolor(2);
    plline(0.0, 20.0);
    plcolor(3);
    plline(0.0, 0.0);
    plcolor(4);
    plline(20.0, 20.0);
    plmove(20.0, 0.0);
    plcolor(5);
    plline(0.0, 20.0);
    plltype(12, 0);
    plmove(4.0, 18.0);
    plcolor(6);
    plline(16.0, 18.0);
    plltype(-6, 0);
    plmove(6.0, 18.0);
    plcolor(7);
    plline(14.0, 18.0);
    
    for (i = 1; i <= 4; i++) {
	plcolor(8+i);
	plltype(i, 1);
        plmove(1.0, 13.0 - i);
        plline(3.0, 13.0 - i);
        plpoint(3.5, 13.0 - i);
	plltype(1, i);
	for (j = 1; j <= 4; j++) {
	    plmove(1.5, 13.0 - i - 0.2*j);
	    plline(1.5 + j, 13.0 - i - 0.2*j);
	}
    }
    plcolor(12);
    plltype(1, 1);
    plcircle(15.0, 9.0, -0.5);
    plcolor(13);
    plcircle(16.0, 9.0, 0.25);
    plcolor(14);
    plcircle(17.0, 9.0, 0.125);
    plcolor(15);
    plcircle(18.0, 9.0, 0.0625);
    plbox(16.0, 8.0, 0.4);
    plbox(17.0, 8.0, 0.2);
    plbox(18.0, 8.0, -0.2);
    plcross(16.0, 7.0, 0.4);
    plcross(17.0, 7.0, 0.2);
    plcross(18.0, 7.0, -0.2);
    plcolor(4);
    pltext("Foo Bar!", 8.0, 5.0, 0.5, 0.0);
    plcolor(5);
    pltext("Fum Bar!", 8.0, 3.0, 0.25, 0.0);
    plcolor(6);
    for (i = 0; i <= 4; i++)
	pltext(" testing angles", 16.0, 10.0, 0.2, 45.0*i);
    plmove(10.0, 8.5);
    plline(10.0, 11.5);
    pljust(-1);
    plcolor(3);
    pltext("left justified",  10.0,  9.0, 0.25, 0.0);
    plcolor(2);
    pljust(0);
    pltext("centered",        10.0, 10.0, 0.25, 0.0);
    plcolor(1);
    pljust(1);
    pltext("right justified", 10.0, 11.0, 0.25, 0.0);
    pljust(0);
    plcolor(7);
    pltext(getparam("headline"),10.0, 19.0, 0.5, 0.0);
    plcolor(1);
    plflush();
    if (*dumpfile)
        pl_screendump(dumpfile);
    if (np>1) {
        plflush();
        plframe();
        plmove(0.0, 0.0);
        plline(20.0, 0.0);
        plline(20.0, 20.0);
        plline(0.0, 20.0);
        plline(0.0, 0.0);

        pljust(0);
        pltext("This is page 2", 10.0,10.0,0.25,0.0);

#define IMAX 100
#define ISTEP 20.0/IMAX

        plmove (0.0,0.0);
        for (i=0; i<IMAX; i++)
           plline(i*ISTEP, i*ISTEP);
    }
    plstop();
}



