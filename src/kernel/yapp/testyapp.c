/*
 *  TESTYAPP:  TESTBED code for any of the yappXXX.c files
 *	This has to be linked with any of the yapp_*.c
 *	files to an executable. The Makefile contains
 *	a direct target for most device drivers...
 *
 *  V1.1 26-oct-90      
 *  V1.2 22-may-92      color  (tested for pgplot)      pjt
 *			-DCOLOR must be used, as not all yapp's have
 *			stubbs for the color functions
 *	 5-may-95       allow > 2 pages..
 */

#include <stdinc.h>
#include <getparam.h>
#include <yapp.h>

string defv[] = {
    "name=***\n     Override for yapp= system keyword (*** = default)",
    "screendump=\n  Filename if screendump requested",
    "headline=yapp\nHeadline for plot",
    "pages=1\n      Number of pages to test",
    "lut=\n         Test this new color LUT (an ascii table)",
    "VERSION=1.3\n  5-may-95 PJT",
    NULL,
};

nemo_main()
{
    int i, j, ip, np, nc;
    string name, dumpfile, headline;
    char label[80];

    name = getparam("name");
    dumpfile = getparam("screendump");
    np = getiparam("pages");
    headline = getparam("headline");
    printf("Testing wth pages=%d  headline=%s\n",np,headline);

    plinit(name, 0.0, 20.0, 0.0, 20.0);     /* open device */
    /* plframe();  */		/* testing */
    plmove(0.0, 0.0);
    plline(20.0, 0.0);
    plline(20.0, 20.0);
    plline(0.0, 20.0);
    plline(0.0, 0.0);
    plline(20.0, 20.0);
    plmove(20.0, 0.0);
    plline(0.0, 20.0);
    plltype(12, 0);
    plmove(4.0, 18.0);
    plline(16.0, 18.0);
    plltype(-6, 0);
    plmove(6.0, 18.0);
    plline(14.0, 18.0);
#ifdef COLOR
    if (hasvalue("lut")) pllut(getparam("lut"),TRUE);	/* load compressed LUT */
    plltype(10, 0);
    nc = plncolors();
    printf("Found %d colors\n",nc);
    for (i = 1; i < nc; i++) {
	plcolor(i);
	plmove(6.0 + 4 * (i - 1.0) / nc, 16.0);
	plline(6.5 + 4 * (i - 1.0) / nc, 17.0);
    }
    plltype(2, 0);
    for (i = 1; i < nc; i++) {
	plcolor(i);
	plmove(12.0 + 8 * (i - 1.0) / nc, 16.0);
	plline(12.5 + 8 * (i - 1.0) / nc, 17.0);
    }
    plcolor(1); /* select foreground color again */
#endif
    for (i = 1; i <= 4; i++) {
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
    plltype(1, 1);
    plcircle(15.0, 9.0, 0.5);
    plcircle(16.0, 9.0, 0.25);
    plcircle(17.0, 9.0, 0.125);
    plcircle(18.0, 9.0, 0.0625);
    plbox(16.0, 8.0, 0.4);
    plbox(17.0, 8.0, 0.2);
    plbox(18.0, 8.0, -0.2);
    plcross(16.0, 7.0, 0.4);
    plcross(17.0, 7.0, 0.2);
    plcross(18.0, 7.0, -0.2);
    pltext("Foo Bar!", 8.0, 5.0, 0.5, 0.0);
    pltext("Fum Bar!", 8.0, 3.0, 0.25, 0.0);
    for (i = 0; i <= 4; i++)
	pltext(" testing angles", 16.0, 10.0, 0.2, 45.0*i);
    plmove(10.0, 8.5);
    plline(10.0, 11.5);
    pljust(-1);
    pltext("left justified",  10.0,  9.0, 0.25, 0.0);
    pljust(0);
    pltext("centered",        10.0, 10.0, 0.25, 0.0);
    pljust(1);
    pltext("right justified", 10.0, 11.0, 0.25, 0.0);
    pljust(0);
    pltext(headline,10.0, 19.0, 0.5, 0.0);
    if (*dumpfile)
        pl_screendump(dumpfile);
    for (ip=2; ip<=np; ip++) {
	dprintf(0,"Creating page %d\n",ip);
        plflush();
        plframe();
        plmove(0.0, 0.0);
        plline(20.0, 0.0);
        plline(20.0, 20.0);
        plline(0.0, 20.0);
        plline(0.0, 0.0);

        pljust(0);
        sprintf(label,"This is page %d",ip);
        pltext(label, 10.0,10.0,0.25,0.0);

#define IMAX 100
#define ISTEP 20.0/IMAX

        plmove (0.0,0.0);
        for (i=0; i<IMAX; i++)
           plline(i*ISTEP, i*ISTEP);
    }
    plstop();
}


