
#include <stdio.h>
#include <math.h>
#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Xaw/Form.h>
#include <X11/Xaw/Command.h>
#include "vodevice.h"

long	thing = 2, filledthing = 1, outlinething = 2;

Widget	canvas;

#define SIZE	512


int	iii;
int	back = 0;
int	doubleb = 1;
int	fill = 0;
int	hatch = 0;

int	vinited = 0;

#define TRANS           25.0
#define SCAL            0.1
#define FILLED          2
#define OUTLINE         3

float   tdir = TRANS;
float   scal = 1.0 + SCAL;
int     but, nplanes;

Display	*dpy;
Window	win;

do_plus()
{
	tdir = TRANS;
}

do_minus()
{
	tdir = -tdir;

	if (scal < 1.0)
		scal = 1.0 + SCAL;
	else
		scal = 1.0 - SCAL;
}

do_back()
{
	back = !back;
	backface(back);
}

do_fill()
{
	fill = !fill;

	thing = (fill ? filledthing : outlinething);
}

do_scale()
{
	scale(scal, scal, scal);
}

do_x()
{
	translate(tdir, 0.0, 0.0);
}
do_y()
{
	translate(0.0, tdir, 0.0);
}

do_z()
{
	translate(0.0, 0.0, tdir);
}

do_double()
{
	doubleb = !doubleb;

	if (doubleb) {
		fprintf(stderr, "Double buffer on\n");
		doublebuffer(1);
	} else
		singlebuffer();
}

quit()
{
	gexit();
	exit(0);
}

resize()
{
	Dimension       w, h;
	Arg             arg[2];

	XtSetArg(arg[0], XtNwidth, &w);
	XtSetArg(arg[1], XtNheight, &h);
	XtGetValues(canvas, arg, 2);

	fprintf(stderr, "resize() %d %d\n", w, h);
	vo_xt_win_size((int)w, (int)h);
	reshapeviewport();
}

repaint()
{
	fprintf(stderr, "repaint() called\n");
	color(0);
	clear();
	
}

XtActionsRec actions[] = {
	{"repaint", 	(XtActionProc)repaint},
	{"resize", 	(XtActionProc)resize}
};

String trans =
	"<Expose>:	repaint()\n \
	 <Configure>:	resize()";


Display		*dpy;
Window		win;
GC		gc;

/*
 * simple program to display a polygon file
 */
main(ac, av)
	int	ac;
	char	**av;
{
	int		w, h;
	Widget		toplevel, 
			panel,
			qbut,
			bbut,
			fbut,
			dbut,
			xbut,
			ybut,
			zbut,
			sbut,
			pbut,
			mbut;

	Arg		wargs[5];
	void		drawscene();
	Dimension	ww, wh, x, y;
	XtTranslations	trans_table;

	vinited = 0;

	ww = wh = 100;
	x = 0;
	y = 0;
	toplevel = XtInitialize(av[0], "xtlcube", NULL, 0, &ac, av);

	panel = XtCreateManagedWidget("panel",
			formWidgetClass,
			toplevel,
			NULL,
		0);


	XtSetArg(wargs[0], XtNlabel, "Quit");
	qbut = XtCreateManagedWidget("quit", 
			commandWidgetClass,
			panel,
			wargs,
		1);

	XtAddCallback(qbut, XtNcallback, quit, NULL);

	XtSetArg(wargs[0], XtNlabel, "Backface");
	XtSetArg(wargs[1], XtNfromHoriz, qbut);
	bbut = XtCreateManagedWidget("backface",
			commandWidgetClass,
			panel,
			wargs,
		2);
	XtAddCallback(bbut, XtNcallback, do_back, NULL);

	XtSetArg(wargs[0], XtNlabel, "Fill");
	XtSetArg(wargs[1], XtNfromHoriz, bbut);
	fbut = XtCreateManagedWidget("fill",
			commandWidgetClass,
			panel,
			wargs,
		2);
	XtAddCallback(fbut, XtNcallback, do_fill, NULL);

	XtSetArg(wargs[0], XtNlabel, "DoubleBuffer");
	XtSetArg(wargs[1], XtNfromHoriz, fbut);
	dbut = XtCreateManagedWidget("doubleb",
			commandWidgetClass,
			panel,
			wargs,
		2);
	XtAddCallback(dbut, XtNcallback, do_double, NULL);

	XtSetArg(wargs[0], XtNwidth, 512);
	XtSetArg(wargs[1], XtNheight, 512);
	XtSetArg(wargs[2], XtNfromVert, qbut);
	canvas = XtCreateManagedWidget("canvas", 
			simpleWidgetClass,
			panel,
			wargs,
		3);

	XtAddActions(actions, XtNumber(actions));
	trans_table = XtParseTranslationTable(trans);
	XtAugmentTranslations(canvas, trans_table);

	XtRealizeWidget(toplevel);


	dpy = XtDisplay(canvas);
	win = XtWindow(canvas);

	vo_xt_window(dpy, win, 512, 512);
	ginit();
	vinited = 1;
	setup_lcube();

	while(1) {
		XEvent	event;

		while(XtPending()) {
			XtNextEvent(&event);
			XtDispatchEvent(&event);
		}
		drawscene();

	}
}

void
drawscene()
{

	short   x, y;

	x = 500 - (int)getvaluator(MOUSEX);
	y = 500 - (int)getvaluator(MOUSEY);
	x *= 3;
	y *= 3;
	pushmatrix();
		rotate(x, 'y');
		rotate(y, 'x');
		color(0);
		clear();
		callobj(thing); 
	popmatrix();

	if (doubleb)
		swapbuffers();

}

