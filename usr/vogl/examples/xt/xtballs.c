#include <stdio.h>
#include <math.h>
#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Xaw/Form.h>
#include <X11/Xaw/Command.h>

#define SIZE	512


int	vinited = 0;

int	iii;

Widget		toplevel, panel, canvas, qbut, rbut;

quit()
{
	gexit();
	exit(0);
}

resize()
{
	Dimension	w, h;
	Arg		arg[2];

	XtSetArg(arg[0], XtNwidth, &w);
	XtSetArg(arg[1], XtNheight, &h);
	XtGetValues(canvas, arg, 2);

	vo_xt_win_size((int)w, (int)h);
	reshapeviewport();

	fprintf(stderr, "resize()\n");
}

repaint()
{
	fprintf(stderr, "repaint()\n");
	draw_balls();
	
}

XtActionsRec actions[] = {
	{"resize", 	(XtActionProc)resize},
};

String trans =
	 "<Configure>:	resize()";

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
	Arg		wargs[6];
	void		drawscene();
	Dimension	ww, wh, x, y;
	XtTranslations	trans_table;

	vinited = 0;

	ww = wh = 100;
	x = 0;
	y = 0;
	toplevel = XtInitialize(av[0], "xtballs", NULL, 0, &ac, av);


	panel = XtCreateManagedWidget("panel", formWidgetClass, toplevel, wargs, 1);

	XtSetArg(wargs[0], XtNlabel, "Quit");
	qbut = XtCreateManagedWidget("quit", commandWidgetClass, panel, wargs, 1);
	XtAddCallback(qbut, XtNcallback, quit, NULL);

	XtSetArg(wargs[0], XtNlabel, "Repaint");
	XtSetArg(wargs[1], XtNfromHoriz, qbut);
	rbut = XtCreateManagedWidget("redraw", commandWidgetClass, panel, wargs, 2);
	XtAddCallback(rbut, XtNcallback, repaint, NULL);

	XtAddActions(actions, XtNumber(actions));

	XtSetArg(wargs[0], XtNwidth, 512);
	XtSetArg(wargs[1], XtNheight, 512);
	XtSetArg(wargs[2], XtNx, 0);
	XtSetArg(wargs[3], XtNy, 0);
	XtSetArg(wargs[4], XtNfromVert, qbut);
	trans_table = XtParseTranslationTable(trans);
	XtSetArg(wargs[5], XtNtranslations, trans_table);
	canvas = XtCreateManagedWidget("canvas", simpleWidgetClass, panel, wargs, 6);
	XtAddEventHandler(canvas, ExposureMask, FALSE, repaint, NULL);

	XtRealizeWidget(toplevel);

	dpy = XtDisplay(canvas);
	win = XtWindow(canvas);

	vo_xt_window(dpy, win, 512, 512);
	ginit();
	vinited = 1;

	XtMainLoop();
}
