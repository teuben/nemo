#include <stdio.h>
#include <math.h>
#include <xview/xview.h>
#include <xview/canvas.h>
#include <xview/panel.h>
#include "vogl.h"
#include "vodevice.h"

#define SIZE	512

unsigned	flags = 1;
int	vinited = 0;

#define DOUBLEB		1
#define BACKFACE	2
#define FILL		4

#define TRANS           25.0
#define SCAL            0.1
#define FILLED          2
#define OUTLINE         3

Object	thing = 2, filledthing = 1, outlinething = 2;

float   tdir = TRANS;
float   scal = 1.0 + SCAL;
int     but, nplanes;
int	w, h;

Display	*dpy;
Window	win, win1, win2;

Panel_item 	toggle;

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
	backface(flags & BACKFACE);
}

do_fill()
{
	if (flags & FILL)
		thing = filledthing;
	else
		thing = outlinething;

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
	if (flags & DOUBLEB) {
		doublebuffer(1);
	} else {
		singlebuffer();
	}
}

quit()
{
	gexit();
	exit(0);
}


resize(win, event, arg, type)
	Xv_window       win;
	Event           *event;
	Notify_arg      arg;
	Notify_event_type type;
{
	int	w, h;
	
	w = xv_get(win, XV_WIDTH);
        h = xv_get(win, XV_HEIGHT);

	fprintf(stderr, "Resize proc: 0x%x, 0x%x (%d %d)\n", dpy, win, w, h);

        vo_xt_win_size(w, h);
	reshapeviewport();

}

menu_proc(menu, menu_item)
	Menu		menu;
	Menu_item	menu_item;
{
	char	*choice = (char *)xv_get(menu_item, MENU_STRING);

	if (!strcmp(choice, "DoubleBuffer")) {
		if (flags & DOUBLEB)
			flags &= ~DOUBLEB;
		else
			flags |= DOUBLEB;

		do_double();
	} else if (!strcmp(choice, "Backface")) {
		if (flags & BACKFACE)
			flags &= ~BACKFACE;
		else
			flags |= BACKFACE;

		do_back();
	} else if (!strcmp(choice, "Filled")) {
		if (flags & FILLED)
			flags &= ~FILL;
		else
			flags |= FILL;

		do_fill();
	} else if (!strcmp(choice, "QUIT"))
		quit();

	xv_set(toggle, PANEL_VALUE, flags, NULL);
}

/*
 * Call menu_show() to display menu on right mouse button push.
 */
void
my_event_proc(window, event)
	Xv_Window window;
	Event *event;
{
	if (event_action(event) == ACTION_MENU && event_is_down(event)) {
		Menu menu = (Menu)xv_get(window, WIN_CLIENT_DATA);
		menu_show(menu, window, event, NULL);
	}
}



toggle_selected(item, value, event)
	Panel_item item;
	unsigned value;
	Event *event;
{
	char	buf[32];
	Frame	frame = xv_get(item, PANEL_CLIENT_DATA);

	fprintf(stderr, "Toggle(%d) ... \n", value);

	buf[0] = 0;
	if (event_id(event) == MS_LEFT) {
		flags = value;
		do_fill();
		do_back();
		do_double();

		xv_set(item, PANEL_VALUE, flags, NULL);

		return XV_OK;
	}
	return XV_ERROR;
}

do_other()
{
	if (win == win1)
		win = win2;
	else
		win = win1;

	vo_xt_window(dpy, win, w, h);
}

main(ac, av)
	int	ac;
	char	**av;
	{
	Frame		frame;
	Canvas		canvas1, canvas2;
	Panel		panel;
	Menu		menu;
	Notify_value	drawscene();

	frame = xv_create(
		0, FRAME,
		FRAME_LABEL, av[1],
		WIN_HEIGHT, SIZE,
		WIN_WIDTH, SIZE,
	0);

	canvas1 = xv_create(
		frame, CANVAS,
		CANVAS_RESIZE_PROC, resize,
		WIN_HEIGHT, SIZE,
		WIN_WIDTH, SIZE,
	0);

	canvas2 = xv_create(
		frame, CANVAS,
		CANVAS_RESIZE_PROC, resize,
		WIN_HEIGHT, SIZE,
		WIN_WIDTH, SIZE,
	0);

	menu = xv_create(0, MENU, 
		MENU_TITLE_ITEM, "Gunge",
		MENU_STRINGS, "DoubleBuffer",
			      "Backface",
			      "Filled",
			      "QUIT",
			      NULL,
		MENU_NOTIFY_PROC, menu_proc,
		NULL
	);


	xv_set(canvas_paint_window(canvas1),
		WIN_CONSUME_EVENTS,  WIN_MOUSE_BUTTONS, NULL,
		WIN_EVENT_PROC,  my_event_proc,
		/* associate the menu to the canvas win for easy retrieval */
		WIN_CLIENT_DATA,  menu,
		NULL
	);

	xv_set(canvas_paint_window(canvas2),
		WIN_CONSUME_EVENTS,  WIN_MOUSE_BUTTONS, NULL,
		WIN_EVENT_PROC,  my_event_proc,
		/* associate the menu to the canvas win for easy retrieval */
		WIN_CLIENT_DATA,  menu,
		NULL
	);

	panel = xv_create(
		frame, PANEL, 
		WIN_BELOW, canvas1,
		WIN_X, 0,
	0);


	toggle = (Panel_item)xv_create(panel, PANEL_TOGGLE,
		PANEL_FEEDBACK,         PANEL_MARKED,
		PANEL_LABEL_STRING,     "Choices",
		PANEL_VALUE,            flags, /* choice 1 */
		PANEL_CHOICE_STRINGS,   "Double Buffer",
					"Backface",
					"Filled",
					 NULL,
		PANEL_NOTIFY_PROC,      toggle_selected,
		PANEL_CLIENT_DATA,      frame,
	0);

	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_STRING, "Scale",

		PANEL_NEXT_ROW, -1, 
		PANEL_NOTIFY_PROC, do_scale,
	0);

	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_STRING, "X-trans",

		PANEL_NOTIFY_PROC, do_x,
	0);

	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_STRING, "Y-trans",

		PANEL_NOTIFY_PROC, do_y,
	0);

	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_STRING, "Z-trans",

		PANEL_NOTIFY_PROC, do_z,
	0);

	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_STRING, "QUIT",

		PANEL_NOTIFY_PROC, quit,
	0);

	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_STRING, " + ",

		PANEL_NOTIFY_PROC, do_plus,
	0 );

	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_STRING, " - ",

		PANEL_NOTIFY_PROC, do_minus,
	0);

	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_STRING, " Other window ",

		PANEL_NOTIFY_PROC, do_other,
	0);

	window_fit(panel);
	window_fit(frame);

	w = (int)window_get(canvas1, WIN_WIDTH);
	h = (int)window_get(canvas1, WIN_HEIGHT);
	win = win1 = (Window)xv_get(canvas_paint_window(canvas1), XV_XID);
	win2 = (Window)xv_get(canvas_paint_window(canvas2), XV_XID);
	dpy = (Display *)xv_get(frame, XV_DISPLAY);

	vo_xt_window(dpy, win, w, h);
	ginit();

	setup_lcube();

	vinited = 1;

#define USE_TIMER 1
#ifdef USE_TIMER
/* The following sets the timer */

	/* FAST AS POSSIBLE */

	notify_set_itimer_func(frame,
		(Notify_func)drawscene, ITIMER_REAL, &NOTIFY_POLLING_ITIMER, NULL);

	xv_main_loop(frame);

#else
	xv_set(frame, WIN_SHOW, TRUE, NULL);

	while(1) {
		XFlush(dpy);
		notify_dispatch();
		drawscene((Notify_client)canvas, 1);
	}
#endif
}

Notify_value
drawscene(c, fd)
	Notify_client c;
	int	fd;
{
	Angle	x, y;

	x = 500 - (int)getvaluator(MOUSEX);
	y = 500 - (int)getvaluator(MOUSEY);
	x *= 3;
	y *= 3;
	pushmatrix();
		rotate(x, 'y');
		rotate(y, 'x');
		color(BLACK);
		clear();
		callobj(thing);	/* The filled or hatched one */
	popmatrix();

	if (flags & DOUBLEB)
		swapbuffers();

	return NOTIFY_DONE;
}

