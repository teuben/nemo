#include <stdio.h>
#include <math.h>
#include <suntool/sunview.h>
#include <suntool/canvas.h>
#include <suntool/panel.h>
#include "vogl.h"
#include "vodevice.h"

#define SIZE	512

Object	thing = 2, filledthing = 1, outlinething = 2;

int	back = 0;
int	doubleb = 1;
int	fill = 0;
int	hatch = 0;

#define TRANS           25.0
#define SCAL            0.1
#define FILLED          2
#define OUTLINE         3

float   tdir = TRANS;
float   scal = 1.0 + SCAL;
int     but, nplanes;

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
		singlebuffer(1);
}

quit()
{
	gexit();
	exit(0);
}

resize(canvas, w, h)
	Canvas	canvas;
	int	w, h;
{

	fprintf(stderr, "Resize proc\n");
	vo_sunview_canvas_size(w, h);
	reshapeviewport();
}

main(ac, av)
	int	ac;
	char	**av;
{
	Frame		frame;
	Canvas		canvas;
	Panel		panel;
	int		w, h;
	Notify_value	drawscene();

	frame = window_create(
			0, FRAME,
			FRAME_LABEL, av[1],
		0);

	canvas = window_create(
			frame, CANVAS,
			CANVAS_AUTO_EXPAND, TRUE, 
			CANVAS_AUTO_SHRINK, TRUE,
			CANVAS_RESIZE_PROC, resize,
			WIN_HEIGHT, SIZE,
			WIN_WIDTH, SIZE,
		0);



	panel = window_create(
			frame, PANEL, 
			WIN_BELOW, canvas,
			WIN_X, 0,
		0);


	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_IMAGE, panel_button_image(panel, "Backface", 0, 0),

		PANEL_NOTIFY_PROC, do_back,
	0);
	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_IMAGE, panel_button_image(panel, "Fill", 0, 0),

		PANEL_NOTIFY_PROC, do_fill,
	0);

	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_IMAGE, panel_button_image(panel, "Scale", 0, 0),

		PANEL_NOTIFY_PROC, do_scale,
	0);

	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_IMAGE, panel_button_image(panel, "X-trans", 0, 0),

		PANEL_NOTIFY_PROC, do_x,
	0);

	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_IMAGE, panel_button_image(panel, "Y-trans", 0, 0),

		PANEL_NOTIFY_PROC, do_y,
	0);

	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_IMAGE, panel_button_image(panel, "Z-trans", 0, 0),

		PANEL_NOTIFY_PROC, do_z,
	0);

	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_IMAGE, panel_button_image(panel, "Double buffer", 0, 0),

		PANEL_NOTIFY_PROC, do_double,
	0);

	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_IMAGE, panel_button_image(panel, "QUIT", 0, 0),

		PANEL_NOTIFY_PROC, quit,
	0);

	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_IMAGE, panel_button_image(panel, " + ", 0, 0),

		PANEL_NOTIFY_PROC, do_plus,
	0);
	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_IMAGE, panel_button_image(panel, " - ", 0, 0),

		PANEL_NOTIFY_PROC, do_minus,
	0);

	window_fit(panel);
	window_fit(frame);
	w = (int)window_get(canvas, WIN_WIDTH);
	h = (int)window_get(canvas, WIN_HEIGHT);
	vo_sunview_canvas(canvas, w, h);
	ginit();

	setup_lcube();

#undef USE_TIMER
#ifdef USE_TIMER
/* The following sets the timer */

	/* FAST AS POSSIBLE */

	notify_set_itimer_func(frame,
		(Notify_func)drawscene, ITIMER_REAL, &NOTIFY_POLLING_ITIMER, NULL);

	window_main_loop(frame);

#else
	window_set(frame, WIN_SHOW, TRUE, NULL);

	while(1) {
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

	if (doubleb)
		swapbuffers();

	return NOTIFY_DONE;
}

