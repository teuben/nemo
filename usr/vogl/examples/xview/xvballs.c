#include <stdio.h>
#include <math.h>
#include <xview/xview.h>
#include <xview/canvas.h>
#include <xview/panel.h>
#include <xview/xv_xrect.h>
#include "vogl.h"

#define SIZE	512

int	vinited = 0;

Display	*dpy;
Window	win;

quit()
{
	gexit();
	exit(0);
}


repaint(canvas, paint_win, dpy, xwin, area)
	Canvas		canvas;
	Xv_Window	paint_win;
	Display		*dpy;
	Window		xwin;
	Xv_xrectlist	*area;
{
	int	w, h;
	
	w = (int)xv_get(paint_win, XV_WIDTH);
        h = (int)xv_get(paint_win, XV_HEIGHT);

	if (!vinited) {
		vinited = 1;
		vo_xt_window(dpy, xwin, w, h);
		ginit();
	}
	
	fprintf(stderr, "Repaint proc: 0x%x, 0x%x [%d %d]\n", dpy, win, w, h);
	color(BLACK);
	clear();

	draw_balls();
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

main(ac, av)
	int	ac;
	char	**av;
{
	Frame		frame;
	Canvas		canvas;
	Panel		panel;
	int		w, h;

	frame = xv_create(
			0, FRAME,
			FRAME_LABEL, av[1],
			WIN_HEIGHT, SIZE,
			WIN_WIDTH, SIZE,
		0);





	canvas = xv_create(
			frame, CANVAS,
			CANVAS_RESIZE_PROC, resize,
			CANVAS_REPAINT_PROC, repaint,
			CANVAS_X_PAINT_WINDOW, TRUE,
			WIN_HEIGHT, SIZE,
			WIN_WIDTH, SIZE,
		0);


	panel = xv_create(
			frame, PANEL, 
			WIN_X, 0,
			WIN_BELOW, canvas,
		0);

	panel_create_item(
		panel, PANEL_BUTTON,
		PANEL_LABEL_STRING, "QUIT",

		PANEL_NOTIFY_PROC, quit,
	0);

	window_fit(panel);
	window_fit(frame);

	/*
	 * The resiz/repaint procs are going to initialise
 	 * it all for us.
	 */

	xv_main_loop(frame);
}
