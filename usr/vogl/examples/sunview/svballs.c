#include <stdio.h>
#include <math.h>
#include <suntool/sunview.h>
#include <suntool/canvas.h>
#include <suntool/panel.h>
#include "vogl.h"

#define SIZE	512

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
	draw_balls();
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
		PANEL_LABEL_IMAGE, panel_button_image(panel, " QUIT ", 0, 0),

		PANEL_NOTIFY_PROC, quit,
	0);

	window_fit(panel);
	window_fit(frame);
	w = (int)window_get(canvas, WIN_WIDTH);
	h = (int)window_get(canvas, WIN_HEIGHT);
	vo_sunview_canvas(canvas, w, h);
	ginit();

	window_main_loop(frame);
}
