/*
 * Movietool: display a succession of rasterfile frames in real-time
 *
 * Copyright 1989, 1990 by Ole H. Nielsen
 *
 * Author: Ole H. Nielsen
 *         Lab of Applied Physics, Bygn. 307
 *         Technical University of Denmark, DK-2800 Lyngby
 *         E-mail: ohnielse@ltf.dth.dk
 */
#include "movietool.h"
#include <stdio.h>
#include <suntool/canvas.h>
#include <suntool/scrollbar.h>
#include <suntool/panel.h>
#include <suntool/alert.h>
#include <sys/types.h>
#include <sys/timeb.h>

/* global declarations */

extern Frame	window_frame, panel_frame;
extern Window 	panel;
extern Canvas	canvas;
extern Pixwin	*canvaspw;
extern Image	*first_frame;
extern Canvas_par canvas_par;
extern int	num_images;
extern int	identical_colormaps;
extern int	xor_flag;
extern void	my_pw_rop ();
extern int	encoded_flag;
extern int	zoom;
extern int	background_color;
extern int	cdplayer_pid;

char	*title = "Movietool-1.2";	/* Title appearing in tool's top bar */
Image	*current_frame = (Image *)NULL;	/* The image being displayed */
int	playing = FALSE;	/* True if Play button was pressed */
int	clipping = TRUE;	/* False if copying encoded images */
				/* onto framebuffer without clipping */

/* local declarations */

#define MAX_FRAMES_PER_SEC 100
#define FORWARD 0
#define BACKWARD 1

static void	display_frame ();
static int	current_number = 0;	/* Current image's number */
static int	step = 1;		/* Stepping +1 or -1 */
static int	auto_reverse = FALSE;	/* Auto-reverse movie */
static int	repeat = FALSE;		/* Repeat movie */
static int	xor_rop;		/* Use XOR raster-ops */
static int	xor_show = FALSE;	/* Show the XOR'ed images */
static int	Halt = FALSE;		/* True if STOP button was pressed */
static unsigned frames_per_sec = MAX_FRAMES_PER_SEC;
static Panel_item framenumber_item, play_item, step_item, eject_item,
		xor_rop_item, xor_show_item, frame_item, filename_item,
		direction_item, reverse_item, repeat_item, per_sec_item,
		clear_item, clipping_item, repaint_item;

static void filename_proc ()
{
	if (current_frame)
		(void) panel_set (filename_item,
			PANEL_LABEL_STRING, current_frame->filename,
			0);
}

static void direction_proc (item, value, event)
Panel_item	item;
unsigned int	value;
Event		*event;
{
	if (value == BACKWARD)
		step = -1;
	else if (value == FORWARD)
		step = 1;
}

static void xor_rop_proc (item, value, event)
Panel_item	item;
unsigned int	value;
Event		*event;
{
	/* Use XOR RasterOP ? */
	xor_rop = value;
}

static void xor_show_proc (item, value, event)
Panel_item	item;
unsigned int	value;
Event		*event;
{
	xor_show = value;
	if (xor_show)	/* Disable the "XOR_use" button */
		(void) panel_set (xor_rop_item, PANEL_SHOW_ITEM, FALSE, 0);
	else
		(void) panel_set (xor_rop_item, PANEL_SHOW_ITEM, TRUE, 0);
	if (current_frame)
		(void) display_frame (current_frame);
}

static void reverse_proc (item, value, event)
Panel_item	item;
unsigned int	value;
Event		*event;
{
	auto_reverse = value;
	if (auto_reverse)	/* Disable the "Repeat" button */
		(void) panel_set (repeat_item, PANEL_SHOW_ITEM, FALSE, 0);
	else
		(void) panel_set (repeat_item, PANEL_SHOW_ITEM, TRUE, 0);
}

static void repeat_proc (item, value, event)
Panel_item	item;
unsigned int	value;
Event		*event;
{
	repeat = value;
}

static void clipping_proc (item, value, event)
Panel_item	item;
unsigned int	value;
Event		*event;
{
	clipping = value;
}

static Image *locate_frame (number)
int number;	/* Assumed between 1 and num_images */
{
	Image *frame;
	if (current_number == 0)
		(void) next_frame (&current_number);
	if (current_frame == (Image *)NULL)
		frame = first_frame;
	else
		frame = current_frame;

	while (number > frame->frameno && frame->right != (Image *)NULL)
		frame = frame->right;
	while (number < frame->frameno && frame->left != (Image *)NULL)
		frame = frame->left;
	return frame;
}

static void display_frame (frame)
Image *frame;
{
	int xoff, yoff;
	long int op = PIX_SRC;
	struct pixrect *pix, *xor, *pr;

	(void) panel_set_value (framenumber_item, frame->frameno);

	pix = frame->pix;
	if (xor_rop || xor_show) {
		if (xor_rop && playing && step < 0)
			/* When playing backwards, the XOR-pixrect required
			 * is the one of the frame to the right */
			xor = frame->right->xor;
		else
			xor = frame->xor;
		if (xor != PIXRECT_NULL) {
			if (xor_show)
				pix = xor;
			else if (xor_rop && playing) {
				pix = xor;
				op = PIX_DST ^ PIX_SRC;	/* XOR operation */
			}
		}
	}

	if (zoom && frame->ras_type == RT_STANDARD) {
		xoff = (canvas_par.width  - pix->pr_size.x) / 2;
		yoff = (canvas_par.height - zoom * pix->pr_size.y) / 2;
	} else if (zoom && frame->ras_type == RT_BYTE_ENCODED) {
		xoff = (canvas_par.width  - zoom * pix->pr_size.x) / 2;
		yoff = (canvas_par.height - zoom * pix->pr_size.y) / 2;
	} else {
		xoff = (canvas_par.width  - pix->pr_size.x) / 2;
		yoff = (canvas_par.height - pix->pr_size.y) / 2;
	}
	xoff = xoff<0 ? 0 : xoff;
	yoff = yoff<0 ? 0 : yoff;

	if (! clipping && playing && xoff > 0 && yoff > 0)
		op |= PIX_DONTCLIP;		/* Safe to omit clipping */

	/* Set colormap */
	if (identical_colormaps == FALSE && frame->map) {
		(void)pw_putcolormap (canvaspw, 0, CMS_SIZE, frame->map->map[0],
			frame->map->map[1], frame->map->map[2]);
	}
	/* Write pixrect to canvas-pixwin */
	(void) my_pw_rop (canvaspw, xoff, yoff,
			pix->pr_size.x, pix->pr_size.y, op, pix, 0, 0, zoom);
}

static void message (msg)
char *msg;
{
	if (*msg) {
		char text[100];
		(void) sprintf (text, "%s - %s", title, msg);
		(void) window_set (window_frame, FRAME_LABEL, text, 0);
	} else
		(void) window_set (window_frame, FRAME_LABEL, title, 0);
}

void step_proc ()
{
	(void) message ("");
	if (next_frame (&current_number) == FALSE)
		return;
	current_frame = locate_frame (current_number);
	(void) filename_proc ();
	(void) display_frame (current_frame);
}

static void clear_proc ()
{
	register int color;

	if (canvas_par.depth == 1) {		/* Monochrome */
		if (background_color)
			color = PIX_SET;
		else
			color = PIX_CLR;
	} else {				/* Color */
		color = PIX_SRC | PIX_COLOR(background_color);
	}

	/* Need CANVAS_WIDTH, CANVAS_HEIGHT to clear entire canvas */
	(void) pw_writebackground (canvaspw, 0, 0,
		canvas_par.canvas_width, canvas_par.canvas_height, color);
}

static void frame_proc (item, frame, event)
Panel_item item;
unsigned int frame;
struct inputevent *event;
{
	int tmp = step;
	current_number = (int)frame;
	step = 0;
	(void) step_proc();
	step = tmp;
}

static void per_sec_proc (item, per_sec, event)
Panel_item item;
unsigned int per_sec;
struct inputevent *event;
{
	frames_per_sec = per_sec;
}

Notify_value halt_proc ()
{
	/* Pause the CD-player, if available */
	if (cdplayer_pid > 0)
		(void) kill (cdplayer_pid, SIGUSR2);

	Halt = TRUE;
	return (NOTIFY_DONE);
}

static int next_frame (number)
int *number;
{
	*number += step;
	if (auto_reverse) {
		if (*number > num_images) {
			*number = num_images - 1;
			step = -1;
			(void) panel_set_value (direction_item, BACKWARD);
		} else if (*number < 1) {
			*number = 2;
			step = 1;
			(void) panel_set_value (direction_item, FORWARD);
		}
	} else if (repeat) {
		if (*number > num_images)
			*number = 1;
		else if (*number < 1)
			*number = num_images;
	} else {
		if (*number > num_images) {
			(void) message ("At end of tape");
			*number = num_images;
			return FALSE;
		} else if (*number < 1) {
			(void) message ("At beginning of tape");
			*number = 1;
			return FALSE;
		}
	}
	return TRUE;
}

#define UPDATE_FREQUENCY 10
static void sleep_awhile (tp0)
struct timeb *tp0;
{
	struct timeb tp;
	int delta_time, delay;
	static int count = 0;
	static unsigned fps_sum = 0;
	unsigned fps, value;

	/* Cannot time any faster because resolution is about 1/50-th second.
	 * This means that if frames_per_sec == MAX_FRAMES_PER_SEC,
	 * the movie will run as fast as possible */
	if (frames_per_sec >= MAX_FRAMES_PER_SEC)
		return;
	(void) ftime (&tp);
	delta_time = (tp.time - tp0->time)*1000 +
		(unsigned int) (tp.millitm - tp0->millitm);
	delay = 1000 / frames_per_sec - delta_time;
	if (delay > 5) {
		/* sleep that many microseconds */
		(void) usleep ((unsigned)delay*1000);
		fps = frames_per_sec;
	} else
		fps = 1000 / delta_time;
	fps_sum += fps;
	if (++count % UPDATE_FREQUENCY == 0) {	/* Update at every Nth frame */
		value = (unsigned) panel_get (per_sec_item, PANEL_VALUE);
		fps = fps_sum / UPDATE_FREQUENCY;
		if (fps != value)
			(void) panel_set_value (per_sec_item, fps);
		fps_sum = 0;
	}
}

static short hglass_array[] = {
#include <images/hglass.cursor>
};
mpr_static (hglass, 16, 16, 1, hglass_array);
static Cursor hourglass = (Cursor)NULL, pointer = (Cursor)NULL;

static void play_proc ()
{
	int next = TRUE;
	struct timeb tp;
	extern Notify_error notify_dispatch ();
	char handle[10];
	Rect *rect;	/* For locking */

	(void) sprintf (handle, "%d", getpid ());
	(void) notify_set_signal_func (handle, halt_proc, SIGURG, NOTIFY_ASYNC);

	/* Start the CD-player, if available */
	if (cdplayer_pid > 0)
		(void) kill (cdplayer_pid, SIGUSR1);

	/* Enable a repeat, even if at beginning/end of tape,
	 * when "Play" is pressed again */
	if (! auto_reverse && ! repeat) {
		if (current_number + step > num_images)
			current_number = 1;
		else if (current_number + step < 1)
			current_number = num_images;
	}
	Halt = FALSE;
	(void) message ("Press Stop/L1 to interrupt");
	(void) panel_set (filename_item, PANEL_SHOW_ITEM, FALSE, 0);
	current_frame = locate_frame (current_number);
	(void) window_set (canvas,	/* Disable retained canvas */
		CANVAS_RETAINED, FALSE,
		CANVAS_AUTO_CLEAR, FALSE,
		0);
	rect = (Rect *) window_get(canvas, WIN_RECT);
	if (! hourglass)
		hourglass = cursor_create (CURSOR_IMAGE, &hglass, 0);
	if (! pointer)
		pointer = cursor_copy (window_get (panel, WIN_CURSOR));
	(void) window_set (panel, WIN_CURSOR, hourglass, 0);
	/* (void) clear_proc(); /* Only if desired */
	(void) ftime (&tp);	/* Timing */
	if (! xor_rop)	/* Double buffering does not work with XOR */
		(void) pw_dbl_access (canvaspw);	/* Double buffering */

	/* Main display loop */

	playing = TRUE;
	while (! Halt) {
		(void) pw_lock (canvaspw, rect);
		(void) display_frame (current_frame);
		if (! xor_rop)
			(void) pw_dbl_flip (canvaspw);
		(void) pw_unlock (canvaspw); if (Halt) break;
		if ((next=next_frame (&current_number)) == FALSE || Halt) break;
		current_frame = locate_frame (current_number); if (Halt) break;
		(void) sleep_awhile (&tp); if (Halt) break;
		(void) ftime (&tp);
	}

	if (! xor_rop)		/* Release double buffer */
		(void) pw_dbl_release (canvaspw);
	current_frame = locate_frame (current_number);
	/* Enable retained canvas */
	(void) window_set (canvas, CANVAS_RETAINED, TRUE, 0);
	(void) display_frame (current_frame);
	(void) window_set (panel, WIN_CURSOR, pointer, 0);
	(void) panel_set_value (per_sec_item, frames_per_sec);
	(void) panel_set (filename_item, PANEL_SHOW_ITEM, TRUE, 0);
	(void) filename_proc ();
	if (next == TRUE)
		(void) message ("");
	playing = FALSE;
}

void eject_proc ()
{
	if (num_images <= 1 ||
		alert_prompt ( (Frame)window_frame, (Event *)NULL,
		ALERT_MESSAGE_STRINGS, "Are you finished with Movietool ?", 0,
		ALERT_BUTTON_YES,	"Yes",
		ALERT_BUTTON_NO,	"Continue", 0)
			== ALERT_YES)
	{

		/* Eject the CD-player, if available */
		if (cdplayer_pid > 0)
			(void) kill (cdplayer_pid, SIGTERM);

		(void) window_set (panel_frame,
			FRAME_NO_CONFIRM, TRUE, 0);
		(void) window_destroy (panel_frame);
		(void) window_set (window_frame,
			FRAME_NO_CONFIRM, TRUE, 0);
		(void) window_destroy (window_frame);
	}
}

void canvas_repaint_proc (can, pw, r)
Canvas can;
Pixwin *pw;
Rectlist *r;
{
	/* Update display of control panel */
	if (window_get (panel_frame, FRAME_CLOSED) ||
		! window_get (panel_frame, WIN_SHOW))
		window_set (panel_frame, WIN_SHOW, TRUE, 0);

	/* Update width and height */
	canvas_par.width  = (int) window_get (canvas, WIN_WIDTH);
	canvas_par.height = (int) window_get (canvas, WIN_HEIGHT);
	canvas_par.canvas_width  = (int) window_get (canvas, CANVAS_WIDTH);
	canvas_par.canvas_height = (int) window_get (canvas, CANVAS_HEIGHT);

	(void) clear_proc();
	if (current_frame > 0)
		(void) display_frame (current_frame);
	else
		step_proc ();		/* First time around */
}

static short on_array[] = {
#include <images/on.pr>
};
mpr_static (on, 64, 16, 1, on_array);

static short off_array[] = {
#include <images/off.pr>
};
mpr_static (off, 64, 16, 1, off_array);

void create_panel_items ()
{
	int slider_width;

	if (num_images <= 1)	/* 1 frame - no fancy stuff needed */
		return;

	eject_item = panel_create_item (panel, PANEL_BUTTON,
		PANEL_ITEM_X, ATTR_COLS(40),
		PANEL_ITEM_Y, ATTR_ROWS(1)-9,
		PANEL_LABEL_IMAGE, panel_button_image(panel, "Eject", 5, 0),
		PANEL_NOTIFY_PROC, eject_proc,
		0);

	play_item = panel_create_item (panel, PANEL_BUTTON,
		PANEL_ITEM_X, ATTR_COLS(1),
		PANEL_ITEM_Y, ATTR_ROWS(1)-9,
		PANEL_LABEL_IMAGE, panel_button_image(panel, "Play", 5, 0),
		PANEL_NOTIFY_PROC, play_proc,
		0);

	step_item = panel_create_item (panel, PANEL_BUTTON,
		PANEL_ITEM_X, ATTR_COLS(9),
		PANEL_ITEM_Y, ATTR_ROWS(1)-9,
		PANEL_LABEL_IMAGE, panel_button_image(panel, "Step", 5, 0),
		PANEL_NOTIFY_PROC, step_proc,
		0);

	clear_item = panel_create_item (panel, PANEL_BUTTON,
		PANEL_ITEM_X, ATTR_COLS(17),
		PANEL_ITEM_Y, ATTR_ROWS(1)-9,
		PANEL_LABEL_IMAGE, panel_button_image(panel, "Clear", 5, 0),
		PANEL_NOTIFY_PROC, clear_proc,
		0);

	repaint_item = panel_create_item (panel, PANEL_BUTTON,
		PANEL_ITEM_X, ATTR_COLS(25),
		PANEL_ITEM_Y, ATTR_ROWS(1)-9,
		PANEL_LABEL_IMAGE, panel_button_image(panel, "Repaint", 7, 0),
		PANEL_NOTIFY_PROC, canvas_repaint_proc,
		0);

	if (xor_flag) {
		xor_rop = TRUE;
		xor_show_item = panel_create_item (panel, PANEL_CYCLE,
			PANEL_ITEM_X, ATTR_COLS(40),
			PANEL_ITEM_Y, ATTR_ROWS(1)-14,
			PANEL_LABEL_STRING, "Image:",
			PANEL_MARK_IMAGES, 0,
			PANEL_NOMARK_IMAGES, 0,
			PANEL_VALUE, 0,
			PANEL_CHOICE_STRINGS, "Regular", "XORed", 0,
			PANEL_NOTIFY_PROC, xor_show_proc,
			0);

		xor_rop_item = panel_create_item (panel, PANEL_CYCLE,
			PANEL_ITEM_X, ATTR_COLS(40),
			PANEL_ITEM_Y, ATTR_ROWS(1),
			PANEL_LABEL_STRING, "XOR Play",
			PANEL_CHOICE_IMAGES, &off, &on, 0,
			PANEL_MARK_IMAGES, 0,
			PANEL_NOMARK_IMAGES, 0,
			PANEL_VALUE, xor_rop,
			PANEL_NOTIFY_PROC, xor_rop_proc,
			0);
	}
	/* frame number display for the current image on show */
	if (num_images < 10)
		slider_width = 20;
	else if (num_images > 100)
		slider_width = 200;
	else
		slider_width = 2 * num_images;

	framenumber_item = panel_create_item (panel, PANEL_SLIDER,
		PANEL_ITEM_X, ATTR_COLS(1),
		PANEL_ITEM_Y, ATTR_ROWS(2),
		PANEL_NOTIFY_LEVEL, PANEL_DONE,
		PANEL_LABEL_STRING, "Frame:",
		PANEL_SLIDER_WIDTH, slider_width,
		PANEL_VALUE, 1,
		PANEL_MIN_VALUE, 1,
		PANEL_MAX_VALUE, num_images,
		PANEL_SHOW_VALUE, TRUE,
		PANEL_NOTIFY_PROC, frame_proc,
		0);

	filename_item = panel_create_item (panel, PANEL_MESSAGE,
		PANEL_ITEM_X, ATTR_COLS(24)+slider_width,
		PANEL_ITEM_Y, ATTR_ROWS(2),
		PANEL_NOTIFY_LEVEL, PANEL_DONE,
		PANEL_SHOW_VALUE, TRUE,
		PANEL_NOTIFY_PROC, filename_proc,
		0);

	direction_item = panel_create_item (panel, PANEL_CYCLE,
		PANEL_ITEM_X, ATTR_COLS(1),
		PANEL_ITEM_Y, ATTR_ROWS(3),
		PANEL_LABEL_STRING, "Direction:",
		PANEL_VALUE, 0,
		PANEL_MARK_IMAGES, 0,
		PANEL_NOMARK_IMAGES, 0,
		PANEL_CHOICE_STRINGS, "Forward", "Backward", 0,
		PANEL_NOTIFY_PROC, direction_proc,
		0);

	reverse_item = panel_create_item (panel, PANEL_CYCLE,
		PANEL_ITEM_X, ATTR_COLS(1),
		PANEL_ITEM_Y, ATTR_ROWS(4),
		PANEL_LABEL_STRING, "Auto-reverse",
		PANEL_VALUE, 0,
		PANEL_CHOICE_IMAGES, &off, &on, 0,
		PANEL_MARK_IMAGES, 0,
		PANEL_NOMARK_IMAGES, 0,
		PANEL_SHOW_MENU, FALSE,
		PANEL_NOTIFY_PROC, reverse_proc,
		0);

	repeat_item = panel_create_item (panel, PANEL_CYCLE,
		PANEL_ITEM_X, ATTR_COLS(24),
		PANEL_ITEM_Y, ATTR_ROWS(4),
		PANEL_LABEL_STRING, "Repeat",
		PANEL_VALUE, 0,
		PANEL_CHOICE_IMAGES, &off, &on, 0,
		PANEL_MARK_IMAGES, 0,
		PANEL_NOMARK_IMAGES, 0,
		PANEL_SHOW_MENU, FALSE,
		PANEL_NOTIFY_PROC, repeat_proc,
		0);

	clipping_item = panel_create_item (panel, PANEL_CYCLE,
		PANEL_ITEM_X, ATTR_COLS(41),
		PANEL_ITEM_Y, ATTR_ROWS(4),
		PANEL_LABEL_STRING, "Clipping",
		PANEL_VALUE, 1,
		PANEL_CHOICE_IMAGES, &off, &on, 0,
		PANEL_MARK_IMAGES, 0,
		PANEL_NOMARK_IMAGES, 0,
		PANEL_SHOW_MENU, FALSE,
		PANEL_NOTIFY_PROC, clipping_proc,
		0);

	per_sec_item = panel_create_item (panel, PANEL_SLIDER,
		PANEL_ITEM_X, ATTR_COLS(1),
		PANEL_ITEM_Y, ATTR_ROWS(5) + 6,
		PANEL_NOTIFY_LEVEL, PANEL_DONE,
		PANEL_LABEL_STRING, "Frames/sec:",
		PANEL_SLIDER_WIDTH, 100,
		PANEL_VALUE, MAX_FRAMES_PER_SEC,
		PANEL_MIN_VALUE, 1,
		PANEL_MAX_VALUE, MAX_FRAMES_PER_SEC,
		PANEL_SHOW_VALUE, TRUE,
		PANEL_NOTIFY_PROC, per_sec_proc,
		0);
}
