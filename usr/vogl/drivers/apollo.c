/*
 *	Apollo driver for vogl.
 */

#include "/sys/ins/base.ins.c"
#include "/sys/ins/gpr.ins.c"
#include "/sys/ins/pad.ins.c"
#include "/sys/ins/kbd.ins.c"
#include <stdio.h>
#include "vogl.h"

#define FONTBASE "/sys/dm/fonts/"

#define MAXCOLORS		256
#define MIN(x,y)		((x) < (y) ? (x) : (y))
#define COLOR_ENTRY(r, g, b)	(gpr_$pixel_value_t) ((r << 16) | (g << 8) | b)

static ios_$id_t		stream_id;
static gpr_$pixel_value_t  	color_value[MAXCOLORS];
static gpr_$pixel_value_t  	old_color_value[MAXCOLORS];

static gpr_$offset_t        	init_bitmap_size;   
static gpr_$bitmap_desc_t   	current_bitmap, front_bitmap, back_bitmap;       
static gpr_$attribute_desc_t   	attribs;
static gpr_$window_t	   	source;
static gpr_$position_t	   	dest_pos;

static gpr_$display_mode_t  	mode = gpr_$direct;
static gpr_$plane_t         	hi_plane_id = 0;    
static boolean              	delete_display;    
static status_$t            	status;
static gpr_$keyset_t        	keys, mouse_keys;

static pad_$window_list_t 	window_info;
static short 			n_windows;
static gpr_$disp_char_t		disp;
static short			disp_len;

static int 			first_time, back_used;
static int 			current_color;
static short 			font_id;

/*
 * Do nothing
 *
 */
int
noop()
{
	return(-1);
}

/*
 * APOLLO_init
 *
 *	initialises window to occupy current window
 */
APOLLO_init()
{
	pad_$window_desc_t 	window;
	int 	size, prefx, prefy, prefxs, prefys, x, y, w, h;
	short 	i;

	pad_$set_scale(stream_$stdout,
		       1,
		       1,
		       status);

	pad_$inq_windows(stream_$stdout,
			 window_info,
				  10,
			   n_windows,
			     status);

	w = window_info[0].width;
	h = window_info[0].height;
	x = window_info[0].left;
	y = window_info[0].top;

	getprefposandsize(&prefx, &prefy, &prefxs, &prefys);

	if (prefx > -1) {
		x = prefx;
		y = prefy;
	}

	if (prefxs > -1) {
		w = prefxs;
		h = prefys;
	}

	size = MIN(w, h);
	vdevice.sizeX = vdevice.sizeY = size - 1;
	vdevice.sizeSx = w - 1;
	vdevice.sizeSy = h - 1;

        init_bitmap_size.x_size = w;
        init_bitmap_size.y_size = h;

	source.window_base.x_coord = source.window_base.y_coord = 0;
	source.window_size.x_size = init_bitmap_size.x_size;
	source.window_size.y_size = init_bitmap_size.y_size;
	dest_pos.x_coord = dest_pos.y_coord = 0;

	/*
	 * Inquire about the actual display ....
	 */

	gpr_$inq_disp_characteristics(mode,
				      stream_$stdout,
				      (short)60,
				      disp,
				      disp_len,
				      status);

	vdevice.depth = disp.n_planes;
	hi_plane_id = disp.n_planes - 1;

	if (prefx == -1 && prefxs == -1) {
		stream_id = stream_$stdout;
	} else {
		window.top = y;
		window.left = x;
		window.width = w;
		window.height = h;

		pad_$create_window("", 0,
				pad_$transcript, 
				1, 
				window,
				stream_id,
				status
		);

		pad_$set_auto_close(stream_id, 1, true, status);
	}

        gpr_$init(mode,
                  stream_id,
                  init_bitmap_size,
                  hi_plane_id,
                  front_bitmap,
                  status);

	current_bitmap = front_bitmap;

	gpr_$set_auto_refresh(true, status);
  
	gpr_$set_cursor_active(false,status);

	/*  Set up all the character stuff  */

	first_time = 1;

	/*  create a key set for the event interupts  */

	lib_$init_set(keys, (short)256);
	lib_$init_set(mouse_keys, (short)6);

	lib_$add_to_set(mouse_keys, (short)6, KBD_$M1D);
	lib_$add_to_set(mouse_keys, (short)6, KBD_$M2D);
	lib_$add_to_set(mouse_keys, (short)6, KBD_$M3D);
	lib_$add_to_set(mouse_keys, (short)6, KBD_$M1U);
	lib_$add_to_set(mouse_keys, (short)6, KBD_$M2U);
	lib_$add_to_set(mouse_keys, (short)6, KBD_$M3U);

	for (i = 0; i < 128; i++) 
		lib_$add_to_set(keys, (short)256, (short)i);

	gpr_$enable_input(gpr_$keystroke, keys, status);
	gpr_$enable_input(gpr_$buttons, keys, status);
	gpr_$enable_input(gpr_$locator, keys, status);

	/*  set default color (colour)   */

	if (disp.n_planes > 1) {
		gpr_$inq_color_map(0L, (short)MAXCOLORS, old_color_value, status);
		gpr_$inq_color_map(0L, (short)MAXCOLORS, color_value, status);

		color_value[0] = COLOR_ENTRY(0,0,0);         /* color--black   */
		color_value[1] = COLOR_ENTRY(255,0,0);       /* color--red     */
		color_value[2] = COLOR_ENTRY(0,255,0);       /* color--green   */
		color_value[3] = COLOR_ENTRY(255,255,0);     /* color--yellow  */
		color_value[4] = COLOR_ENTRY(0,0,255);       /* color--blue    */
		color_value[5] = COLOR_ENTRY(255,0,255);     /* color--magenta */
		color_value[6] = COLOR_ENTRY(0,255,255);     /* color--cyan    */
		color_value[7] = COLOR_ENTRY(255,255,255);   /* color--white   */

		/* modify color table */
		gpr_$acquire_display(status);
		gpr_$set_color_map((long)0,
				    (short)MAXCOLORS,
				    color_value,
				    status);
		gpr_$release_display(status);
	}

	back_used = 0;

	return(1);
}


/*
 * APOLLO_draw
 *
 *	Draw a lines from the current graphics position to (x, y).
 */
APOLLO_draw(x, y)
	int	x, y;
{
	gpr_$acquire_display(status);
	gpr_$move((short)vdevice.cpVx, (short)(vdevice.sizeSy - vdevice.cpVy), status);
	gpr_$line((short)x, (short)(vdevice.sizeSy - y), status);
	gpr_$release_display(status);
	vdevice.cpVx = x;
	vdevice.cpVy = y;
}

/*
 * APOLLO_getkey
 *
 *	grab a character from the keyboard
 */
int
APOLLO_getkey()
{
	gpr_$event_t     et;
	char		 ed;
	gpr_$position_t  pos;


	et = gpr_$no_event;
	while (et != gpr_$keystroke)
		(void)gpr_$event_wait(et, ed, pos, status);

	/*
	 * What for this stupid gpr event system map the return key
	 * to a SYN character
	 */

	if ((ed & 0x7f) == 0x16)
		ed = 13;

	return(ed);
}

/*
 * APOLLO_checkkey
 *
 *	checks if there is key waiting in the keyboard
 */
int
APOLLO_checkkey()
{
	gpr_$event_t     et;
	char             ed;
	gpr_$position_t  pos;

	et = gpr_$no_event;

	(void)gpr_$cond_event_wait(et, ed, pos, status);

	return((et == gpr_$keystroke ? 1 : 0));
}

/*
 * APOLLO_locator
 *
 *      return the window location of the cursor, plus which mouse button,
 * if any, is been pressed.
 *
 *	LOCATOR - needs to return straight away.
 */
int
APOLLO_locator(wx, wy)
	int	*wx, *wy;
{
	gpr_$event_t     et;
	char             ed;
	gpr_$position_t  pos;
	gpr_$position_t  origin;

	gpr_$bitmap_desc_t	curs_pat;
	gpr_$raster_op_array_t	curs_raster_op;
	boolean			active;
	


	(void)gpr_$cond_event_wait(et, ed, pos, status);

	gpr_$inq_cursor(curs_pat, curs_raster_op, active, pos, origin, status);

	*wx = (int)pos.x_coord;
	*wy = (int)vdevice.sizeSy - (int)pos.y_coord;

	if (ed < 'a')			/* absorb button up */
		return(0);

	return(1 << ((int)ed - 'a'));
}

/*
 * APOLLO_clear
 *
 *	clear the window to the current color.
 */
APOLLO_clear()
{
        int     x[4], y[4];

        if (vdevice.maxVx != vdevice.sizeSx
           || vdevice.maxVy != vdevice.sizeSy
           || vdevice.minVx != 0
           || vdevice.minVy != 0) {
		x[0] = x[3] = vdevice.minVx;
		y[0] = y[1] = vdevice.maxVy;
		y[2] = y[3] = vdevice.minVy;
		x[1] = x[2] = vdevice.maxVx;

                APOLLO_fill(4, x, y);
        } else {
		gpr_$acquire_display(status);
		gpr_$clear((long)current_color, status);
		gpr_$release_display(status);
        }

}

/*
 * APOLLO_exit
 *
 *	reset the window back to normal mode (sigh)
 */
APOLLO_exit()
{
	if (disp.n_planes > 1) {
		gpr_$acquire_display(status);
		gpr_$set_color_map((long)0,
				    (short)MAXCOLORS,
				    old_color_value,
				    status);
		gpr_$release_display(status);
	}
	gpr_$terminate(delete_display,status); /*Terminate gpr.*/
}

/*
 * APOLLO_color
 *
 *	change the drawing color.
 */
APOLLO_color(ind)
	int	ind;
{
	if (disp.n_planes <= 1) {
		if (ind > 0)
			current_color = 1;

		if (disp.invert)
			current_color = !current_color;

	} else {
		current_color = ind % MAXCOLORS;
	}

	gpr_$set_draw_value((int)current_color, status); 
	gpr_$set_text_value((int)current_color, status);
	gpr_$set_fill_value((int)current_color, status);
	/*
	 * GPR manual says that this sets text background to
	 * 'transparent'
	 */
	gpr_$set_text_background_value((int)-1, status);
}

/*
 * APOLLO_mapcolor
 *
 *	set a colormap entry.
 */
APOLLO_mapcolor(in, r, g, b)
	int	in;
	short	r, g, b;
{
	if (in < 0 || in > MAXCOLORS)
		verror("APOLLO_mapcolor: index out of range");
	
	color_value[in] = COLOR_ENTRY(r, g, b);

	/* modify color table */

	gpr_$acquire_display(status);
	gpr_$set_color_map((long)0,
			    (short)MAXCOLORS,
			    color_value,
			    status);
	gpr_$release_display(status);
}

/*
 * APOLLO_fill
 *
 *	draw a filled polygon
 */
APOLLO_fill(n, x, y)
	int	n;
	int	*x, *y;
{
	short 	nx[256], ny[256], startx, starty;
	int	i;

	for (i = 1; i < n; i++) {
		ny[i-1] = vdevice.sizeSy - y[i];
		nx[i-1] = x[i];
	}

	startx = x[0];
	starty = vdevice.sizeSy - y[0];

	gpr_$acquire_display(status);
	gpr_$start_pgon(startx, starty, status);

	gpr_$pgon_polyline(nx, ny, (short)(n - 1), status);
	gpr_$close_fill_pgon(status);
	gpr_$release_display(status);
}

/*
 * APOLLO_font
 *
 *	load in a hardware font
 */
APOLLO_font(fontfile)
	char	*fontfile;
{
	gpr_$offset_t start;
	short xy_end;
	char	buf[256];

	if (first_time) {
		gpr_$unload_font_file(font_id, status);
		first_time = 0;
	}

	strcpy(buf, FONTBASE);
	strcat(buf, fontfile);
	gpr_$load_font_file(buf, (short)strlen(buf), font_id, status);

	if (status.all) 
		return (0);

	gpr_$set_text_font(font_id, status);

	/* see above comment about &'s */

	gpr_$inq_text_offset("H", (short)1, start, xy_end, status);

	vdevice.hwidth = xy_end;
	vdevice.hheight = start.y_size;

	return (1);
}

/*
 * APOLLO_char
 *
 *	output a character
 */
APOLLO_char(c)
	char	c;
{
	char s[2];
	s[1] = '\0';
	s[0] = c;
	gpr_$acquire_display(status);
	gpr_$move((short)vdevice.cpVx, (short)(vdevice.sizeSy - vdevice.cpVy), status);
	gpr_$text(s[0], (short)1, status);
	gpr_$release_display(status);
}

/*
 * APOLLO_string
 *
 *	display a string
 */
APOLLO_string(s)
	char	s[];
{
	short	len;

	len = (short)strlen(s);
	gpr_$acquire_display(status);
	gpr_$move((short)vdevice.cpVx, (short)(vdevice.sizeSy - vdevice.cpVy), status);
	gpr_$text(s[0], len, status);
	gpr_$release_display(status);
}

/*
 * APOLLO_backb
 *
 *	Allocates and sets the back display buffer.
 */
int
APOLLO_backb()
{
	if (!back_used) {
		gpr_$allocate_attribute_block(attribs, status);
		gpr_$allocate_bitmap(init_bitmap_size, (short)hi_plane_id, attribs, back_bitmap, status);
		gpr_$enable_input(gpr_$keystroke, keys, status);
		gpr_$enable_input(gpr_$buttons, keys, status);
		gpr_$enable_input(gpr_$locator, keys, status);
		back_used = 1;
		if (status.all) {
			fprintf(stderr, "APOLLO_backb: problem (status = %d)\n", status);
			error_$print(status);
			return(-1);
		}
	}

	gpr_$set_bitmap(back_bitmap, status);
	gpr_$set_color_map((long)0,
			    (short)MAXCOLORS,
			    color_value,
			    status);

	current_bitmap = back_bitmap;

	return(1);
}

/*
 * APOLLO_swapb
 *
 *	Swap the front and back buffers' role - actually just copy
 *	the back buffer to the display buffer.
 */
APOLLO_swapb()
{
	gpr_$acquire_display(status);
	gpr_$set_bitmap(front_bitmap, status);
	gpr_$pixel_blt(back_bitmap, source, dest_pos, status);
	gpr_$set_bitmap(back_bitmap, status);
	gpr_$release_display(status);
}

/*
 * APOLLO_frontb
 *
 *	Make sure we are drawing in the front (display) buffer
 */
APOLLO_frontb()
{
	current_bitmap = front_bitmap;
	gpr_$set_bitmap(front_bitmap, status);
}


static DevEntry apdev = {
	"apollo",
	"f16.b",
	"f7x13.b",
	APOLLO_backb,
	APOLLO_char,
	APOLLO_checkkey,
	APOLLO_clear,
	APOLLO_color,
	APOLLO_draw,
	APOLLO_exit,
	APOLLO_fill,
	APOLLO_font,
	APOLLO_frontb,
	APOLLO_getkey,
	APOLLO_init,
	APOLLO_locator,
	APOLLO_mapcolor,
	noop,
	noop,
	APOLLO_string,
	APOLLO_swapb,
	noop
};

/*
 * _APOLLO_devcpy
 *
 *      copy the apollo device into vdevice.dev.
 */
_APOLLO_devcpy()
{
    vdevice.dev = apdev;
}
