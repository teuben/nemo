/*
 * This file is part of:
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
#include <suntool/panel.h>
#include <suntool/scrollbar.h>

#define ICONIMAGE	"movietool.icon"
#define WIDTH	300
#define HEIGHT	300
#define MAXWIDTH	1000
#define MAXHEIGHT	800
#define PANEL_FONT	"/usr/lib/fonts/fixedwidthfonts/screen.b.12"
#define AUDIOFILE	"movietool.au"

/* global declarations */

extern int	collect_images();
extern void	audio();
extern void	cdrom();
extern void	create_panel_items(), canvas_repaint_proc();
extern char	*title;
extern Image	*first_frame;

Frame		window_frame, panel_frame;
Window 		panel; /* Declaring as Window to get colormap capability */
Canvas		canvas;
static Cursor	cursor;
Pixwin		*canvaspw, *panelpw, *window_frame_pw, *panel_frame_pw;
Pixfont		*panelfont;
Canvas_par	canvas_par;
void		put_color ();
void		usage();
char		colormapname[20];
int		num_images = 0;		/* Total no. of images loaded */
int		background_color = 0;	/* Canvas background color */
		/* True if all images share the same colormap */
int		identical_colormaps	= FALSE;
		/* True if using RT_BYTE_ENCODED rasters */
int		encoded_flag	= FALSE;
		/* True if using XOR when showing successive images */
int		xor_flag	= FALSE;
		/* Zoom factor.  Zero means no zoom */
int		zoom		= 0;
int		mono_panel	= FALSE; /* If control panel is monochrome*/
int		panel_exists	= FALSE; /* If control panel exists */
int		panel_show	= FALSE; /* If control panel is shown */

/* local declarations */

static int	audio_flag	= FALSE; /* Whether audio will be output */
static int	cdrom_flag	= FALSE; /* Whether CD-ROM tool will appear */

static short	icon_image[]={		/* Movietool icon */
#include ICONIMAGE
};
mpr_static(icon_pixrect, 64, 64, 1, icon_image);


main(argc, argv)
int argc; char **argv;
{
	extern char *optarg;
	extern int optind, opterr;
	Icon icon;
	static unsigned char red[CMS_SIZE], green[CMS_SIZE], blue[CMS_SIZE];
	static colormap_t window_cmap = { RMT_NONE, 0, { red, green, blue } };
	int background_flag = FALSE;
	int i, j, argcoffset=1, width, height;
	int dx=0, dy=0, dw=0, dh=0;
	int c, errflag = 0;

	if (getenv ("WINDOW_PARENT") == NULL) {
		fprintf (stderr, "%s can only be run under SunView.\n",argv[0]);
		exit (1);
	}

	icon = icon_create (ICON_IMAGE, &icon_pixrect, 0);

	/* Open the base frame.  Note that window_create consumes SunView args*/
	window_frame = window_create (0, FRAME, 
		FRAME_ARGC_PTR_ARGV, &argc, argv,
		FRAME_CMDLINE_HELP_PROC, usage,
		FRAME_LABEL, title,
		FRAME_ICON, icon,
		FRAME_INHERIT_COLORS, TRUE,
		0);

	/* Get command line arguments */

	while ((c = getopt(argc, argv, "iz:s:XeabmCc")) != -1)
		switch (c) {
		case 'i':
			identical_colormaps = TRUE;
			printf ("Colormap of 1st frame used throughout.\n");
			break;
		case 'z':
			(void) sscanf (optarg, "%d", &zoom);
			if (zoom > 10 || zoom < 2)  {
				fprintf (stderr, "Bad arg=%d to -zoom.\n",zoom);
				zoom = 0;
				errflag++;
			} else
				(void) printf ("Zoom by factor %d.\n", zoom);
			break;
		case 's':
			j =	sscanf (optarg,		"%d", &dx) +
				sscanf (argv[optind],	"%d", &dy) +
				sscanf (argv[optind+1],	"%d", &dw) +
				sscanf (argv[optind+2],	"%d", &dh);
			if (j != 4 || dx < 0 || dy < 0 || dw < 2 || dh < 2) {
				(void) fprintf (stderr, "Bad args to -sub.\n");
				errflag++;
			} else {
				(void) printf (
				"Displaying subregion at %d,%d of size %d,%d.\n"
					, dx, dy, dw, dh);
				optind += 3;	/* Warning in "man 3 getopt" */
			}
			break;
		case 'X':
			xor_flag = TRUE;
			printf ("XOR operations will be used when playing.\n");
			break;
		case 'e':
			encoded_flag = TRUE;
			printf ("RT_BYTE_ENCODED images will be used.\n");
			break;
		case 'a':
#if defined AUDIODIR && defined AUDIOPLAYER
			audio_flag = TRUE;
			printf ("Listen...\n");
#endif
			break;
		case 'C':
		case 'c':
#if defined CDPLAYER
			cdrom_flag = TRUE;
			printf ("CD-ROM player tool will appear\n");
#endif
			break;
		case 'b':
			background_flag = TRUE;
			printf ("Canvas background colour will be set\n");
			break;
		case 'm':
			mono_panel = TRUE;
			printf ("Panel will appear in monochrome\n");
			break;
		case '?':
			errflag++;
			break;
		}

	if (encoded_flag && (dx || dy || dw || dh)) {
		printf ("-e switch is incompatible with subregion.\n");
		errflag++;
	}

	if (errflag) {
		(void) usage (argv[0]);
		exit (1);
	}

	/* Initial canvas parameters */
	canvas_par.width		= WIDTH;	
	canvas_par.height		= HEIGHT;
	canvas_par.depth		= 1;
	canvas_par.canvas_width		= WIDTH;	
	canvas_par.canvas_height	= HEIGHT;

	num_images = collect_images (argc - optind, argv + optind,
		zoom, window_cmap, dx, dy, dw, dh);

	(void) printf ("Number of frames loaded in: %d.\n", num_images);
	if (num_images <= 0) {
		(void) usage (argv[0]);
		exit (1);
	}
	if (background_flag)
		background_color = get_upper_left_color (first_frame);

	/* create a panel for the tool to hold control panels */
	if ((panelfont = pf_open(PANEL_FONT)) == NULL) {
		fprintf(stderr, "Cannot open font %s.\n", PANEL_FONT);
		exit(1);
	}
	if (num_images > 1) {
		panel_exists = TRUE;

		/* Open the panel frame */
		panel_frame = window_create (window_frame, FRAME, 
			FRAME_LABEL, "Movietool control panel",
			FRAME_SHOW_LABEL, TRUE,
			FRAME_ICON, icon,
			WIN_X, 0,
			WIN_Y, -154,
			0);
		panel = window_create (panel_frame, PANEL,
			WIN_FONT, panelfont, 
			0);
		(void) create_panel_items ();
		(void) window_fit (panel);
		(void) window_fit (panel_frame);
	}

	/* Create canvas for display purposes */

	/* Adjust window size */
	width = canvas_par.width < WIDTH ? WIDTH : canvas_par.width;
	width = width > MAXWIDTH ? MAXWIDTH : width;
	height = canvas_par.height < HEIGHT ? HEIGHT : canvas_par.height;
	height = height > MAXHEIGHT ? MAXHEIGHT : height;

	canvas = window_create (window_frame, CANVAS,
		CANVAS_RETAINED,	TRUE,
		CANVAS_FIXED_IMAGE,	FALSE,
		CANVAS_AUTO_CLEAR,	FALSE,
		CANVAS_AUTO_SHRINK,	FALSE,
		CANVAS_AUTO_EXPAND,	TRUE,
		WIN_WIDTH,		width,
		WIN_HEIGHT,		height,
		CANVAS_WIDTH,		canvas_par.width,
		CANVAS_HEIGHT,		canvas_par.height,
		CANVAS_DEPTH,		canvas_par.depth,
		CANVAS_REPAINT_PROC,	canvas_repaint_proc,
		0);

	/* Initialize canvas properties */

	/* Get pixwin for displaying movie */
	canvaspw = (Pixwin *)window_get (canvas, CANVAS_PIXWIN);

	if (canvas_par.depth == 1) {		/* MONOCHROME */
		(void) window_set (canvas,
			CANVAS_FAST_MONO, TRUE,
			0);
		if (window_get (canvas, CANVAS_FAST_MONO))
			fprintf (stderr, "Using frame-buffer's fast overlay plane\n");
		else
			fprintf (stderr, "No frame-buffer overlay plane available\n");
	} else if (canvas_par.depth == 8) {		/* COLOR */
		(void) sprintf (colormapname, "movietool%d", getpid () );
		window_frame_pw = (Pixwin *)window_get (window_frame,
			WIN_PIXWIN);
		(void) pw_setcmsname (window_frame_pw, colormapname);
		if (panel_exists) {
			panelpw = (Pixwin *)window_get (panel, WIN_PIXWIN);
			if (! mono_panel) {
				(void) pw_setcmsname (panelpw, colormapname);
				panel_frame_pw = (Pixwin *)window_get
					(panel_frame, WIN_PIXWIN);
				(void) pw_setcmsname (panel_frame_pw,
					colormapname);
			}
		}
		(void) pw_setcmsname (canvaspw, colormapname);
		(void) put_color (&window_cmap);
		if (window_frame_pw->pw_pixrect->pr_depth < 8)
			fprintf (stderr, "%s - WARNING: color frames will not make sense on a monochrome display.\n", argv[0]);
	} else {
		(void) fprintf (stderr, "%s: I don't know canvas depth=%d.\n",
			argv[0], canvas_par.depth);
		exit(1);
	}

	/* Scrollbars must be set after CANVAS_FAST_MONO (SunView, p. 125) */
	if (canvas_par.height > height)
		(void) window_set (canvas,
		WIN_VERTICAL_SCROLLBAR,
			scrollbar_create (SCROLL_LINE_HEIGHT, 20, 0),
		0);
	if (canvas_par.width > width)
		(void) window_set (canvas,
		WIN_HORIZONTAL_SCROLLBAR, 
			scrollbar_create (SCROLL_LINE_HEIGHT, 20, 0),
		0);

	/* Disable mouse inside canvas */
	cursor = window_get (canvas, WIN_CURSOR);
	(void) cursor_set (cursor,
		CURSOR_SHOW_CURSOR, FALSE,
		0);
	(void) window_set (canvas,
		WIN_CURSOR, cursor,
		0);

	(void) window_fit (canvas);
	(void) window_fit (window_frame);
/*	(void) step_proc (); */

	/* if (panel_exists)
		(void) paneledit_init (window_frame);	/* Panel editing */

	if (audio_flag)		/* Fanfare */
		(void) audio (AUDIOFILE);

	if (cdrom_flag)		/* Forking off a CD-ROM player */
		(void) cdrom ();

	(void) window_main_loop (window_frame);

	exit (0);
}

#define SHOW(string) (void) fprintf (stderr, (string))

static void usage(progname)
char *progname;
{
	(void) fprintf (stderr, "Usage: %s\n", progname);
	SHOW("       -z (or -zoom) factor - Zoom picture by factor\n");
	SHOW("       -s (or -sub) x y w h - Pick out subregion of raster\n");
	SHOW("       -i (or -ident)       - Identical colormaps used\n");
	SHOW("       -X (or -XOR)         - Play with XOR operations\n");
	SHOW("       -e (or -encoded)     - Use encoded rasters (saves RAM)\n");
	SHOW("       -a (or -audio)       - Play audio file\n");
	SHOW("       -C (or -CD-ROM)      - CD-ROM tool wanted\n");
	SHOW("       -b (or -background)  - Set background color to 1st pixel\n");
	SHOW("       -m (or -mono_panel)  - Make panel monochrome\n");
	SHOW("       file [file...]       - Sun raster files to be displayed\n");
}

/* Set the colormap for the current frame */
void put_color (cmap)
colormap_t *cmap;
{
	if ( ! cmap) return;

	(void) pw_putcolormap (window_frame_pw, 0, CMS_SIZE,
		cmap->map[0], cmap->map[1], cmap->map[2]);
	(void) pw_putcolormap (canvaspw, 0, CMS_SIZE,
		cmap->map[0], cmap->map[1], cmap->map[2]);

	/* Colors in the control panel */
	if (panel_exists && (! mono_panel)) {
		(void) pw_putcolormap (panel_frame_pw, 0, CMS_SIZE,
			cmap->map[0], cmap->map[1], cmap->map[2]);
		(void) pw_putcolormap (panelpw, 0, CMS_SIZE,
			cmap->map[0], cmap->map[1], cmap->map[2]);
	}
}

#define FLAG 0x80

/* Select the background color as the color of uppermost left pixel
 * of the first frame.  Might be made more sophisticated in the future */

static int get_upper_left_color (frame)
Image *frame;
{
	register unsigned char *image;
	int color;

	if (frame->ras_type == RT_STANDARD) {
		image = (unsigned char *) mpr_d(frame->pix)->md_image;
		color = *image;			/* First pixel of image */
	} else if (frame->ras_type == RT_BYTE_ENCODED) {
		image = (unsigned char *) frame->pix->pr_data;
		if (*image == FLAG) {		/* FLAG, Repeat count, color */
			if (*(image+1))
				color = *(image+2);
			else
				color = FLAG;	/* Count 0 means pixel==FLAG */
		} else
			color = *image;		/* First pixel of image */
	} else
		color = 0;

	return color;
}
