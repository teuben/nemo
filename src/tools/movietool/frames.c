#include "movietool.h"
#include <stdio.h>
#include <pixrect/pr_io.h>

/* Global variables */

extern Canvas_par canvas_par;
extern int identical_colormaps;
extern int encoded_flag;
extern int xor_flag;
extern struct pixrect *pr_mag(), *pr_load_encoded_image();

Image *first_frame;	/* The first image of the movie sequence */

/* Local variables */

static void copy_colormap ();
static struct pixrect * pr_xor ();

int collect_images (num_files, filenames, zoom, window_cmap, dx, dy, dw, dh)
int num_files;			/* Number of file names */
char **filenames;		/* File names */
int zoom;			/* Zoom factor */
colormap_t window_cmap;		/* window's default colormap */
int dx, dy, dw, dh;		/* selecting pixrect region */
{
	Image *frame, *left_frame;
	struct pixrect *pix, *region_pix;
	colormap_t *colormap, *global_cmap;
	int num_images = 0, image;
	struct rasterfile raster;
	FILE *fp;
	
	printf("Number of files to be checked on: %d\n", num_files);

	first_frame = (Image *)malloc (sizeof (Image));
	first_frame->frameno = 1;
	first_frame->left = NULL;
	first_frame->right = NULL;
	frame = first_frame;	/* Initial frame */
	canvas_par.depth = 1;	/* Assume monochrome from the outset */

	/* Collect images in the sequence in which they are provided */

	for (image=0; image < num_files; image++) {

		printf ("%s: ", filenames[image]);
		fflush (stdout);

		if ( (fp=fopen (filenames[image], "r")) == NULL) {
			printf ("Cannot open file (skipped)\t\t\n");
			fclose (fp); continue;
		}
		if (pr_load_header (fp, &raster) == PIX_ERR) {
			printf ("Bad header (skipped)\t\t\n");
			fclose (fp); continue;
		}

		colormap = (colormap_t *) NULL;
		/* Ignore colormap (by colormap=NULL) for monochrome pictures 
	 	 * and when assuming identical colormaps */
		if (raster.ras_maptype == RMT_NONE)
			;
		else if (raster.ras_maptype == RMT_EQUAL_RGB) {
			if (identical_colormaps && canvas_par.depth == 8)
				;
			else if (raster.ras_depth == 8)	/* Allocate */
				colormap = (colormap_t *)calloc
					(1, sizeof (colormap_t));
			if (pr_load_colormap (fp, &raster, colormap)==PIX_ERR) {
				printf ("Bad colormap (skipped)\t\t\n");
				fclose (fp); free(colormap); continue;
			}
			/* Set window colormap to 1st colormap read */
			if (identical_colormaps && canvas_par.depth == 8)
				colormap = global_cmap;
			else if (colormap && canvas_par.depth != 8) {
				global_cmap = colormap;
				/* window_cmap has entries 0 and 255 modified */
				copy_colormap (colormap, &window_cmap);
				canvas_par.depth = 8;
			}
		} else {
			printf ("I don't know colormap type=%d (skipped)\t\t\n",
				raster.ras_maptype);
			fclose (fp); continue;
		}

		printf ("loading...");
		fflush (stdout);
		if (encoded_flag && raster.ras_type == RT_BYTE_ENCODED) {
			frame->ras_type = RT_BYTE_ENCODED;
			pix = pr_load_encoded_image (fp, &raster);
		} else {
			frame->ras_type = RT_STANDARD;
			pix = pr_load_image (fp, &raster, colormap);
		}
		fclose (fp);

		if (pix == PIXRECT_NULL) {
			printf ("Trouble loading image (skipped)\t\t\n");
			free(colormap); continue;
		}

		/* Save file name */
		frame->filename = filenames[image];

		/* Select pixrect subregion */

		if (dx || dy || dw || dh) {
			printf ("subregion...");
			fflush (stdout);
			region_pix = mem_create (dw, dh, pix->pr_depth);
			if (pr_rop (region_pix, 0, 0, dw, dh,
				PIX_SRC | PIX_DONTCLIP, pix, dx, dy)==PIX_ERR) {
				printf ("(couldn't copy subregion)");
				fflush (stdout);
			} else {
				pr_destroy (pix);
				pix = region_pix;
			}
		}

		/* Zoom picture */

		if (zoom && frame->ras_type == RT_STANDARD) {
			printf ("zooming...");
			fflush (stdout);
			/* Zoom in the x-direction, only */
			frame->pix = pr_mag (pix, zoom, 1);
			if (frame->pix == PIXRECT_NULL) {
				printf ("(bad zoom)");
				fflush (stdout);
				frame->pix = pix;
			} else
				pr_destroy (pix);
		} else
			frame->pix = pix;

		/* XOR'ed pixrect */

		if (xor_flag && num_images > 0) {
			printf ("XOR-ing...");
			fflush (stdout);
			frame->xor = pr_xor (frame->left->pix, frame->pix);
			if (frame->xor == PIXRECT_NULL) {
				printf ("(couldn't XOR)");
				fflush (stdout);
			}
		} else
			frame->xor = PIXRECT_NULL;

		/* Determine required canvas size */
		/* Take zooming into account ! */
		if (pix->pr_size.x*(zoom?zoom:1) > canvas_par.width)
			canvas_par.width = pix->pr_size.x*(zoom?zoom:1);
		if (pix->pr_size.y*(zoom?zoom:1) > canvas_par.height)
			canvas_par.height = pix->pr_size.y*(zoom?zoom:1);

		frame->map = colormap;
		frame->frameno = ++num_images;
		printf (" (size %dx%dx%d)",
			pix->pr_size.x * (zoom ? zoom : 1),
			pix->pr_size.y * (zoom ? zoom : 1),
			pix->pr_depth);
		printf ("   \r");	/* Carriage return, no newline */
		fflush (stdout);

		/* Allocate the next frame */

		if (image < (num_files-1)) {
			frame->right = (Image *)malloc (sizeof (Image));
			left_frame = frame;
			frame = frame->right;
			frame->left = left_frame;
		}
		frame->right = NULL;
	}
	printf ("\n");
	return (num_images);
}

static void copy_colormap (c1, c2)
colormap_t *c1, *c2;
{
	register int i;

	i = c1->length;
	c2->type   = c1->type;
	c2->length = i;
	bcopy (c1->map[0], c2->map[0], i);
	bcopy (c1->map[1], c2->map[1], i);
	bcopy (c1->map[2], c2->map[2], i);
	/* Set 1st and last entry to black and white, respectively */
	/* Disabled:
	c2->map[0][0] = c2->map[1][0] = c2->map[2][0] = 0;
	i = c1->length - 1;
	c2->map[0][i] = c2->map[1][i] = c2->map[2][i] = 255;
	*/
}

static struct pixrect *pr_xor (pr1, pr2)	/* pr1 ^ pr2 (XOR op) */
struct pixrect *pr1, *pr2;
{
	struct pixrect *xor;
	int x = pr2->pr_size.x;
	int y = pr2->pr_size.y;

	if (pr1->pr_size.x != x || pr1->pr_size.y != y ||
		pr1->pr_depth != pr2->pr_depth)
		return (PIXRECT_NULL);
	if ((xor = mem_create (x, y, pr2->pr_depth)) == PIXRECT_NULL)
		return (PIXRECT_NULL);
	if (pr_rop (xor, 0, 0, x, y, PIX_SRC | PIX_DONTCLIP, pr1, 0, 0)
		== PIX_ERR) {
		pr_destroy (xor);
		return (PIXRECT_NULL);
	}
	if (pr_rop (xor, 0, 0, x, y, PIX_SRC^PIX_DST | PIX_DONTCLIP, pr2, 0, 0)
		== PIX_ERR) {
		pr_destroy (xor);
		return (PIXRECT_NULL);
	}
	return (xor);
}
