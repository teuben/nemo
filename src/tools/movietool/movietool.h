#include <suntool/sunview.h>
#define PIXRECT_NULL (Pixrect *)NULL
#define CMS_SIZE 256	/* Colormap segment size */

typedef struct image {
	int frameno;		/* Frame's sequence number */
	char *filename;		/* Name of file containing image */
	struct pixrect *pix;	/* The basic picture */
	struct pixrect *xor;	/* The XOR of previous and this picture */
	int ras_type;		/* As in rasterfile(5): RT_BYTE_ENCODED ? */
	struct image *left;	/* Points left */
	struct image *right;	/* Points right */
	colormap_t *map;	/* Frame's colormap */
} Image;

typedef struct {	/* Canvas parameters */
	/* WIN_* attributes */
	int width;		/* Width  of region displayed within window */
	int height;		/* Height of region displayed within window */
	int depth;		/* Depth  of canvas (1 or 8 bits) */
	/* CANVAS_* attributes */
	int canvas_width;	/* Width  of entire canvas */
	int canvas_height;	/* Height of entire canvas */
} Canvas_par;
