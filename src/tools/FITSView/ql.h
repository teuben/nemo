#include <suntool/sunview.h>
#include <suntool/canvas.h>
#include <math.h>
#include <stdio.h>
#include "fits.h"

#define VERSION 3.0
#define CMS_SIZE 256
#if !defined(WEDGE_HEIGHT)
#	define WEDGE_HEIGHT 32
#endif

typedef struct {
	Pixrect *pr;
	char fn[64];
	char label[60];
} PICTURE;
extern PICTURE pic[];
extern PICTURE pic24;
extern Pixfont *pf;
extern int is24Bit;
extern int colorMapSize;
extern int curPic;
extern int zoom, xsize, ysize;	/* linear zoom, size of full pixrect */
extern float maxv, minv;	/* image pixel range.  If set to zero,
				 * DATAMAX and DATAMIN from the FITS
				 * header will be used */
extern int display_type;	/* Cycle choice between grey scale (0) and
				 * 3-color(1) */


/* loadimage.c*/

extern void LoadImage(/*p, use24Bits, offset*/);
extern void make_wedge(/*pr, x, y, w, h, minLabel, maxLabel*/);
extern void StartDisplay();
extern int Is24Bit();
extern Pixrect *CreatePr24();
extern ShowCurrent(/*pw*/);
extern ConvertToPr24(/*srcPic, how*/);
int invert_pixels;

#define RED 	3
#define BLUE	2
#define GREEN	1
#define GREY	0
