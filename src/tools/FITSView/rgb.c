/* Display three FITS files as red/green/blue 'planes' in an 8-bit image */

#include "ql.h"

/*#define WEDGE_HEIGHT 24*/
#define MAXPIC 4
#define WINWIDTH 1025
#define WINHEIGHT 851

PICTURE pic[MAXPIC]; 
PICTURE pic24;
int curPic = 0;

Frame frame;
Canvas canvas;
Pixfont *pf;
int is24Bit;
int zoom = 1, xsize, ysize;	/* linear zoom, size of full pixrect */
float maxv = 0, minv = 0;	/* image pixel range.  If set to zero,
				 * DATAMAX and DATAMIN from the FITS
				 * header will be used */
/* These map the available 3 or 2 bits into 8-bit values and try to
 * compensate for non linearities in the display.  They AREN'T correct yet
 */
char colorLevel3[8] = { 0, 76, 102, 127, 152, 178, 203, 255 };
char colorLevel2[8] = { 0, 110, 170, 255 };

void MakeRgb();
void repaint();

main(argc, argv)
int argc;
char **argv;
{
	register Pixwin *pw;
	Pixrect *pr;
	int i;
	int winwidth = WINWIDTH, winheight = WINHEIGHT;
	unsigned char grey[CMS_SIZE],
		rgb_red[CMS_SIZE],
		rgb_green[CMS_SIZE],
		rgb_blue[CMS_SIZE],
		rgb_black[CMS_SIZE];

	char cmd[32];
	char *cp;

        fprintf(stdout,"This is rgb V%2.1f\n\n",VERSION);
	if(argc != 11)
	   	error("rgb usage: rgb zoom red_image l h green_image l h blue_image l h\n");

	zoom = atoi(argv[1]);
	if(zoom < 1 || zoom > 16)
		error("bad zoom = %d", zoom);

	pf = pf_default();
	is24Bit=Is24Bit();

	/* put first image in red plane  */
	minv = atof(argv[3]);
	maxv = atof(argv[4]);
        sscanf(argv[2],"%s", pic[0].fn);
	LoadImage(&pic[0], 0, RED);
	pr = pic[curPic].pr;
	xsize = pr->pr_size.x;
	ysize = pr->pr_size.y;

	/* make frame, etc. */
	if(imageHdr[0]->title[0])
		cp = imageHdr[0]->title;
	else
		cp = argv[1];
	frame = window_create(0,FRAME,
		FRAME_LABEL, cp,
/*		FRAME_ARGS, argc, argv, */
		0);
	
	
	canvas = window_create(frame, CANVAS,
		CANVAS_AUTO_SHRINK, FALSE,
		CANVAS_WIDTH, xsize,
		CANVAS_HEIGHT, ysize,
		CANVAS_RETAINED, FALSE,
		CANVAS_REPAINT_PROC, repaint,
		0);
	if(xsize > winwidth) {
		window_set(canvas,
			WIN_WIDTH, winwidth,
			WIN_HORIZONTAL_SCROLLBAR, scrollbar_create(0),
			0);
	} else {
		window_set(canvas,
			WIN_WIDTH, xsize + ((ysize > winheight)? 14: 0),
			0);
	}
	if(ysize > winheight) {
		window_set(canvas,
			WIN_HEIGHT, winheight,
			WIN_VERTICAL_SCROLLBAR, scrollbar_create(0),
			0);
	} else {
		window_set(canvas,
			WIN_HEIGHT, ysize + ((xsize > winwidth)? 14: 0),
			0);
	}
	window_fit(frame);
	pw = canvas_pixwin(canvas);

	/* Set up the grey and red-green-blue arrays.  The rgb map assumes
	 * that red is in bits 7-5, green in 4-2, and blue in 1-0.  Blue
	 * is given fewer because it is not perceived well.  Black is used
	 * to enable display of a single color.
	 */
	for(i = 0; i < CMS_SIZE; i++) {
		grey[i] = i;
		rgb_red[i] = colorLevel3[i >> 5];
		rgb_green[i] = colorLevel3[(i >> 2) & 7];
		rgb_blue[i] = colorLevel2[(i & 3)];
		rgb_black[i] = 0;
	}
	pw_setcmsname(pw, "3color");
	pw_putcolormap(pw, 0, CMS_SIZE, grey, grey, grey);

	pw_rop(pw, 0, 0, pr->pr_size.x, pr->pr_size.y, PIX_SRC, pr, 0, 0);
	StartDisplay();

	/* put green image in second picture */
	minv = atof(argv[6]);
	maxv = atof(argv[7]);
        sscanf(argv[5],"%s", pic[1].fn);
	LoadImage(&pic[1], 0, GREEN);
	pr = pic[1].pr;
	pw_rop(pw, 0, 0, pr->pr_size.x, pr->pr_size.y, PIX_SRC, pr, 0, 0);

	/* put blue image in third picture */
	minv = atof(argv[9]);
	maxv = atof(argv[10]);
        sscanf(argv[8],"%s", pic[2].fn);
	LoadImage(&pic[2], 0, BLUE);
	pr = pic[2].pr;
	pw_rop(pw, 0, 0, pr->pr_size.x, pr->pr_size.y, PIX_SRC, pr, 0, 0);
	curPic = 2;

	/* Make and display the 3-color image */
	MakeRgb();
	pr = pic[curPic = 3].pr;
	pw_rop(pw, 0, 0, pr->pr_size.x, pr->pr_size.y, PIX_SRC, pr, 0, 0);
	pw_putcolormap(pw, 0, CMS_SIZE, rgb_red, rgb_green, rgb_blue);

	for(;;) {
	    printf("> ");
	    scanf("%31s", cmd);
	    switch(*cmd) {
	    case 'd':
		scanf("%d", &curPic);
		if(curPic < 0 || curPic > 3) {
			printf("pic 0 is red, 1 green, 2 blue, 3 rgb\n");
			break;
		}
		pr = pic[curPic].pr;
		pw_rop(pw, 0, 0, pr->pr_size.x, pr->pr_size.y, PIX_SRC,
			pr, 0, 0);
		if(curPic == 3) {
			pw_putcolormap(pw, 0, CMS_SIZE, rgb_red, rgb_green,
				rgb_blue);
		} else {
			pw_putcolormap(pw, 0, CMS_SIZE, grey, grey, grey);
		}
		break;
	    case 'r':
		pw_putcolormap(pw, 0, CMS_SIZE, rgb_red, rgb_black, rgb_black);
		break;
	    case 'g':
		pw_putcolormap(pw, 0, CMS_SIZE, rgb_black, rgb_green,rgb_black);
		break;
	    case 'b':
		pw_putcolormap(pw, 0, CMS_SIZE, rgb_black, rgb_black, rgb_blue);
		break;
	    case 'q':
		exit(0);
	    case 'w':
		{
		    FILE *out;
		    colormap_t colormap;
		    char fileName[64];

		    colormap.type = RMT_EQUAL_RGB;
		    colormap.map[0] = rgb_red;
		    colormap.map[1] = rgb_green;
		    colormap.map[2] = rgb_blue;
	printf("In write routine\n");
		    scanf("%s", fileName);
		    if((out = fopen(fileName, "w+")) < 0) {
			printf("Can't open ");
			perror(fileName);
			break;
		    }
		    pr_dump(pr, out, &colormap, RT_STANDARD, 0);
		    break;
		}
	    default:
		printf("cmds are:\nd n\tLoad picture n\
		\nr\tshow red 'plane' only\ng\tgreen\nb\tblue\nq\tquit\n");
	   }
	}
}


/* This takes a red image in pic0, green in pic1, and blue in pic2 and makes
 * a 3-color image in pic3 using dithering to avoid contouring.
 */
#define DITHERSIZE 4
void MakeRgb()
{
	unsigned char *rel, *gel, *bel, *rgbel;	/* pointers to pixels */
	int x, y;				/* current position */
	int linebytes;
	static char dither[DITHERSIZE][DITHERSIZE] = {
		{ 0,  8,  2, 10},
		{12,  4, 14,  6},
		{ 3, 11,  1,  9},
		{15,  7, 13,  5}
	};
	char *ditherrow, *ditherel;	/* row and element of dither matrix */
	int d;				/* current dither matrix value */
	int r, g, b, t;			/* temporaries */

	rel = (unsigned char *)mpr_d(pic[0].pr)->md_image;
	gel = (unsigned char *)mpr_d(pic[1].pr)->md_image;
	bel = (unsigned char *)mpr_d(pic[2].pr)->md_image;
	if(!(pic[3].pr = mem_create(xsize, ysize, 8)) )
		error("Can't create rgb Pixrect");
	CheckSame(0, 3);
	rgbel = (unsigned char *)mpr_d(pic[3].pr)->md_image;
	linebytes = mpr_mdlinebytes(pic[0].pr);
	ditherrow = dither[0];
	for(y = 0; y < ysize; y++) {
	    ditherel = ditherrow;
	    for(x = 0; x < xsize; x++) {
		d = *ditherel;
		/* Given three bits, the maximum value is 0xe0, not 255 */
		t = (rel[x] * 7) >> 3;
		/* Red = the high three bits.  If the next three bits
		 * are > the dither value, increase the high part */
		if((r = t & 0xe0) != 0xe0 &&
		    ((t >> 1) & 0x0f) > d)
			r += 0x20;
		
		/* same for green, but shift the three bits down to match the
		 * colormap. */
		t = (gel[x] * 7) >> 3;
		if((g = t >> 3 & 0x1c) != 0x1c &&
		    ((t >> 1) & 0x0f) > d)
			g += 4;

		/* Blue is similar, but we use only two bits */
		t = (bel[x] * 3) >> 2;
		if((b = t >> 6 & 3) != 3 &&
		    ((t >> 2) & 0x0f) > d)
			b++;

		rgbel[x] = r | g | b;

		if(ditherel > ditherrow + DITHERSIZE - 2)
			ditherel = ditherrow;
		else
			ditherel++;
	    }
	    rel += linebytes;
	    gel += linebytes;
	    bel += linebytes;
	    rgbel += linebytes;
	    if(ditherrow > dither[DITHERSIZE - 2])
		ditherrow = dither[0];
	    else
	    	ditherrow += DITHERSIZE;
	}
}

/* Check that pic1 through pic2 have the same size */
CheckSame(pic1, pic2)
int pic1, pic2;
{
	Pixrect *pr;

	for( ; pic1 <= pic2; pic1++) {
		pr = pic[pic1].pr;
		if(pr->pr_size.x != xsize)
			error("Image %d: wid = %d, xsize = %d\n",pic1,
				pr->pr_size.x, xsize);
		if(pr->pr_size.y != ysize)
			error("Image %d: height = %d, ysize = %d\n",pic1,
				pr->pr_size.y, ysize);
	}
}

static void repaint(canvas, pw, repaint_area)
Canvas canvas;
Pixwin *pw;
Rectlist *repaint_area;
{
	Pixrect *pr = pic[curPic].pr;

	pw_rop(pw, 0, 0, pr->pr_size.x, pr->pr_size.y, PIX_SRC, pr, 0, 0); 
}
