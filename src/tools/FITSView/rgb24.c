/* Display three FITS files as red/green/blue in a 24-bit image */

#include <math.h>
#include "ql.h"

#define MAXPIC 1
#define WINWIDTH 1025
#define WINHEIGHT 851

PICTURE pic[MAXPIC];
int curPic = 0;

Frame frame;
Canvas canvas;
Pixfont *pf;
int is24Bit;
PICTURE pic24;
int zoom = 1, xsize, ysize;	/* linear zoom, size of full pixrect */
float maxv = 0, minv = 0;	/* image pixel range.  If set to zero,
				 * DATAMAX and DATAMIN from the FITS
				 * header will be used */
void repaint();

main(argc, argv)
int argc;
char **argv;
{
	register Pixwin *pw;
	Pixrect *pr;
	double d;
	int winwidth = WINWIDTH, winheight = WINHEIGHT;
	char cmd[32];
	char *cp;

        fprintf(stdout,"This is rgb24 V%2.1f\n\n",VERSION);
	if(argc != 11)
		error("rgb usage: rgb zoom fn1 l h fn2 l h fn3 l h\n");

	zoom = atoi(argv[1]);
	if(zoom < 1 || zoom > 16)
		error("bad zoom = %d", zoom);

	pf = pf_default();
	if( !(is24Bit = Is24Bit()))
		error("rgb24 only runs on 24 bit frame buffers");

	/* put first image in blue plane */
	minv = atof(argv[3]);
	maxv = atof(argv[4]);
        sscanf(argv[2],"%s", pic24.fn);
	LoadImage(&pic24,1,1);
	pr = pic24.pr;
	xsize = pr->pr_size.x;
	ysize = pr->pr_size.y;

	/* make frame, etc. */
	if(imageHdr[0]->title[0])
		cp = imageHdr[0]->title;
	else
		cp = argv[1];
	frame = window_create(0,FRAME,
		FRAME_LABEL, cp,
		0);
	
	
	canvas = window_create(frame, CANVAS,
		CANVAS_AUTO_SHRINK, FALSE,
		CANVAS_WIDTH, xsize,
		CANVAS_HEIGHT, ysize,
		CANVAS_RETAINED, FALSE,
		CANVAS_REPAINT_PROC, repaint,
		CANVAS_COLOR24, TRUE,
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

	gammacorrect(pw->pw_pixrect, .37);

	ShowCurrent(pw);
	StartDisplay();

	/* put second image in green plane */
	minv = atof(argv[6]);
	maxv = atof(argv[7]);
        sscanf(argv[5],"%s", pic24.fn);
	LoadImage(&pic24,1, 2);
	ShowCurrent(pw);

	/* put third image in red picture */
	minv = atof(argv[9]);
	maxv = atof(argv[10]);
        sscanf(argv[8],"%s", pic24.fn);
	LoadImage(&pic24,1, 3);
	ShowCurrent(pw);

	for(;;) {
	    printf("> ");
	    scanf("%31s", cmd);
	    switch(*cmd) {
	    case 'g':
		scanf("%lf", &d);
		gammacorrect(pw->pw_pixrect, d);
		break;
	    case 't':
		writeTGA();
		break;
	    case 'q':
		exit(0);
		    break;
	    default:
		printf("cmds are:\ng\tgamma correct\nt\twrite TGA file (tga.out)\nq\tquit\n");
	   }
	}
}

gammacorrect(pr, power)
Pixrect *pr;
double power;
{
	unsigned char red[256];
	int i;

	for(i = 0; i < 256; i++) {
		if(power > 0)
			red[i] = 255 * pow(i / 255., power);
		else if(power == 0)
			red[i] = i;
		else
			red[i] = - power;
	}
	pr_putlut(pr, 0, 255, red, red, red);
}

static void repaint(canvas, pw, repaint_area)
Canvas canvas;
Pixwin *pw;
Rectlist *repaint_area;
{
	ShowCurrent(pw);
}
