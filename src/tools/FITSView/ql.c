/* This has a simple (no data compression) PostScript hardcopy routine - mwp */
/* also has TGA output routine */

#include "ql.h"

#define MAXPIC 2
#define WINWIDTH 1025
#define WINHEIGHT 851

PICTURE pic[MAXPIC];
int curPic=0;

Frame frame;
Canvas canvas;
Pixfont *pf;
int is24Bit;
PICTURE pic24;
int zoom = 1, xsize, ysize;	/* linear zoom, size of full pixrect */
float maxv = 0, minv = 0;	/* image pixel range.  If set to zero,
				 * DATAMAX and DATAMIN from the FITS
				 * header will be used */
double gammaCor = .37;		/* Actually 1/Gamma */
void repaint();

main(argc, argv)
int argc;
char **argv;
{
	register Pixwin *pw;
	Pixrect *pr;
	int i;
	int winwidth = WINWIDTH, winheight = WINHEIGHT;

	char cmd[32];
	char *cp;
	char c;

        fprintf(stdout,"This is ql V%2.1f\n\n",VERSION);
	if(argc != 5)
		error("ql usage: ql zoom fn l h\n");
	invert_pixels=1;
	zoom = atoi(argv[1]);
	if(zoom < 1 || zoom > 16)
		error("bad zoom = %d", zoom);

	pf = pf_default();
	is24Bit = Is24Bit();

	/* put image in first picture */
	minv = atof(argv[3]);
	maxv = atof(argv[4]);
	if(is24Bit) {
		strncpy(pic24.fn, argv[2], 64);
		LoadImage(&pic24, 1, GREY);
		pr = pic24.pr;
	} else {
		colorMapSize = 128;
		strncpy(pic[0].fn, argv[2], 64);
		LoadImage(&pic[0], 0, GREY);
		pr = pic[0].pr;
	}
	xsize = pr->pr_size.x;
	ysize = pr->pr_size.y;

	/* make frame, etc. */
	if(imageHdr[0]->title[0])
		cp = imageHdr[0]->title;
	else
		cp = argv[1];
	if(is24Bit)
		strcpy(pic24.label, cp);
	else
		strcpy(pic[0].label, cp);
	frame = window_create(0,FRAME,
		FRAME_LABEL, cp,
		0);
	
	canvas = window_create(frame, CANVAS,
		CANVAS_AUTO_SHRINK, FALSE,
		CANVAS_WIDTH, xsize,
		CANVAS_HEIGHT, ysize,
		CANVAS_RETAINED, FALSE,
		CANVAS_REPAINT_PROC, repaint,
		0);
	if(is24Bit)
		window_set(canvas,
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
	if(!is24Bit) {
		pw_setcmsname(pw, "halfgrey");
	}
	gammacorrect(pw, gammaCor);
	ShowCurrent(pw);
	StartDisplay();

	for(;;) {
	    printf("> ");
	    scanf("%31s", cmd);
	    switch(*cmd) {
	    case 'g':
		scanf("%lf", &gammaCor);
		gammacorrect(pw, gammaCor);
		break;
	    case 'h':
		hardcopy(pw);
		break;
	    case 'i':
		invert_pixels = -invert_pixels;
		ShowCurrent(pw);
		break;
	    case 'p':
		scanf("%d", &i);
		pseudocolor(pw, gammaCor, i);
		break;
	    case 's':
		scanf("%f %f", &minv, &maxv);
		if(is24Bit) {
			LoadImage(&pic24, 1, GREY);
		} else {
			LoadImage(&pic[0], 0, GREY);
		}
		ShowCurrent(pw);
		break;
	    case 't':
		writeTGA();
		break;
	    case 'q':
		exit(0);
	    default:
		printf("cmds are:\
\ng gamma(.38)\tChange display gamma correction\
\nh hardcopy\tgenerate a PostScript hardcopy\
\ni invert\tinvert pixels\
\np nsteps(128)\tChange to pseudo color display\
\ns l h\t\tRescale picture\
\nt TGA\t\tgenerate a TRUEVISION TGA file (tga.out) from image\
\nq quit\n");
	   }
	}
}

gammacorrect(pw, power)
Pixwin *pw;
double power;
{
	unsigned char map[256];
	int i;
	double scale = 1. / (colorMapSize - 1);

	if(power == 0)
		scale *= 255;
	for(i = 0; i < colorMapSize; i++) {
		if(power > 0)
			map[i] = 255 * pow(i * scale, power);
		else if(power == 0)
			map[i] = i * scale;
		else
			map[i] = - power;
	}
	if(is24Bit)
		pr_putlut(pw->pw_pixrect, 0, 255, map, map, map);
	else
		pw_putcolormap(pw, 0, colorMapSize, map, map, map);
}

pseudocolor(pw, power, nsteps)
Pixwin *pw;
double power;
int nsteps;
{
#define HMIN .01
#define HMAX .90
	double h, s, v, r, g, b;	/* hue, sat, value, red, green, blue */
	double factor;
	int i, cms, level;
	unsigned char red[256], green[256], blue[256];

	if(power <= 0)
		power = 1.0;
	cms = (is24Bit)? 256: colorMapSize;
	factor = 1.0 / (double)(nsteps - 1);
	for(i = 0; i < cms; i++) {
		level = (cms - 1 - i) * (nsteps) / cms;
		h = level * factor;
		if(h < HMIN) {
			s = h / HMIN;
			v = 1.0;
			h = 0;
		} else if(h > HMAX) {
			v = (1 - h) / (1 - HMAX);
			s = h = 1.0;
		} else {
			h = (h - HMIN) / (HMAX - HMIN);
			s = v = 1.0;
		}
		rainbow(h, s, v, &r, &g, &b);
		if(i < 10 || i > cms - 10) {
		}
		red[i] = 255 * pow(r, power);
		green[i] = 255 * pow(g, power);
		blue[i] = 255 * pow(b, power);
	}
	if(is24Bit)
		pr_putlut(pw->pw_pixrect, 0, 255, red, green, blue);
	else
		pw_putcolormap(pw, 0, colorMapSize, red, green, blue);
}

static void repaint(canvas, pw, repaint_area)
Canvas canvas;
Pixwin *pw;
Rectlist *repaint_area;
{
	ShowCurrent(pw);
}

/* convert the image to postscript  and print out */

#define INCH 72
#define PAGEX 8.5   /* Useable page is 8 x 10.5 inches (.25 inch margins) */
#define PAGEY 11    /* Useable area is 576 x 756 pixels			  */ 
#define MARGIN 0.25

hardcopy(pw)
Pixwin *pw;
{
	FILE *fout;
	Pixrect *pr;	
	register unsigned char *start,*image, pixval;
	int n,i,j,xs,ys,linebytes,fd, status;
	float scale,yscale, xscale, xoffset, yoffset, aspect;
	char hbuf[512];
	char tempfile[16];
	char dispose[80];
	char out[4];
	char *op, *s, *getenv();
	
	if(is24Bit) {
		fprintf(stderr,"hardcopy only works on 8-bit screens. sorry\n");
/*		fout = fopen("/dev/null","a");
		fprintf(stderr,"Output to /dev/null..."); */
		return(-1); 
	}

	for(i=0; i< 4;out[i++]='\0')
		;

	if((s = getenv("R_DISPOSE")) == NULL) {
		fprintf(stderr,"You must define R_DISPOSE in your environment.\nSee FITSView documentation file for details.\n");
		return(-1);
	}

        /* Create a temporary file to hold PostScript program. */

        strcpy (tempfile, "/tmp/psXXXXXX");
        if ((fd = mkstemp (tempfile)) == -1) {
            	fprintf (stderr, "cannot create temporary file %s\n", tempfile);
            	return (-1);
        } else
            	fout = fdopen (fd, "a");

	sprintf(dispose,s,tempfile);

       /* Initialize the hbuf array, which contains the hex encoded
        * representations of the 256 possible grey-level binary byte values.
        */
        for (n=0, op=hbuf;  n < 256 ;n++) {
            i = ((n >> 4) & 017);
            *op++ = (i < 10) ? i + '0' : (i-10) + 'A';
            i = (n & 017);
            *op++ = (i < 10) ? i + '0' : (i-10) + 'A';
        }

	pr = pic[0].pr;
/* lock the frame buffer during readout */
/*	pw_lock(pw,&pr)*/

	start = (unsigned char *) mpr_d(pr)->md_image;
	linebytes = mpr_d(pr)->md_linebytes;
	xs = pr->pr_size.x;
	ys = pr->pr_size.y;
	/*fprintf(stderr,"start %d linebytes: %d\n",start,linebytes);*/
	
	/* scale the page according by smallest pagesize/imagesize ratio 
 	 * and center image on page, keeping the aspect ratio.
	 */ 

	xscale = INCH*(PAGEX-MARGIN*2)/(float)xs;
	yscale = INCH*(PAGEY-MARGIN*2)/(float)ys;
	aspect = (float) xs/ys;
	if(xscale < yscale) {
	   yscale = xscale/aspect;
	}
	else {
	   xscale = yscale*aspect;
	}
	xoffset = (PAGEX - xscale)/2;
	yoffset = (PAGEY - yscale)/2;
/*	fprintf(stderr,"xoffset=%4.2f inches; yoffset=%4.2f inches\n",xoffset,yoffset);
	fprintf(stderr,"image is %d by %d pixels: aspect=%4.2f\n",xs,ys,aspect); 
	fprintf(stderr,"image is %4.2f inches by %4.2f inches\n",xscale,yscale);*/
	fprintf(fout,"%%! ql hardcopy\n");
	fprintf(fout,"/inch {72 mul} def\n");
	fprintf(fout,"%4.2f inch %4.2f inch translate\n",xoffset, yoffset);
	fprintf(fout,"%4.2f inch %4.2f inch scale\n",xscale,yscale); 
	fprintf(fout,"%d %d 8 [%d 0 0 -%d 0 %d] { <",xs,ys,xs,ys,ys);

/* output the image,  converting pixel values to hex  				*/
/* multiply pixval by 2 because colormap has 128 values but postscript uses 256 */

	for(i=0;i<ys;i++) {
	  	for(j=0;j<xs;j++) {
			image=start+j+i*linebytes;
			pixval = (*image)*2;

			if(invert_pixels==-1)
				pixval=254-pixval;	

/* add 2 to array address so that white is FF and not FE. "0" will be 02 instead
 * of 00, but its still pretty damn black on the page.
 */ 
			op = hbuf + 2 + (pixval*2);
			out[0] = *op++;
			out[1] = *op; 
			fprintf(fout,"%2s",out);
		}
	}
	fprintf(fout,"> } image \n");

/* put the image label (from FITS header) on the page */

	fprintf(fout,"stroke\n");
	fprintf(fout,"/Times-Roman findfont 15 scalefont setfont\n");
	fprintf(fout,"0 %4.2f inch moveto\n",yscale+.5);
	/*fprintf(fout,"0 setgray (%s) show\n",pic[0].label);*/
	fprintf(fout,"0 setgray (MBM12) show\n");
	fprintf(fout,"showpage\n");
	fclose(fout);
	close(fd);

/* now send it to the printer */
	if( (status=system(dispose)) !=0)
		fprintf(stderr,"ql: Can't execute command: %s\n",dispose);

	return(status);

}
