/* movie  -  run a FITS movie from individual images */

#define WEDGE_HEIGHT 256
/*#define VERSION 3.0*/
#include "ql.h"
#include <suntool/panel.h>
#include <suntool/icon.h>
#include <suntool/textsw.h>
#include <sys/file.h>

#define MAXPIC 200
#define WINWIDTH 1025
#define WINHEIGHT 851

static void RunMovie();
static void PrevFrame();
static void NextFrame();
static void quit_proc();
static void action_set();
static void choose_color();
void create_help_popup();
void done_proc(); 
static void MouseEvent();
static void Slid1Proc();
static void Slid2Proc();
static Panel_setting MinGet();
static Panel_setting MaxGet();
static Panel_setting GamGet();
static Panel_setting PseudoGet();
void repaint();
static void Dither(/*R-image,G-image,B-image*/);
static void reverse_1(); 
static void itoa();

static short icon_image[] = {
#include "movie.icon"
};
mpr_static(ii,64,64,1,icon_image); /* memory pixrect containing prog icon */

static short    dat_mouse_left[] = {
#include <images/confirm_left.pr>
};
mpr_static(mouse_left, 16, 16, 1, dat_mouse_left); /* mouse icons */

static short    dat_mouse_middle[] = {
#include <images/confirm_middle.pr>
};
mpr_static(mouse_middle, 16, 16, 1, dat_mouse_middle);

static short    dat_mouse_right[] = {
#include <images/confirm_right.pr>
};
mpr_static(mouse_right, 16, 16, 1, dat_mouse_right);

Frame frame, help_frame;
Canvas canvas;
Panel panel, help_panel; 
Textsw textsw;
Icon icon;
Pixfont *pf;
PICTURE pic[MAXPIC];
PICTURE pic24;
int is24Bit;

char tab[2]  = "\t"; /* these are used in MinGet and MaxGet */
char null[2] = "\0";
char pt[2]=".";      /* this used for filename assembly      */

/* These map the available 3 or 2 bits into the 8-bit values and try
 * to correct for the non-linearity of the display. They AREN'T correct yet.
 */
char colorLevel3[8] = { 0, 76, 102, 127, 152, 178, 203, 255};
char colorLevel2[8] = { 0, 110, 170, 255 };

static unsigned char rgb_red[256], rgb_green[256],
	      rgb_blue[256], rgb_black[256];

int curPic = 0;
int zoom = 1, xsize, ysize;	/* linear zoom, size of full pixrect 		*/
static int npic;	        /* number of images to load          		*/
unsigned usecs=0;               /* microseconds to sleep between frames,
				 * changed by slider 				*/
int movie_repeat=1;             /* number of times to cycle movie, changed
				 * by slider 					*/
static int clicks=1;  		/* fix double clicking bug--a kluge! 		
				 * image frame interprets a single mouse click 
				 * as a  double click for some reason.		*/
float maxv = 0, minv = 0;	/* image pixel range.  If set to zero,
				 * DATAMAX and DATAMIN from the FITS
				 * header will be used 				*/
int action_flag=0;		/* toggle for Min/Max Action:
				 * 0 = change current image
				 * 1 = change all images
				 */
int display_type = 0;		/* Cycle choice between grey scale (0) and
				 * pseudocolor (1) and 3-color (2) 		*/
int nsteps	= 128;		/* default contour level for pseudocolor 	*/
double gammaCor = .37; 	        /* actually 1/gamma) 				*/
int Dithered=0;  	/* flag turned on when using dithering on 8-bit display */

Panel_item  pseudo_item, wait_item;

main(argc, argv)
int argc;
char **argv;
{
	register Pixwin *pw;
	Pixrect *pr;
	Panel_item min_item, max_item, gam_item;
	void show_help(); 

	int nextPic, start, i;
	int winwidth = WINWIDTH, winheight = WINHEIGHT;
	char fs[30]; 
	char fn[40];
	char s[4];
	char *cp;

	fprintf(stdout,"This is movie V%2.1f\n\n",VERSION);
	if(argc < 5) {
		fprintf(stderr,"\nusage: movie npic i zoom filespec [min max]\n");
		fprintf(stderr," min and max are optional.\nfilenames will be assembled thusly:\n");
		error("filename = filespec + '.' + 'n', n=i to npic");
	}
	if(argc == 7) {
		minv=atof(argv[5]);
		maxv=atof(argv[6]);
	}
	if(minv==maxv) 
	  argc=5;
        
	if( (npic = atoi(argv[1])) > MAXPIC)
	        error("too many files. %d maximum",MAXPIC);
	zoom = atoi(argv[3]);
        start = atoi(argv[2]);
/*        fprintf(stderr, "\nzoom=%d  start=%d\n",zoom,start);*/
	if(zoom < 1 || zoom > 16)
		error("bad zoom = %d", zoom);

        sscanf(argv[4],"%s",fs);	
	strcat(fs,pt);

	pf = pf_default(); /*use default font*/
	;
	icon = icon_create(ICON_IMAGE, &ii, 0);

/* load 1st image */
         nextPic=0;
	 strcpy(fn,fs);
	 itoa(start,s);
 	 strcat(fn,s);
	 strcpy(pic[nextPic].fn,fn);
/* exit prog if first file not accessible */
          if( (access(fn,F_OK)!=0)||(access(fn,R_OK)!=0) ) {
	    fprintf(stderr,"\07\nFirst file %s not accessible.\nHave you got the filespec right?\n",fn);
            exit(1);
	  }		

	 is24Bit = Is24Bit();
	 /*if( !(is24Bit = Is24Bit()) )
		colorMapSize = 128;*/

	 LoadImage(&pic[nextPic], 0, GREY);
	 pr = pic[curPic=nextPic].pr;  /* nextPic or curPic alone DON'T WORK! */
	 xsize = pr->pr_size.x;
	 ysize = pr->pr_size.y;

	 if(imageHdr[0]->title[0])
		cp = imageHdr[0]->title;
	 else
		cp = pic[curPic].fn;
	 strcpy(pic[curPic].label, cp);
	 if(is24Bit) {
		pic24.pr = CreatePr24();
	 	ConvertToPr24(&pic[curPic],0);
	 }
/* make frame etc */
	 frame = window_create(0,FRAME,
		FRAME_LABEL, cp,
		FRAME_ARGS, argc, argv,
		FRAME_ICON, icon,
/* 		FRAME_PROPS_ACTION_PROC, show_help, */
/*		FRAME_PROPS_ACTIVE,  TRUE, */
/*what do these do??*/
		0);
	 canvas = window_create(frame, CANVAS,
		CANVAS_AUTO_SHRINK, FALSE,
		CANVAS_WIDTH, xsize,
		CANVAS_HEIGHT, ysize,
		CANVAS_RETAINED, FALSE,
		CANVAS_REPAINT_PROC, repaint,
		WIN_EVENT_PROC, MouseEvent,
		0);
	    if(is24Bit)
		window_set(canvas, CANVAS_COLOR24, TRUE, 0);
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

	   create_help_popup(); 

/* create all the fancy panel stuff */

	   panel = window_create(frame, PANEL,
/*		WIN_EVENT_PROC, MouseEvent, */ /* this screws up buttons */
		WIN_COLUMNS, 20, WIN_ROWS, 25, 0);

	   panel_create_item(panel,PANEL_MESSAGE,
		PANEL_LABEL_IMAGE, &mouse_left,
		PANEL_ITEM_X, ATTR_COL(1),
		PANEL_ITEM_Y, ATTR_COL(1),
		0 );
   	   panel_create_item(panel,PANEL_MESSAGE,
		PANEL_LABEL_STRING, "Previous Frame",
		PANEL_ITEM_X, ATTR_COL(4),
		PANEL_ITEM_Y, ATTR_COL(1),
		0 ); 
	   panel_create_item(panel,PANEL_MESSAGE,
		PANEL_LABEL_IMAGE, &mouse_middle,
		PANEL_ITEM_X, ATTR_COL(1),
		PANEL_ITEM_Y, ATTR_COL(2),
		0 );
	   panel_create_item(panel,PANEL_MESSAGE,
		PANEL_LABEL_STRING, "Next Frame",
		PANEL_ITEM_X, ATTR_COL(4),
		PANEL_ITEM_Y, ATTR_COL(2),
		0 ); 
	   panel_create_item(panel,PANEL_MESSAGE,
		PANEL_LABEL_IMAGE, &mouse_right,
		PANEL_ITEM_X, ATTR_COL(1),
		PANEL_ITEM_Y, ATTR_COL(3),
		0 );
	   panel_create_item(panel,PANEL_MESSAGE,
		PANEL_LABEL_STRING, "Run Movie",
		PANEL_ITEM_X, ATTR_COL(4),
		PANEL_ITEM_Y, ATTR_COL(3),
		0 ); 
	   panel_create_item(panel,PANEL_MESSAGE,
		PANEL_LABEL_STRING, "Frame Speed: ",
		PANEL_ITEM_X, ATTR_COL(1),
		PANEL_ITEM_Y, ATTR_COL(5),
		0 ); 
	   panel_create_item(panel, PANEL_SLIDER, 
		PANEL_VALUE, 		0,
		PANEL_MIN_VALUE,	0,
		PANEL_MAX_VALUE,	10,
		PANEL_SLIDER_WIDTH,	60,
		PANEL_ITEM_X,		ATTR_COL(1),
		PANEL_ITEM_Y,		ATTR_COL(6),
		PANEL_NOTIFY_PROC,	Slid1Proc,
		0);
	   panel_create_item(panel,PANEL_MESSAGE,
		PANEL_LABEL_STRING, "Movie Cycles: ",
		PANEL_ITEM_X, ATTR_COL(1),
		PANEL_ITEM_Y, ATTR_COL(7),
		0 ); 
	   panel_create_item(panel, PANEL_SLIDER,
		PANEL_VALUE, 		1,
		PANEL_MIN_VALUE,	1,
		PANEL_MAX_VALUE,	10,
		PANEL_SLIDER_WIDTH,	60,
		PANEL_ITEM_X,		ATTR_COL(1),
		PANEL_ITEM_Y,		ATTR_COL(8),
		PANEL_NOTIFY_PROC,	Slid2Proc,
		0);
	   panel_create_item(panel,PANEL_MESSAGE,
		PANEL_LABEL_STRING, "Min/Max Action:",
		PANEL_ITEM_X, ATTR_COL(1),
		PANEL_ITEM_Y, ATTR_COL(10),
		0 ); 
	   panel_create_item(panel,PANEL_CYCLE,
		PANEL_ITEM_X, 		ATTR_COL(1),
		PANEL_ITEM_Y, 		ATTR_COL(11),
/*		PANEL_LABEL_STRING, 	"Min/Max Action:",*/
		PANEL_CHOICE_STRINGS, 	"Change current image ", 
                			"Change all images", 
					0,
		PANEL_DISPLAY_LEVEL, 	PANEL_CURRENT,
		PANEL_NOTIFY_PROC, 	action_set,
		0);
	   min_item=panel_create_item(panel, PANEL_TEXT,
		PANEL_LABEL_STRING,  	"Min: ",
		PANEL_VALUE_DISPLAY_LENGTH, 10, 
		PANEL_ITEM_X,		ATTR_COL(1),
		PANEL_ITEM_Y,		ATTR_COL(12),
		PANEL_NOTIFY_PROC,	MinGet,
		0);
	   max_item=panel_create_item(panel, PANEL_TEXT,
		PANEL_LABEL_STRING,  	"Max: ",
		PANEL_VALUE_DISPLAY_LENGTH, 10, 
		PANEL_ITEM_X,		ATTR_COL(1),
		PANEL_ITEM_Y,		ATTR_COL(13),
		PANEL_NOTIFY_PROC,	MaxGet,
		0);
	   gam_item=panel_create_item(panel, PANEL_TEXT,
		PANEL_LABEL_STRING,  	"Gamma: ",
		PANEL_VALUE_DISPLAY_LENGTH, 10, 
		PANEL_ITEM_X,		ATTR_COL(1),
		PANEL_ITEM_Y,		ATTR_COL(14),
		PANEL_NOTIFY_PROC,	GamGet,
		0);
	   pseudo_item=panel_create_item(panel, PANEL_TEXT,
		PANEL_LABEL_STRING,  	"Number of Contours: ",
		PANEL_VALUE_DISPLAY_LENGTH, 5, 
		PANEL_ITEM_X,		ATTR_COL(1),
		PANEL_ITEM_Y,		ATTR_COL(15),
		PANEL_NOTIFY_PROC,	PseudoGet,
		PANEL_SHOW_ITEM,	FALSE,
		0);
	   wait_item=panel_create_item(panel, PANEL_TEXT,
		PANEL_LABEL_STRING,  	"Please wait...",
		PANEL_VALUE_DISPLAY_LENGTH, 1, 
		PANEL_ITEM_X,		ATTR_COL(1),
		PANEL_ITEM_Y,		ATTR_COL(24),
		PANEL_SHOW_ITEM,	FALSE,
		0);
	   panel_create_item(panel,PANEL_CYCLE,
		PANEL_ITEM_X, 		ATTR_COL(1),
		PANEL_ITEM_Y, 		ATTR_COL(17),
		PANEL_CHOICE_STRINGS, 	"Grey scale", 
					"Pseudo-color",
                			"3-Color", 
					0,
		PANEL_DISPLAY_LEVEL, 	PANEL_CURRENT,
		PANEL_NOTIFY_PROC, 	choose_color,
		0);
	   panel_create_item(panel,PANEL_BUTTON,
		PANEL_LABEL_IMAGE,
		panel_button_image(panel," Help ",0,0),
		PANEL_ITEM_X, ATTR_COL(1),
		PANEL_ITEM_Y, ATTR_COL(19),
		PANEL_NOTIFY_PROC, show_help,
		0 ); 
	   panel_create_item(panel,PANEL_BUTTON,
		PANEL_LABEL_IMAGE,
		panel_button_image(panel," Quit ",0,0),
		PANEL_ITEM_X, ATTR_COL(1),
		PANEL_ITEM_Y, ATTR_COL(20),
		PANEL_NOTIFY_PROC, quit_proc,
		0 ); 
	   panel_create_item(panel,PANEL_CYCLE,
		PANEL_ITEM_X, 		ATTR_COL(1),
		PANEL_ITEM_Y, 		ATTR_COL(17),
		PANEL_CHOICE_STRINGS, 	"Grey scale", 
					"Pseudo-color",
                			"3-Color", 
					0,
		PANEL_DISPLAY_LEVEL, 	PANEL_CURRENT,
		PANEL_NOTIFY_PROC, 	choose_color,
		0);
		
	  pw = canvas_pixwin(canvas);
	  window_fit_width(panel);
	  window_fit(frame);

/* set up grayscale "colormap" and colormap for dithering */

       if(!is24Bit) {
		pw_setcmsname(pw, "halfgrey");
		for(i=0; i < colorMapSize; i++) {
			rgb_red[i] =   colorLevel3[i >> 5];
			rgb_green[i] = colorLevel3[(i >> 2) & 7];
			rgb_blue[i] =  colorLevel2[(i & 3)];
			rgb_black[i] = 0; 
		} 
       }
       gammacorrect(pw,gammaCor);
       ShowCurrent(pw);

       StartDisplay();

/* load up rest of images and display each in turn */


       for(nextPic=1;nextPic<npic;nextPic++) {
	if(argc != 7) {
		minv = maxv = 0;
        }
	strcpy(fn,fs);
	itoa(nextPic+start,s);
 	strcat(fn,s);
	strcpy(pic[nextPic].fn,fn);

/* this avoids core dumps when npic is wrong */
/* check file accessibility by access(fn) */
/* if they don't exist spit out error, set npic=nextPic, break loop and continue*/

          if( (access(fn,F_OK)!=0)||(access(fn,R_OK)!=0) ) {
	    npic=nextPic;	
	    fprintf(stderr,"\07\nfile %s not accessible.\nUsing only images already loaded\n",fn);
            break;
	  }		
	  LoadImage(&pic[nextPic], 0, GREY);
	  if(imageHdr[0]->title[0])
		cp = imageHdr[0]->title;
	  else
		cp = pic[nextPic].fn;
	  strcpy(pic[nextPic].label, cp);
	  xsize = pr->pr_size.x;
	  if(is24Bit)
		ConvertToPr24(&pic[curPic],0);
	  curPic = nextPic;
	  ShowCurrent(pw);
       } 

/* make up a extra pic for possible dithered image but don't display*/
       if(!is24Bit) {
		strcpy(pic[npic].fn,pic[0].fn);
		LoadImage(&pic[npic], 0, GREY);
	 	pr = pic[npic].pr;  
       }

       notify_no_dispatch();
       window_set(frame,WIN_SHOW,FALSE,0); /*unsatisfactory: causes blink*/
       window_main_loop(frame);
}

static void repaint(canvas, pw, repaint_area)
Canvas canvas;
Pixwin *pw;
Rectlist *repaint_area;
{
	ShowCurrent(pw);
}

static void MouseEvent(canvas, event)
Canvas canvas;
Event *event;
{
	switch (event_action(event)) {
	    case MS_LEFT:
		PrevFrame();
		break;
	    case MS_MIDDLE:
		NextFrame();
		break;
	    case MS_RIGHT:
		RunMovie();
		break;
	}
}	

static void Slid1Proc(item,value,event)
Panel_item item;
int value;
Event *event;
{
 	usecs = value*100000; 
}

static void Slid2Proc(item,value,event)
Panel_item item;
int value;
Event *event;
{
	movie_repeat = value;
}


static Panel_setting MinGet(min_item, event)
Event *event;
Panel_item min_item;
{
	int badstring;
	char minbuf[10];

	badstring=0;
	strcpy(minbuf, (char*)panel_get_value(min_item));

	if ( ((strcmp(minbuf,null)) == 0) || ((strcmp(minbuf,tab)) == 0) ) 
		badstring=1;

/* if badstring, leave minv unchanged */
	if (!badstring)   
		minv=atof(minbuf);

	return(panel_text_notify(min_item, event));
}

static Panel_setting MaxGet(max_item, event)
Event *event;
Panel_item max_item;
{
	Pixwin *pw;
	char maxbuf[10];
	int badstring;

	badstring=0;
	pw=canvas_pixwin(canvas);
	strcpy(maxbuf, (char*)panel_get_value(max_item));

	if ( ((strcmp(maxbuf,null)) == 0) || ((strcmp(maxbuf,tab)) == 0) ) 
		badstring=1;

/* if badstring, leave maxv unchanged and don't bother reloading the pic */

	if (!badstring)  {  
	  maxv=atof(maxbuf); 
	  if(action_flag==1) {
	        panel_set(wait_item, PANEL_SHOW_ITEM, TRUE, 0);
		for(curPic=0;curPic<npic;curPic++) {
	  		LoadImage(&pic[curPic], 0, GREY);
			if(is24Bit)
	 			ConvertToPr24(&pic[curPic],0);
			if(Dithered && curPic==0)
				gammacorrect(pw,gammaCor);
			ShowCurrent(pw);
		}
	        panel_set(wait_item, PANEL_SHOW_ITEM, FALSE, 0);
		curPic--;
	  } else {
		LoadImage(&pic[curPic], 0, GREY);
		if(is24Bit)
	 		ConvertToPr24(&pic[curPic],0);
		if(Dithered)
			gammacorrect(pw,gammaCor);
		ShowCurrent(pw);
	  }
	}

	return(panel_text_notify(max_item, event));
}
static Panel_setting GamGet(gam_item, event)
Event *event;
Panel_item gam_item;
{
	Pixwin *pw;

	char gambuf[10];
	int badstring;

	badstring=0;
	pw=canvas_pixwin(canvas);
	strcpy(gambuf, (char*)panel_get_value(gam_item));

	if ( ((strcmp(gambuf,null)) == 0) || ((strcmp(gambuf,tab)) == 0) ) 
		badstring=1;

/* if badstring, leave gamma unchanged */

	if (!badstring)  {  
	  gammaCor = (double) atof(gambuf); 
	  gammacorrect(pw,gammaCor);
	}

	return(panel_text_notify(gam_item, event));
}
static void PrevFrame()  /*display previous frame*/ 
{
	Pixwin *pw;
	Pixrect *pr;
	int Red,Blue, cmap_switch;

	if(clicks == 2) { /* to fix double-clicking bug */
	   pw=canvas_pixwin(canvas);
	   if((--curPic) < 0) 
		curPic = npic-1;
	   if(is24Bit) {
		if(display_type == 2) {
			PICTURE *tPic = &pic[(curPic < 2)? curPic - 2 + npic: curPic - 2];
	 		ConvertToPr24(tPic, 5);
			strcpy(pic24.label, pic[curPic].label);
		} else 
			ConvertToPr24(&pic[curPic], 0);
	   }
	   if(!is24Bit && display_type==2) {
		Dithered=1;
		Blue = curPic-1;
		Red = curPic+1;
		cmap_switch=curPic;
		if(curPic == 0) 
			Blue = curPic;
		if(curPic == npic-1) {
			Red = curPic;
			cmap_switch=-1;
		}
		Dither(&Red,&curPic,&Blue);	
		pr = pic[npic].pr;
		strcpy(pic[npic].label,pic[curPic].label);
		window_set(frame,FRAME_LABEL,pic[npic].label,0);
		pw_rop(pw,0,0, pr->pr_size.x, pr->pr_size.y,PIX_SRC,pr,0,0);
		switch (cmap_switch) {
		   case 0:	
		     pw_putcolormap(pw,0,colorMapSize,rgb_black,rgb_green,rgb_blue);
		     break;
		   case -1:
		     pw_putcolormap(pw,0,colorMapSize,rgb_red,rgb_green,rgb_black);
		     break;
		   default:
		     pw_putcolormap(pw,0,colorMapSize,rgb_red,rgb_green,rgb_blue);
		     break;
		}
	   }
	   else 
		ShowCurrent(pw);
	  clicks=1;
	} else
	    clicks++;
}

static void NextFrame()  /*display next frame*/ 
{
	Pixwin *pw;
	Pixrect *pr;
	int Red,Blue,cmap_switch;

	if(clicks == 2) { /* to fix double-clicking bug */
	  pw=canvas_pixwin(canvas);
	  if( (++curPic) > (npic-1) ) 
		curPic = 0;
	  if(is24Bit) 
		ConvertToPr24(&pic[curPic], (display_type == 2)? 4: 0);
	  if(!is24Bit && display_type==2) {
		Dithered=1;
		Blue = curPic-1;
		Red = curPic+1;
		cmap_switch=curPic;
		if(curPic == 0)
			Blue = curPic;
		if(curPic == npic-1) {
			Red = curPic;
			cmap_switch=-1;
		}
		Dither(&Red,&curPic,&Blue);	
		pr = pic[npic].pr;
		strcpy(pic[npic].label,pic[curPic].label);
		window_set(frame,FRAME_LABEL,pic[npic].label,0);
		pw_rop(pw,0,0, pr->pr_size.x, pr->pr_size.y,PIX_SRC,pr,0,0);
		switch (curPic) {
		   case 0:	
		     pw_putcolormap(pw,0,colorMapSize,rgb_black,rgb_green,rgb_blue);
		     break;
		   case -1:
		     pw_putcolormap(pw,0,colorMapSize,rgb_red,rgb_green,rgb_black);
	  	     break;
		   default:
		     pw_putcolormap(pw,0,colorMapSize,rgb_red,rgb_green,rgb_blue);
		     break;
		}
	   }
	   else 
		ShowCurrent(pw);
	   clicks=1;
	} else
	   clicks++;
}

static void RunMovie() /* run movie movie_repeat times */
{
	int i,Red,Blue,cmap_switch;
	Pixwin *pw;
	Pixrect *pr;
	
	if(clicks == 2) { /* to fix double-clicking bug */
	   pw=canvas_pixwin(canvas);
	   for(i=0;i<movie_repeat;i++) {	
		for(curPic=0;curPic<npic;curPic++) {
	  	   if(is24Bit) 
			ConvertToPr24(&pic[curPic], (display_type == 2)? 4: 0);
		   if(!is24Bit && display_type==2) {
			Dithered=1;
			Blue = curPic-1;
			Red = curPic+1;
			cmap_switch=curPic;
			if(curPic == 0)
				Blue = curPic;
			if(curPic == npic-1) {
				Red = curPic;
				cmap_switch=-1;
			}
			Dither(&Red,&curPic,&Blue);	
			pr = pic[npic].pr;
			strcpy(pic[npic].label,pic[curPic].label);
			window_set(frame,FRAME_LABEL,pic[npic].label,0);
			pw_rop(pw,0,0, pr->pr_size.x, pr->pr_size.y,PIX_SRC,pr,0,0);
			switch (cmap_switch) {
		   	  case 0:	
		     	    pw_putcolormap(pw,0,colorMapSize,rgb_black,rgb_green,rgb_blue);
		     	    break;
		   	  case 1:
		     	    pw_putcolormap(pw,0,colorMapSize,rgb_red,rgb_green,rgb_blue);
		     	    break;
		   	  case -1:
		     	    pw_putcolormap(pw,0,colorMapSize,rgb_red,rgb_green,rgb_black);
		     	    break;
			}
		   }	
		   if((curPic > 1 || display_type != 2) && !Dithered) 
			ShowCurrent(pw);
		   if(usecs > 100)
			usleep(usecs);
		}
		curPic--;
	}
	clicks=1;
       }
       else
	clicks++;
}

static void quit_proc() 
{
	window_set(frame, FRAME_NO_CONFIRM, TRUE, 0);
	window_destroy(frame);
	exit(0);
}

/* convert integer n to character string s[]		*/
/* note use of do-while statement:            		*/
/* the string is generated backwards, then reversed. 	*/

static void itoa(n,s)
char s[];
int n;
{
	int i,sign;
	if((sign=n)<0) 			/* record sign */
		n= -n;
		i=0;
		do {           /* generate digits in reverse order */
			s[i++] = n % 10 + '0';  /* get next digit */
		} while ((n /= 10) > 0);	/* delete it      */
		if (sign < 0) {
			s[i++] = '-';
		}
		s[i] = '\0';
		reverse_1(s);
}

/* SUBROUTINE: reverse string s[] and put it in  back in s[]    */
/* call has form:					        */
/*			reverse_1(s);			        */
/* use subroutine reverse_2 if you don't want s[] overwritten   */

#define MAXCHAR 100
static void reverse_1(str_in)
char *str_in;
{
	int i,j;
	static char out[MAXCHAR];
	static char *str_out;

	/* first find length of str_in, not including the trailing \0 */
	    str_out = out;
	    i=0;
	    while(str_in[i] != '\0')
		 i++; 
	/* now stuff str_out[] with the reverse of str_in[] */
	    j=i-1;
	    i=0;
	    while(j != -1) {
	     str_out[i] = str_in[j];
	     j--;
	     i++;
	    }
	    str_out[i+1] = '\0';	 /* add trailing \0 to str_out[] */
	/* now overwrite str_in[] with str_out[] */
		strcpy(str_in,str_out);    
}

gammacorrect(pw, power)
Pixwin *pw;
double power;
{
	int i;
	double scale = 1. / (colorMapSize - 1);
	unsigned char map[256]; 

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

char *helpmsg[] = {
#include "movie_help.h"
};
void create_help_popup()
{
	void done_proc();

	help_frame = window_create(frame, FRAME, 
		     FRAME_LABEL, " Movie Help ",
		     FRAME_SHOW_LABEL, TRUE,
		     0);

	textsw = (Textsw)window_create(help_frame, TEXTSW,
			TEXTSW_BROWSING, TRUE,
			TEXTSW_CONTENTS, *helpmsg , 
			TEXTSW_DISABLE_LOAD, TRUE,
			0);

	help_panel = window_create(help_frame, PANEL, 0);

	panel_create_item(help_panel, PANEL_BUTTON,
			PANEL_LABEL_IMAGE,  
			panel_button_image(help_panel, "  Done  ", 0,0),
			PANEL_ITEM_X,  ATTR_COL(65),
			PANEL_ITEM_Y,  ATTR_ROW(0),
			PANEL_NOTIFY_PROC,  done_proc, 
			0);

	window_fit_height(help_panel); 
	window_fit(textsw);
	window_fit(help_frame);
}

void show_help()
{
	window_set(help_frame, WIN_SHOW, TRUE, 0);
}

void done_proc()
{
	window_set(help_frame, WIN_SHOW, FALSE, 0);
}

static void action_set(item,value,event)
Panel_item item;
int value;
Event *event;
{
	action_flag = value;
}

static void choose_color(item, value, event)
Panel_item item;
int value;
Event *event;
{
	Pixwin *pw;

	pw = canvas_pixwin(canvas);
	display_type = value;

	if(display_type==2 && !is24Bit) 
		Dithered=1;
	else
		Dithered=0;

	if(display_type==1) 
	    panel_set(pseudo_item, PANEL_SHOW_ITEM, TRUE, 0);
	else {
	    panel_set(pseudo_item, PANEL_SHOW_ITEM, FALSE, 0);
	    gammacorrect(pw,gammaCor);
	}
}

static Panel_setting PseudoGet(pseudo_item, event)
Event *event;
Panel_item pseudo_item;
{
	Pixwin *pw;
	char pseudobuf[10];
	int badstring;

	badstring=0;
	pw=canvas_pixwin(canvas);
	strcpy(pseudobuf, (char*)panel_get_value(pseudo_item));

	if ( ((strcmp(pseudobuf,null)) == 0) || ((strcmp(pseudobuf,tab)) == 0) ) 
			badstring=1;

	/* if badstring, leave nsteps unchanged */

	if (!badstring)  {  
		 nsteps = atoi(pseudobuf); 
		 pseudocolor(pw,gammaCor,nsteps);
	}

	return(panel_text_notify(pseudo_item, event));
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

/* This takes a red image in R, green in G, and blue in B and makes
 * a 3-color image in pic[npic] using dithering to avoid contouring.
 */
#define DITHERSIZE 4
static void Dither(R,G,B)
int *R,*G,*B;			/* image numbers for R,G,B planes  */
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
	int RR,GG,BB;
	RR=*R;
	GG=*G;
	BB=*B;
	rel = (unsigned char *)mpr_d(pic[RR].pr)->md_image;
	gel = (unsigned char *)mpr_d(pic[GG].pr)->md_image;
	bel = (unsigned char *)mpr_d(pic[BB].pr)->md_image;

	CheckSame(RR,BB); 
	rgbel = (unsigned char *)mpr_d(pic[npic].pr)->md_image;
	linebytes = mpr_mdlinebytes(pic[RR].pr);
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
/* go to next row */
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
