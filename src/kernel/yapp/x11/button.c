#include "config.h"
#include <X11/Xlib.h>
#include "event.h"
/***   Button handling --- very simple, but working ***/

/*
 *  Global funtions:
 *
 *          makebutton(xpos,ypos, width,height,text,eventfunction)
 *                     -   creates a button
 *          drawbuttons() - display buttons on screen
 *          buttonloop() - waits for button events
 *          exitbuttonloop() - it should be clear :-)
 *
 *
 */

#define MAXBUTTONS 50

typedef struct {
    real x;
    real y;
    real w;
    real h;
    char text[64];
    void (*func)();
} buttonrec, *buttonrecptr;

static buttonrec x_buttons[MAXBUTTONS];

static int buttoncounter = 0;
static bool breakloop = FALSE;

extern real x_xmin, x_ymax;
extern real x_xfactor, x_yfactor;

#define MAPX(x) ( (real)(x)/x_xfactor + x_xmin )
#define MAPY(y) ( x_ymax - (real)(y)/x_yfactor )

#define INBOX(b,x,y) ( (x>((b).x-(b).w)) && (y>((b).y-(b).h)) && (x<((b).x+(b).w)) && (y < ((b).y+(b).h)) )

static void lookup(ix, iy)
int ix,iy;
{
    int i;
    real x,y;
    void flicker();
    
    x = MAPX(ix);
    y = MAPY(iy);

    for(i = 0; i < buttoncounter; i++)
    {
	if(INBOX(x_buttons[i],x,y))
	{
	    if(x_buttons[i].func != NULL)
	    {
		x_buttons[i].func();
	    }
	    break;
	}	
    }
}

static void drawbutton(b)
buttonrecptr b;
{
    plmove(b->x - b->w, b->y - b->h);
    plline(b->x + b->w, b->y - b->h);
    plline(b->x + b->w, b->y + b->h);
    plline(b->x - b->w, b->y + b->h);
    plline(b->x - b->w, b->y - b->h);

    pljust(0);
    pltext(b->text, b->x, b->y, b->h*6.0 , 0.0);
    pljust(-1);    
}

makebutton(bx, by, bw, bh, btext, bfunc)
real bx, by, bw, bh;
string btext;
void (*bfunc)();
{
    int b;

    if(buttoncounter == MAXBUTTONS)
	error("cannot handle more than %d buttons (sorry!)\n",
	      MAXBUTTONS);
    
    bw /= 2.; bh /= 2.;

    b=buttoncounter;
    x_buttons[b].x = bx;
    x_buttons[b].y = by;
    x_buttons[b].w = bw;
    x_buttons[b].h = bh;
    strcpy(x_buttons[b].text, btext);
    x_buttons[b].func = bfunc;
    
    buttoncounter++;
}

drawbuttons()
{
    buttonrecptr bp;

    for(bp = x_buttons; bp < x_buttons + buttoncounter; bp++)
	drawbutton(bp);    
}

buttonloop()
{
    XEvent ev;
    XButtonEvent *bev;

    breakloop = FALSE;
    while(!breakloop)
    {
	if (x_ask_event(&ev))
	{
	    if(ev.type == ButtonPress)
	    {
		bev = (XButtonEvent *) &ev;
		lookup(bev->x, bev->y);
	    }
	}
    }
}

exitbuttonloop()
{
    breakloop = TRUE;
}
