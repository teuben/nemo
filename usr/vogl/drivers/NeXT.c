/*
 * VOGL/VOGLE driver for NeXTStep.
 * shamelessly hacked from X11 driver 9/91 by
 *
 * Patrick J. Flynn
 * School of EECS
 * Washington State University
 * Pullman, WA 99164-2752 USA cha cha cha
 * flynn@eecs.wsu.edu
 *
 * modified 9/22/91 to use NeXToid defaults (owner: "VOGL") for
 *
 * + preferred small & large font names/styles/sizes:
 *    `SmallFont' `LargeFont'
 *
 *   If no defaults are present, SmallFont is "Ohlfs 9.0" and
 *                               LargeFont is "Ohlfs 18.0"
 *
 *   The space in the font string is significant (see NeXT_font below)
 *
 * + window size & placement: `WindowX' `WindowY'
 *                            `WindowW' `WindowH'
 * 
 * If no defaults are present, the window is anchored (LL) at (10,10)
 * and is 792 pixels wide & high.  If the user has set the preferred
 * position & size using the VOGL prefXXX() calls, they override the
 * defaults.
 *
 * THIS CODE HAS ONLY BEEN TESTED ON A DOUBLE-HEADED CUBE (2 bit + ND)
 * There are a couple of color-screen-centrisms and 24-bit-centrisms in
 * the code, but NextStep is pretty good about drawing as good as it
 * can on the available hardware.
 *
 */

/*
 * set VOGLE if this is really for the VOGLE library.
#define VOGLE 1
 */

#include <stdio.h>
#include <math.h>
#include <strings.h> /* strfoo() */

#ifdef VOGLE
#include "vogle.h"
#else
#include "vogl.h"
#endif

#import <appkit/Application.h> /* [NXApp colorScreen] */
#import <appkit/Window.h>
#import <appkit/View.h>
#import <appkit/Font.h>
#import <dpsclient/wraps.h>
#import <dpsclient/event.h>
#import <appkit/color.h>
#import <appkit/graphics.h>
#import <appkit/nextstd.h>
#import <appkit/NXImage.h>     /* double buffer */
#import <appkit/afm.h>         /* hack to get metrics */
#import <appkit/defaults.h>

#define	CMAPSIZE	4096   /* perhaps we'll go RGBmode someday. */

Window *winder;
View *view;
Font *phont;
NXFontMetrics *fm;
NXImage *backbuf = nil;
NXColor colormap[CMAPSIZE]; /* yeah, I know it's ugly. */
int back_used;

id drawable; /* the thing we're drawing in (a View or an NXImage) */

/*
 * NeXT_init()
 * initializes the NeXT display.
 */
NeXT_init()
{
  int x0,y0,xs,ys;
  NXRect r;
  char name[80];
  const char *ptr;

  NXApp = [Application new];
  getprefposandsize(&x0,&y0,&xs,&ys);

  if (x0<0) {
#ifdef VOGLE
    ptr=NXReadDefault("VOGLE","WindowX");
#else
    ptr=NXReadDefault("VOGL","WindowX");
#endif
    if (!ptr)
      x0=10;
    else
      if (sscanf(ptr,"%d",&x0) != 1) x0=10;
    };

  if (y0<0) {
#ifdef VOGLE
    ptr=NXReadDefault("VOGLE","WindowY");
#else
    ptr=NXReadDefault("VOGL","WindowY");
#endif
    if (!ptr)
      y0=10;
    else 
      if (sscanf(ptr,"%d",&y0) != 1) y0=10;
    };
  
  if (xs<0) {
#ifdef VOGLE
    ptr=NXReadDefault("VOGLE","WindowW");
#else
    ptr=NXReadDefault("VOGL","WindowW");
#endif
    if (!ptr)
      xs=792;
    else
      if (sscanf(ptr,"%d",&xs) != 1) xs=792;
    };
  
  if (ys<0) {
#ifdef VOGLE
    ptr=NXReadDefault("VOGLE","WindowH");
#else
    ptr=NXReadDefault("VOGL","WindowH");
#endif
    if (!ptr)
      ys=792;
    else
      if (sscanf(ptr,"%d",&ys) != 1) ys=792;
    };
  
  NXSetRect(&r,x0,y0,xs,ys);
  winder = [[Window alloc] initContent:&r
                                 style:NX_TITLEDSTYLE
                               backing:NX_RETAINED
                            buttonMask:0
                                 defer:NO
                                screen:[NXApp colorScreen]];
  [winder setDepthLimit:NX_TwentyFourBitRGBDepth];
  vdevice.depth = 24;
  colormap[0]=NX_COLORBLACK;
  colormap[1]=NX_COLORRED;
  colormap[2]=NX_COLORGREEN;
  colormap[3]=NX_COLORYELLOW;
  colormap[4]=NX_COLORBLUE;
  colormap[5]=NX_COLORMAGENTA;
  colormap[6]=NX_COLORCYAN;
  colormap[7]=NX_COLORWHITE;

#ifdef VOGLE
    sprintf(name, "vogle %d", getpid());
    [winder setTitle:name];
#else
  if (!vdevice.wintitle) {
    sprintf(name, "vogl %d", getpid());
    [winder setTitle:name];
    }
  else {
    [winder setTitle:vdevice.wintitle];
    };
#endif

  drawable = view = [winder contentView];

  phont = nil;

  vdevice.sizeX = vdevice.sizeY = MIN(xs,ys) - 1;
  vdevice.sizeSx = xs;
  vdevice.sizeSy = ys;
  back_used = 0;
  [[winder makeKeyAndOrderFront:0] display];
  [drawable lockFocus];
  return(1);
}

/*
 * NeXT_exit
 *
 *	cleans up before returning the window to normal.
 */
NeXT_exit()
{
  [winder free];
  if (backbuf) [backbuf free];
  [NXApp free];
  return(1);
}

/*
 * NeXT_draw
 *
 *	draws a line from the current graphics position to (x, y).
 *
 */
NeXT_draw(x, y)
  int	x, y;
{
  PSsetlinewidth(0.0);
  PSnewpath();
  PSmoveto(vdevice.cpVx,vdevice.cpVy);
  PSlineto(x,y);
  PSstroke();
  if (vdevice.sync)
	  NXPing();
}

/*
 * NeXT_getkey
 *
 *	grab a character from the keyboard - blocks until one is there.
 */
int
NeXT_getkey()
{
  NXEvent e;
  NXEvent *ee = NXGetOrPeekEvent(DPSGetCurrentContext(),&e,
                                 NX_KEYDOWNMASK,NX_FOREVER,NX_BASETHRESHOLD,0);
  if (ee == NULL) fprintf(stderr,"error: NeXT_getKey failed.\n");
  return ee->data.key.charCode;
}
 
/*
 * NeXT_checkkey
 *
 *	Check if there has been a keyboard key pressed.
 *	and return it if so.
 */
int
NeXT_checkkey()
{
  NXEvent e;
  NXEvent *ee = NXGetOrPeekEvent(DPSGetCurrentContext(),&e,
                                 NX_KEYDOWNMASK,0.0,NX_BASETHRESHOLD,0);
  if (ee == NULL) return 0;
  return ee->data.key.charCode;
}

/*
 * NeXT_locator
 *
 * return the window location of the cursor, plus which mouse button,
 * if any, is been pressed.
 * 
 * the right mouse button is the LSB.  There is no middle button on a NeXT
 * so the 2LSB is always 0.  The left mouse button is the 3LSB.
 */
int
NeXT_locator(wx, wy)
	int	*wx, *wy;
{
  int msk=0,flg;
  NXPoint p;
  [winder getMouseLocation:&p];
  [view convertPoint:&p fromView:nil];
  *wx = p.x;
  *wy = p.y;
  PSbuttondown(&flg);
  if (flg) msk=4;
  PSrightbuttondown(&flg);
  if (flg) msk+=1;
  return msk;
}

#ifdef VOGLE
/*
 * NeXT_clear
 *
 * Clear the screen (or current buffer )to current colour
 */
NeXT_clear()
{
  NXRect r;
  [view getBounds:&r];
  NXRectFill(&r);
  if (vdevice.sync)
	  NXPing();
}

#else

/*
 * NeXT_clear
 *
 * Clear the viewport (or current buffer )to current colour
 */
NeXT_clear()
{
  NXRect r;
  float w=vdevice.maxVx - vdevice.minVx;
  float h=vdevice.maxVy - vdevice.minVy;
  NXSetRect(&r,vdevice.minVx,vdevice.minVy,w,h);
  [view getBounds:&r];
  NXRectFill(&r);
  if (vdevice.sync)
	  NXPing();
}

#endif

/*
 * NeXT_color
 *
 *	set the current drawing color index.
 */
NeXT_color(ind)
        int	ind;
{
  NXSetColor(colormap[ind]);
}

/*
 * NeXT_mapcolor
 *
 *	change index i in the color map to the appropriate r, g, b, value.
 */
NeXT_mapcolor(i, r, g, b)
	int	i;
	int	r, g, b;
{
  if (i >= CMAPSIZE) return(-1);
  colormap[i]=NXConvertRGBToColor(1.0*r/255.0,1.0*g/255.0,1.0*b/255.0);
}

/*
 * NeXT_font
 *
 *   Set up a hardware font. Return 1 on success 0 otherwise.
 *
 * This is system-dependent.  I assume that the fontfile parameter
 * has the font family name followed by a blank, followed by the size 
 * in points, e.g. "Ohlfs 384.7", "Helvetica-BoldOblique 1.0".
 * Note that the size can be floating-point.
 *
 */
char *strdup(const char *c)  /* ugly blech yech */
{
  char *d=malloc(1+strlen(c));
  bcopy(c,d,1+strlen(c));
  return d;
}

/* input: NS = "Name size", i.e. "Ohlfs 32.9".
 * output: *s = 32.9, return value = "Ohlfs";
 * returns NULL if string doesn't have both name and size.
 */
char *_getFontNameNSize(const char *ns,float *s)
{
  char *p,*q;
  if (!ns) return 0;
  p=strdup(ns); /* barf */
  q=index(p,' ');
  if (!q) return 0;
  *q++ = '\0'; /* null-terminate name */
  if (sscanf(q,"%f",s) != 1) { free(p); return 0;};
  return p;
}
  
  
NeXT_font(fontfile)
        char	*fontfile;
{
  char *name;
  float size;
  Font *newfont;
  if (!strcmp(fontfile,"small")) {
#ifdef VOGLE
    name=_getFontNameNSize(NXReadDefault("VOGLE","SmallFont"),&size);
#else
    name=_getFontNameNSize(NXReadDefault("VOGL","SmallFont"),&size);
#endif
    if (!name) { name=strdup("Ohlfs"); size=9.0; };
    }
  else if (!strcmp(fontfile,"large")) {
#ifdef VOGLE
    name=_getFontNameNSize(NXReadDefault("VOGLE","LargeFont"),&size);
#else
    name=_getFontNameNSize(NXReadDefault("VOGL","LargeFont"),&size);
#endif
    if (!name) { name=strdup("Ohlfs"); size=18.0; };
    }
  else {
    name=_getFontNameNSize((const char *)fontfile,&size);
    if (!name) { name=strdup("Ohlfs"); size=9.0; };
    };
  newfont = [Font newFont:name size:size matrix:NX_IDENTITYMATRIX];
  free(name);
  if (newfont) {
    NXFontMetrics *fm = [newfont metrics];
    int i;
    float wmax=0.0;
    [phont free];
    phont = newfont;
    [phont set];

    /* set hheight and hwidth for hardware fonts.  The scale factor
       is apparently equal to the point size.  It isn't well-documented. */

    vdevice.hheight=size*(fm->ascender-fm->descender);

    /* with this li'l loop we can handle var-width fonts but they look
       ugly. Stick to Courier & Ohlfs if you can. */
    for(i=0;i<MIN(fm->widthsLength,256);i++) {
      float w=fm->widths[i];
      if (w>wmax) wmax=w;
      };
    vdevice.hwidth=size*wmax;
    return 1;
    }
  else
    return 0;
}

/* 
 * NeXT_char (outputs one char)
 */
NeXT_char(c)
	char	c;
{
  char	s[2];

  s[0] = c; s[1]='\0';
  PSmoveto(vdevice.cpVx,vdevice.cpVy);
  PSshow(s);
  if (vdevice.sync)
	  NXPing();
}

/*
 * NeXT_string
 *
 *	Display a string at the current drawing position.
 */
NeXT_string(s)
        char	s[];
{
  PSmoveto(vdevice.cpVx,vdevice.cpVy);
  PSshow(s);
  NXPing();
}

/*
 * NeXT_fill
 *
 *	fill a polygon
 */
NeXT_fill(n, x, y)
	int	n, x[], y[];
{
  int	i;
  PSnewpath();
  PSmoveto(x[0],y[0]);
  for(i=1;i<n;i++) PSlineto(x[i],y[i]);
  PSclosepath();
  PSfill();
  if (vdevice.sync)
	  NXPing();
}

/*
 * NeXT_backbuf
 *
 *	Set up double buffering by allocating the back buffer and
 *	setting drawing into it.
 */
NeXT_backbuf()
{
  NXRect r;
  [drawable unlockFocus];
  [view getBounds:&r];
  backbuf = [[NXImage alloc] initSize:&r.size];
  if (![backbuf useCacheWithDepth:NX_TwentyFourBitRGBDepth]) {
    fprintf(stderr,"couldn't create backing buffer.\n");
    return 0;
    };
  drawable = backbuf;
  [drawable lockFocus];
  back_used = 1;
  return(1);
}

/*
 * NeXT_swapbuf
 *
 *	Swap the back and front buffers. (Really, just copy the
 *	back buffer to the screen).
 */
NeXT_swapbuf()
{
  NXPoint origin = {0,0};
  [drawable unlockFocus];
  [view lockFocus];
  [backbuf composite:NX_COPY toPoint:&origin];
  [winder flushWindow];
  [view unlockFocus];
  [drawable lockFocus];
}

/*
 * NeXT_frontbuf
 *
 *	Make sure we draw to the screen.
 */
NeXT_frontbuf()
{
  [drawable unlockFocus];
  drawable = view;
  [drawable lockFocus];
}

/*
 * Syncronise the display with what we think has been sent to it...
 */
NeXT_sync()
{
	NXPing();
}

NeXT_setls(ls)
	int	ls;
{
}

NeXT_setlw(w)
	int	w;
{
}

/*
 * the device entry
 */
static DevEntry NeXTdev = {
	"NeXT",
	"small",
	"large",
	NeXT_backbuf,
	NeXT_char,
	NeXT_checkkey,
	NeXT_clear,
	NeXT_color,
	NeXT_draw,
	NeXT_exit,
	NeXT_fill,
	NeXT_font,
	NeXT_frontbuf,
	NeXT_getkey,
	NeXT_init,
	NeXT_locator,
	NeXT_mapcolor,
	NeXT_setls,
	NeXT_setlw,
	NeXT_string,
	NeXT_swapbuf, 
	NeXT_sync
};

/*
 * _NeXT_devcpy
 *
 *	copy the NeXT device into vdevice.dev.
 */
_NeXT_devcpy()
{
	vdevice.dev = NeXTdev;
}
