#include "options.h"
#ifdef X11
#ifndef lint
static char rcsid[] = "$Header$";
#endif lint


/***********************************************************
Copyright 1987, 1988 by Digital Equipment Corporation, Maynard, Massachusetts,
and the Massachusetts Institute of Technology, Cambridge, Massachusetts.

                        All Rights Reserved

Permission to use, copy, modify, and distribute this software and its 
documentation for any purpose and without fee is hereby granted, 
provided that the above copyright notice appear in all copies and that
both that copyright notice and this permission notice appear in 
supporting documentation, and that the names of Digital or MIT not be
used in advertising or publicity pertaining to distribution of the
software without specific, written prior permission.  

DIGITAL DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING
ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL
DIGITAL BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR
ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS
SOFTWARE.

******************************************************************/

/*
 * Canvas.c - Canvas widget
 *
 */

#define XtStrlen(s)		((s) ? strlen(s) : 0)

#include <X11/IntrinsicP.h>
#include <stdio.h>
#include <X11/Xos.h>
#include <ctype.h>
#include <X11/StringDefs.h>
#include "CanvasP.h"

/****************************************************************
 *
 * Full class record constant
 *
 ****************************************************************/

/* Private Data */

#define offset(field) XtOffset(CanvasWidget, field)
static XtResource resources[] = {
    {XtNforeground, XtCForeground, XtRPixel, sizeof(Pixel),
	offset(canvas.foreground), XtRString, "Black"},
    {XtNfont,  XtCFont, XtRFontStruct, sizeof(XFontStruct *),
	offset(canvas.font),XtRString, "Fixed"},
};

static void Initialize();
static void Realize();
static void Resize();
static void Redisplay();
static Boolean SetValues();
static void ClassInitialize();

CanvasClassRec canvasClassRec = {
  {
/* core_class fields */	
#define superclass		(&simpleClassRec)
    /* superclass	  	*/	(WidgetClass) superclass,
    /* class_name	  	*/	"Canvas",
    /* widget_size	  	*/	sizeof(CanvasRec),
    /* class_initialize   	*/	ClassInitialize,
    /* class_part_initialize	*/	NULL,
    /* class_inited       	*/	FALSE,
    /* initialize	  	*/	Initialize,
    /* initialize_hook		*/	NULL,
    /* realize		  	*/	Realize,
    /* actions		  	*/	NULL,
    /* num_actions	  	*/	0,
    /* resources	  	*/	resources,
    /* num_resources	  	*/	XtNumber(resources),
    /* xrm_class	  	*/	NULLQUARK,
    /* compress_motion	  	*/	TRUE,
    /* compress_exposure  	*/	TRUE,
    /* compress_enterleave	*/	TRUE,
    /* visible_interest	  	*/	FALSE,
    /* destroy		  	*/	NULL,
    /* resize		  	*/	Resize,
    /* expose		  	*/	Redisplay,
    /* set_values	  	*/	SetValues,
    /* set_values_hook		*/	NULL,
    /* set_values_almost	*/	XtInheritSetValuesAlmost,
    /* get_values_hook		*/	NULL,
    /* accept_focus	 	*/	NULL,
    /* version			*/	XtVersion,
    /* callback_private   	*/	NULL,
    /* tm_table		   	*/	NULL,
    /* query_geometry		*/	NULL,
  }
};
WidgetClass canvasWidgetClass = (WidgetClass)&canvasClassRec;
/****************************************************************
 *
 * Private Procedures
 *
 ****************************************************************/

static void ClassInitialize()
{

} /* ClassInitialize */


static void GetnormalGC(lw)
    CanvasWidget lw;
{
    XGCValues	values;

    values.foreground	= lw->canvas.foreground;
    values.font		= lw->canvas.font->fid;

    lw->canvas.normal_GC = XtGetGC(
	(Widget)lw,
	(unsigned) GCForeground | GCFont,
	&values);
}

/* ARGSUSED */
static void Initialize(request, new)
 Widget request, new;
{
    CanvasWidget lw = (CanvasWidget) new;

    GetnormalGC(lw);

} /* Initialize */


static void Realize(w, valueMask, attributes)
    register Widget w;
    Mask *valueMask;
    XSetWindowAttributes *attributes;
{

    (*superclass->core_class.realize) (w, valueMask, attributes);
    
} /* Realize */



/*
 * Repaint the widget window
 */

/* ARGSUSED */
static void Redisplay(w, event, region)
    Widget w;
    XEvent *event;
    Region region;
{
}


static void Resize(w)
    Widget w;
{
}

/*
 * Set specified arguments into widget
 */

/* ARGSUSED */
static Boolean SetValues(current, request, new)
    Widget current, request, new;
{
    CanvasWidget curlw = (CanvasWidget) current;
    CanvasWidget reqlw = (CanvasWidget) request;
    CanvasWidget newlw = (CanvasWidget) new;
    Boolean was_resized = False;

    /* we have to know if the size change is going to take
       before calling Resize() */
    if ((curlw->core.width != newlw->core.width ||
	 curlw->core.height != newlw->core.height) &&
	(XtMakeResizeRequest(current, newlw->core.width, newlw->core.height,
			     &newlw->core.width, &newlw->core.height)
	 == XtGeometryNo)) {
	newlw->core.width = curlw->core.width;
	newlw->core.height = curlw->core.height;
    }

    if (curlw->canvas.foreground != newlw->canvas.foreground
	|| curlw->canvas.font->fid != newlw->canvas.font->fid) {

	XtDestroyGC(curlw->canvas.normal_GC);
	GetnormalGC(newlw);
    }

    return( was_resized );
}
#endif
