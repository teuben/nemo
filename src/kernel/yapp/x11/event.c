/*
 *   event.c - Part of my_x and yapp_x (NEMO) graphics driver programs
 *         Provides generic X11-based plotting functions
 *         
 *   Author:     Patrick (AdM) Frisch ( frisch@walhall.uni-sw.gwdg.de )
 *   Date:       sometimes in 1992
 *   Revised:    March 93
 *   Version:    1.01
 *   Module
 *   Remarks:    event handling (under construction !! Not finished yet ..)
 *
 *
 *   Global Functions:
 *
 *               x_ask_event(XEvent *)
 *               x_handle_key_event(XEvent *)
 *
 */

#include "config.h"
#include <X11/Xlib.h>
#include "event.h"

extern Display *x_dpy;
extern Window x_win;
extern Pixmap x_pixmap;
extern GC x_gc;
extern XWindowAttributes x_wa;

x_ask_event(xev)
XEvent *xev;
{
    long ev_mask = ButtonPressMask | StructureNotifyMask;
    XSelectInput(x_dpy, x_win, ev_mask);
    if(XCheckWindowEvent(x_dpy, x_win, ev_mask, xev))
    {
#ifdef DEBUG
	puts(event_names[xev->type]);
#endif
	return 1;
    }
    else
        return 0;
}

x_handle_key_event(xev)
XEvent xev;
{
    
}



