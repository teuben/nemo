/* $Id$
 *
 * paneledit.c               version 1.0             February 2, 1989
 *
 * Allows dynamic repositioning of all panel items in all panels
 *
 * Authors:
 *                                      procedure move_item() by:
 *     Mark Phillips                    Chuck Musciano
 *     Department of Mathematics        Advanced Technology Department
 *     University of Maryland           Harris Corporation
 *     College Park, Maryland 20742     PO Box 37, MS 3A/1912
 *     (301) 454-2693                   Melbourne, FL 32902
 *     ARPA: mbp@lakisis.umd.edu        (407) 727-6131
 *                                      ARPA: chuck@trantor.harris-atd.com
 *
 * Everything except procedure "move_item":
 *      Copyright 1989 by Mark Phillips and the University of Maryland
 *
 * Procedure "move_item":
 *      Copyright 1988 by Chuck Musciano and Harris Corporation
 *
 * Permission to use, copy, modify, and  distribute this software and
 * its  documentation for  any purpose   and without   fee is  hereby
 * granted, provided that  the above copyright  notices appear in all
 * copies and that both  those  copyright notices and this permission
 * notice appear in supporting  documentation, and  that the names of
 * Mark Phillips, the  University  of Maryland,  Chuck Musciano,  and
 * Harris  Corporation  not be   used  in advertising   or  publicity
 * pertaining  to   distribution of   the software without  specific,
 * written prior   permission.   Mark Phillips,   the   University of
 * Maryland,  Chuck   Musciano,   and  Harris  Corporation  make   no
 * representations about the  suitability  of  this software  for any
 * purpose.   It is  provided  "as  is" without  express   or implied
 * warranty.
 */

#include <stdio.h>
#include <suntool/sunview.h>
#include <suntool/panel.h>

typedef	struct save_block_tag {
  caddr_t	client_data;
  caddr_t	event_proc;
} Save_block;

static Frame base_frame;
static Frame edit_frame;
static Panel edit_panel;
static Panel_item edit_cycle;
static int toggle_edit(), move_item();

#define		NO_DRAG		0
#define		DRAG_X		1
#define		DRAG_Y		2
static int dragging = NO_DRAG;

static int editing = 0;

/*-----------------------------------------------------------------------
 * Function:	paneledit_init
 * Description:	Set up editing of SunView panels
 * Args  IN:	frame: handle of frame
 * Returns:	success status: 0 for success, 1 for failure
 * Notes:	All panel items in all subwindows of frame can be edited.
 */
paneledit_init(frame)
Frame frame;
{
  int h;

  base_frame = frame;

  edit_frame =
    window_create(base_frame, FRAME,
		  FRAME_LABEL, "paneledit - 1.0",
		  FRAME_SHOW_LABEL, TRUE,
		  0);

  edit_panel = window_create(edit_frame, PANEL,0);
  edit_cycle =
    panel_create_item(edit_panel, PANEL_CYCLE,
		      PANEL_LABEL_STRING, "Panel Editing ",
		      PANEL_CHOICE_STRINGS, "Off", "On", 0,
		      PANEL_NOTIFY_PROC, toggle_edit,
		      0);
  window_fit(edit_panel);
  window_fit(edit_frame);
  window_set(edit_frame,
/*
	     WIN_X, (int)window_get(base_frame, WIN_WIDTH),
	     WIN_Y, 0,
*/
	     WIN_SHOW, TRUE,
	     0);
}

/*-----------------------------------------------------------------------
 * Function:	toggle_edit
 * Description:	toggle whether panel items are being edited
 * Notes:	this is the event procedure for the editing cycle item
 */
static int
toggle_edit()
{
  /* We toggle the editing state and apply new state to window tree */
  editing = !editing;
  set_editing_in_window(base_frame, editing);
}

/*-----------------------------------------------------------------------
 * Function:	set_editing_in_window
 * Description:	turn editing on or off in a specific window
 * Args  IN:	window: the window's handle
 *		edit: 0 to turn editing off, 1 to turn it on
 * Returns:	nothing
 * Notes:	This procedure calls itself to operate on the entire
 *		window tree with window as root
 */
static int
set_editing_in_window(window, edit)
     Window window;
     int edit;
{
  int winno;
  Window subwindow;

  winno = 0;
  do {
    /* Get next subwindow */
    subwindow = (Window)window_get(window, FRAME_NTH_WINDOW, winno);
    if (subwindow!=NULL) {
      /* Make sure this is not the edit frame! */
      if (subwindow != edit_frame) {
	switch ( (Window_type)(window_get(subwindow,WIN_TYPE)) ) {
	case FRAME_TYPE:
	  set_editing_in_window(subwindow, edit);
	  break;
	case PANEL_TYPE:
	  set_editing_in_panel(subwindow, edit);
	  break;
	default:
	  /* do nothing */
	  break;
	}
      }
    }
    ++winno;
  } while (subwindow != NULL);
}

/*-----------------------------------------------------------------------
 * Function:	set_editing_in_panel
 * Description:	turn editing on or off for all items in a panel
 * Args  IN:	panel: the panel's handle
 *		edit: 0 to turn editing off, 1 to turn it on
 * Returns:	nothing
 */
static int
set_editing_in_panel(panel, edit)
     Panel panel;
     int edit;
{
  Panel_item item;
  Save_block *b;

  /* Loop thru each item in this panel */
  panel_each_item(panel, item) {

    if (edit) {

      /* allocate a Save_block for item's current info */
      b = (Save_block*)malloc(sizeof(Save_block));
      if (b == NULL) {
	printf(stderr, "Sorry, no memory left for panel_edit !!\n");
	exit(1);
      }
      b->client_data = (caddr_t)panel_get(item, PANEL_CLIENT_DATA);
      b->event_proc = (caddr_t)panel_get(item, PANEL_EVENT_PROC);

      /* Store pointer to this Save_block as item's new client data */
      panel_set(item, PANEL_CLIENT_DATA, b, 0);

      /* For editing, use move_item as event proc */
      panel_set(item, PANEL_EVENT_PROC, move_item, 0);
    }

    else {

      /* restore previous event proc and client data */
      b = (Save_block*)panel_get(item, PANEL_CLIENT_DATA);
      panel_set(item, PANEL_EVENT_PROC, b->event_proc, 0);
      panel_set(item, PANEL_CLIENT_DATA, b->client_data, 0);

      /* free up the Save_block */
      free(b);
    }

  } panel_end_each;

}

/*-----------------------------------------------------------------------
 * Function:	move_item
 * Description:	A notify proc which allows items to move about and then
 *		  report their position.
 * Args  IN:	item: item's handle
 *		event: event which caused notification
 * Returns:	nothing
 * Notes: 	Copyright 1988 by Chuck Musciano and Harris Corporation
 */
static int
move_item(item, event)
     Panel_item	item;
     Event	*event;
{
  static int old_x, old_y;
  static Panel_item old_item;

  if (event_id(event) == MS_LEFT)
    if (event_is_down(event)) {
      dragging = DRAG_X;
      old_x = event_x(event);
      old_item = item;
    }
    else
      dragging = NO_DRAG;
  else if (event_id(event) == MS_MIDDLE)
    if (event_is_down(event)) {
      dragging = DRAG_Y;
      old_y = event_y(event);
      old_item = item;
    }
    else
      dragging = NO_DRAG;
  else if (event_id(event) == MS_RIGHT && event_is_down(event))
    printf("%5d %5d\n", panel_get(item, PANEL_ITEM_X),
	   panel_get(item, PANEL_ITEM_Y));
  else if (event_id(event) == LOC_DRAG)
    if (dragging == DRAG_X) {
      panel_set(old_item, PANEL_ITEM_X,
		panel_get(item, PANEL_ITEM_X) + event_x(event) - old_x, 0);
      old_x = event_x(event);
    }
    else if (dragging == DRAG_Y) {
      panel_set(old_item, PANEL_ITEM_Y,
		panel_get(item, PANEL_ITEM_Y) + event_y(event) - old_y, 0);
      old_y = event_y(event);
    }
    else
      panel_default_handle_event(item, event);
}
