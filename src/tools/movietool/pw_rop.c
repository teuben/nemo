#include "movietool.h"
#include <stdio.h>
#include <suntool/canvas.h>
#include <suntool/scrollbar.h>

/* global declarations */

extern Canvas	canvas;
extern Pixwin	*canvaspw;
extern Image	*current_frame;
extern int	encoded_flag;
extern int	playing;
extern int	clipping;

void my_pw_rop (pw, dx, dy, dw, dh, op, pr, sx, sy, zoom)
Pixwin *pw;
int dx, dy, dw, dh, op, sx, sy, zoom;
struct pixrect *pr;
{
	struct pixrect *pr_tmp;
	register int y, z;

	if (encoded_flag && current_frame->ras_type == RT_BYTE_ENCODED) {
		if (! clipping && playing &&
			pr->pr_depth == canvaspw->pw_pixrect->pr_depth) {

			/* This approach isn't clean since it scribbles directly
		 	 * onto the display without any attempt at clipping.
		 	 * However, the speed it offers is attractive.
		 	 * If you zoom the tool to full-screen size,
		 	 * you (hopefully) won't notice the problems. */

			/* Add pixwin coordinates relative to screen origin */
			dx += canvaspw->pw_clipdata->pwcd_screen_x;
			dy += canvaspw->pw_clipdata->pwcd_screen_y;
			/* Add scrollbar thickness */
			dx += (int) scrollbar_get (window_get
				(canvas, WIN_VERTICAL_SCROLLBAR),
				SCROLL_THICKNESS);
			dy += (int) scrollbar_get (window_get
				(canvas, WIN_HORIZONTAL_SCROLLBAR),
				SCROLL_THICKNESS);
			if (zoom) {
				if (pr_decode_and_zoom (canvaspw->pw_pixrect,
					dx, dy, pr, zoom) == PIX_ERR)
					fprintf(stderr,
					"Error in pr_decode_and_zoom\n");
			} else {
				if (pr_decode (canvaspw->pw_pixrect,
					dx, dy, pr) == PIX_ERR)
					fprintf(stderr, "Error in pr_decode\n");
			}
		} else {
			/* Copy to intermediate pixrect */
			if (zoom) {
				pr_tmp = mem_create (pr->pr_size.x * zoom,
					pr->pr_size.y * zoom, pr->pr_depth);
				(void) pr_decode_and_zoom (pr_tmp, 0, 0, pr,
					zoom);
				dw *= zoom;
				dh *= zoom;
			} else {
				pr_tmp = mem_create(pr->pr_size.x,
					pr->pr_size.y, pr->pr_depth);
				(void) pr_decode (pr_tmp, 0, 0, pr);
			}
			(void) pw_rop (canvaspw, dx, dy, dw, dh,
				op, pr_tmp, sx, sy);
			pr_close (pr_tmp);
		}
	} else
		if (zoom) {
			for (y=0; y<pr->pr_size.y; y++) {
				for (z=0; z<zoom; z++) {
					(void)pw_rop (canvaspw, dx, z+zoom*y+dy,
						dw, 1, op, pr, sx, y+sy);
				}
			}
		} else
			(void)pw_rop (canvaspw, dx, dy, dw, dh, op, pr, sx, sy);
}
