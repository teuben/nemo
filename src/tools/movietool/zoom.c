#include <suntool/sunview.h>

#define ROP PIX_SRC | PIX_DONTCLIP
#define PIXRECT_NULL (struct pixrect *)NULL

/* pr_mag copies a magnified view of spr to dpr using pixel replication  */

struct pixrect *pr_mag (spr, magx, magy)
struct pixrect *spr;	/* source pixrect */
int magx, magy;	 	/* x and y magnifications */
{
	register short xmax = spr->pr_size.x, ymax = spr->pr_size.y - 1;
	register short x, y, i, delta;
	register struct pixrect *dpr;

	if (magx < 1 || magx > 10 || magy < 1 || magy > 10)
		return PIXRECT_NULL;

	if ((dpr = mem_create (xmax*magx, (ymax+1)*magy, spr->pr_depth)) ==
		PIXRECT_NULL)
		return PIXRECT_NULL;

	if (magx == 1)
		pr_rop (dpr, 0, 0, xmax, ymax+1, ROP, spr, 0, 0);
	else {
		for (x = 0; x < xmax; x++) {
			i = x*magx;
			for (delta = 0; delta < magx; delta++, i++)
				pr_rop (dpr, i, 0, 1, ymax, ROP, spr, x, 0);
		}
		xmax *= magx;
	}

	if (magy > 1) {
		for (y = ymax; y >= 0; y--) {
			i = y*magy;
			for (delta = 0; delta < magy; delta++, i++)
				pr_rop (dpr, 0, i, xmax, 1, ROP, dpr, 0, y);
		}
	}

	return dpr;
}
