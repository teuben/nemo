#include "vogl.h"

void
pc_fill(n, x, y)
	int	n, *x, *y;
{
	register	int	i;
	int		tx[3], ty[3], xstart, xend;

	/*
	 * If it's a rectangle then things are even easier....
	 */
	if (n == 5 && x[0] == x[3] && x[1] == x[2] &&
		y[0] == y[1] && y[2] == y[3]) {
		int	low, high;
		xstart = x[0];
		xend = x[1];
		low = y[0];
		high = y[2];
		if (low >= high) {
			low = high;
			high = y[0];
		}
		for (i = low; i <= high; i++) {
			vdevice.cpVx = xstart;
			vdevice.cpVy = i;
			(*vdevice.dev.Vdraw)(xend, i);
		}
	} else {

		for (i = 0; i < n - 2; i++) {
			tx[0] = x[0];
			ty[0] = y[0];
			tx[1] = x[i+1];
			ty[1] = y[i+1];
			tx[2] = x[i+2];
			ty[2] = y[i+2];
			fill_tri(tx, ty);
		}
	}
}

#define ABS(x)		((x) > 0 ? (x) : -(x))
#define SIGN(x)		((x) > 0 ? 1 : ((x) == 0 ? 0 : (-1)))
#define	MAX(x,y)	((x) > (y) ? (x) : (y))

fill_tri(tx, ty)
	int tx[3], ty[3];
{
	register int y;
	register int	dx1, dx2, dy1, dy2, inc1, inc2, sx1, sx2, sy1,
	xl, xr, x1, x2, ix1, ix2, iy1, iy2, doit, doit2, y1, y2;

	/* 
	 * Order in y values
	 */
	y1 = 0;
	while (!y1) {
		y1 = 1;
		for (y2 = 0; y2 < 2; y2++)
			if (ty[y2] > ty[y2+1] || (ty[y2] == ty[y2+1] && tx[y2] < tx[y2+1])) {
				x1 = ty[y2];
				ty[y2] = ty[y2+1];
				ty[y2+1] = x1;
				x1 = tx[y2];
				tx[y2] = tx[y2+1];
				tx[y2+1] = x1;
				y1 = 0;
			}
	}
 	dx1 = tx[1] - tx[0];
	dx2 = tx[2] - tx[0];
	dy1 = ty[1] - ty[0];
	dy2 = ty[2] - ty[0];
	sx1 = SIGN(dx1);
	sx2 = SIGN(dx2);
	sy1 = SIGN(dy1);
	ix1 = ABS(dx1);
	ix2 = ABS(dx2);
	iy1 = ABS(dy1);
	iy2 = ABS(dy2);
	inc1 = MAX(ix1, iy1);
	inc2 = MAX(ix2, iy2);

	x1 = x2 = y1 = y2 = 0;
	xl = xr = tx[0];
	y = ty[0];

	while (y != ty[1]) {
		doit = 0;
		doit2 = 0;
		while (!doit) {		/* Wait until y changes */
			x1 += ix1;
			y1 += iy1;
			if (x1 > inc1) {
				x1 -= inc1;
				xl += sx1;
			}
			if (y1 > inc1) {
				y1 -= inc1;
				y += sy1;
				doit = 1;
			}
		}
		while (!doit2) {	/* Wait until y changes */
			x2 += ix2;
			y2 += iy2;
			if (x2 > inc2) {
				x2 -= inc2;
				xr += sx2;
			}
			if (y2 > inc2) {
				y2 -= inc2;
				doit2 = 1;
			}
		}
		
		vdevice.cpVx = xl;
		vdevice.cpVy = y;
		(*vdevice.dev.Vdraw)(xr, y);
	}
	dx1 = tx[2] - tx[1];
	dy1 = ty[2] - ty[1];
	sx1 = SIGN(dx1);
	sy1 = SIGN(dy1);
	ix1 = ABS(dx1);
	iy1 = ABS(dy1);
	inc1 = MAX(ix1, iy1);
	xl = tx[1];
	x1 = 0;
	while (y != ty[2]) {
		doit = 0;
		doit2 = 0;
		while (!doit) {		/* Wait until y changes */
			x1 += ix1;
			y1 += iy1;
			if (x1 > inc1) {
				x1 -= inc1;
				xl += sx1;
			}
			if (y1 > inc1) {
				y1 -= inc1;
				y += sy1; 
				doit = 1;
			}
		}
		while (!doit2) {	/* Wait until y changes */
			x2 += ix2;
			y2 += iy2;
			if (x2 > inc2) {
				x2 -= inc2;
				xr += sx2;
			}
			if (y2 > inc2) {
				y2 -= inc2;
				doit2 = 1;
			}
		}
		vdevice.cpVx = xl;
		vdevice.cpVy = y;
		(*vdevice.dev.Vdraw)(xr, y);
	}

	return(0);
}
