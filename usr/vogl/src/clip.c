#include "vogl.h"

static int planes[] = {'\01', '\02', '\04', '\010', '\020', '\040'};

float wc[2][6];

/*
 * clip
 *
 * Clips a 3D line using Homogeneous clipping.
 * Reference: Newman and Sproull
 *
 */
void
clip(p0, p1)
	register Vector	p0, p1;
{
	float	t, t1, t2, dx, dy, dz, dw;
	int	vx, vy, c1, c2, i;

	c1 = MakeEdgeCoords(0, p0);
	c2 = MakeEdgeCoords(1, p1);

	if (!(c1 & c2)) {

		t1 = 0.0;
		t2 = 1.0;
		for (i = 0; i < 6; i++)
			if (wc[0][i] < 0.0 || wc[1][i] < 0.0) {
				t = wc[0][i] / (wc[0][i] - wc[1][i]);
				if (wc[0][i] < 0.0) {
					if (t > t1) 
					    t1 = t;
				} else if (t < t2) 
					    t2 = t;
			}
		 
		 if (t2 >= t1) {
			vdevice.cpVvalid = 1;
			dx = p1[V_X] - p0[V_X];
			dy = p1[V_Y] - p0[V_Y];
			dz = p1[V_Z] - p0[V_Z];
			dw = p1[V_W] - p0[V_W];
			if (t2 != 1.0) {
				p1[V_X] = p0[V_X] + t2 * dx;
				p1[V_Y] = p0[V_Y] + t2 * dy;
				p1[V_Z] = p0[V_Z] + t2 * dz;
				p1[V_W] = p0[V_W] + t2 * dw;
				vdevice.cpVvalid = 0;
			}
			if (t1 != 0.0) {
				p0[V_X] = p0[V_X] + t1 * dx;
				p0[V_Y] = p0[V_Y] + t1 * dy;
				p0[V_Z] = p0[V_Z] + t1 * dz;
				p0[V_W] = p0[V_W] + t1 * dw;
			}
			vdevice.cpVx = WtoVx(p0);
			vdevice.cpVy = WtoVy(p0);
			vx = WtoVx(p1);
			vy = WtoVy(p1);
			(*vdevice.dev.Vdraw)(vx, vy);
			vdevice.cpVx = vx;
			vdevice.cpVy = vy;
		}
	}
}

/*
 * MakeEdgeCoords
 *
 * calculates distance from point to each clipping plane in Homogeneous
 * clipping coordinates. Return code if on outside of any clipping plane
 *
 */
int
MakeEdgeCoords(i, p)
	int	i;
	Vector	p;
{
	int	j, k = 0;

	wc[i][0] = p[V_W] + p[V_X];
	wc[i][1] = p[V_W] - p[V_X];
	wc[i][2] = p[V_W] + p[V_Y];
	wc[i][3] = p[V_W] - p[V_Y];
	wc[i][4] = p[V_W] + p[V_Z];
	wc[i][5] = p[V_W] - p[V_Z];
	
	for (j = 0; j < 6; j++)
		if (wc[i][j] < 0.0)
			k |= planes[j];

	return(k);
}

/*
 * quickclip
 *
 * A variation on the above that assumes p0 is a valid position in device coords
 *
 */
void
quickclip(p0, p1)
	register Vector	p0, p1;
{
	register float	t, t1;
	register int	vx, vy, i;

	t1 = 1.0;

	wc[0][0] = p0[V_W] + p0[V_X];
	wc[0][1] = p0[V_W] - p0[V_X];
	wc[0][2] = p0[V_W] + p0[V_Y];
	wc[0][3] = p0[V_W] - p0[V_Y];
	wc[0][4] = p0[V_W] + p0[V_Z];
	wc[0][5] = p0[V_W] - p0[V_Z];
	wc[1][0] = p1[V_W] + p1[V_X];
	wc[1][1] = p1[V_W] - p1[V_X];
	wc[1][2] = p1[V_W] + p1[V_Y];
	wc[1][3] = p1[V_W] - p1[V_Y];
	wc[1][4] = p1[V_W] + p1[V_Z];
	wc[1][5] = p1[V_W] - p1[V_Z];

	for (i = 0; i < 6; i++)
		if (wc[0][i] >= 0.0 && wc[1][i] < 0.0) {
			t = wc[0][i] / (wc[0][i] - wc[1][i]);
			if (t < t1)
				    t1 = t;
		}
	 
	vdevice.cpVvalid = 1;
	if (t1 != 1.0) {
		p1[V_X] = p0[V_X] + t1 * (p1[V_X] - p0[V_X]);
		p1[V_Y] = p0[V_Y] + t1 * (p1[V_Y] - p0[V_Y]);
		p1[V_Z] = p0[V_Z] + t1 * (p1[V_Z] - p0[V_Z]);
		p1[V_W] = p0[V_W] + t1 * (p1[V_W] - p0[V_W]);
		vdevice.cpVvalid = 0;
	}

	vx = WtoVx(p1);
	vy = WtoVy(p1);
	(*vdevice.dev.Vdraw)(vx, vy);
	vdevice.cpVx = vx;
	vdevice.cpVy = vy;
}

