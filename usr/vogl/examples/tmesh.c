/*
 * Demonstrate triangular mesh
 */

#ifdef SGI
#include "gl.h" 
#include "device.h" 
#else
#include "vogl.h"
#include "vodevice.h"
#endif

#include <math.h>

#define NTRIANGLE 20

Coord	cs[NTRIANGLE], sn[NTRIANGLE];

main()
{
	int	i, itest,  dobackface, dofill, dodouble;
	char	ans;
	float	H;
	short	idata;
	int	xr, yr;

	printf("Backfacing ON or OFF (Y/N)? ");
	ans = getchar();
	(void)getchar();
	dobackface = (ans == 'y' || ans == 'Y');

	printf("Fill the polygons (Y/N)? ");
	ans = getchar();
	(void)getchar();
	dofill = (ans == 'y' || ans == 'Y');

 	printf("double buffer (Y/N)? ");
	ans = getchar();
	dodouble = (ans == 'y' || ans == 'Y');

	winopen("piston");
	/* 
	 * Wait for REDRAW event ...
	 */
	while (qread(&idata) != REDRAW)
		;

	if (dodouble)
		doublebuffer();

	gconfig();

	unqdevice(INPUTCHANGE);
	qdevice(QKEY);
	qdevice(ESCKEY);
        qdevice(REDRAW);

	makecyl();

	polymode(PYM_LINE);
	if (dofill)
		polymode(PYM_FILL);

	if (dobackface)
		backface(1);
/*
 * set up a perspective projection with a field of view of
 * 40.0 degrees, aspect ratio of 1.0, near clipping plane 0.1,
 * and the far clipping plane at 1000.0.
 */
	perspective(400, 1.5, 0.1, 600.0);
	lookat(0.0, -6.0, 4., 0.0, 0.0, 0.0, 0);

/*
 * here we loop back here adnaseum until someone hits a key
 */
	xr = yr = 0;

 	while(1) {
		for (i = 0; i < 360; i += 5) {
			color(BLACK);
			clear();
			color(RED);
			H = 1.0 + cos(2.0 * 3.14159265*i / 180.0);

			yr = 500 - (int)getvaluator(MOUSEY);
			xr = 500 - (int)getvaluator(MOUSEX);
			yr = 500 - (int)getvaluator(MOUSEY);
			xr *= 3;
			yr *= 3;

			pushmatrix();
				rotate(xr, 'x');
				rotate(yr, 'y');
				piston(H);
			popmatrix();

			if (dodouble)
				swapbuffers();

			if (qtest()) {
				itest = qread(&idata);
				if( itest == QKEY || itest == ESCKEY) {
					 gexit();
					 exit(0);
				}
			}
		}
	}
}

makecyl()
{
	float	dphi, phi, pi = 3.141592653589;
	int	k;
	

	dphi = 2.*pi/(NTRIANGLE - 1);

	for (k = 0; k < NTRIANGLE; k++) {
		phi = k * dphi;
		cs[k] = cos(phi);
		sn[k] = sin(phi);
	}
}

piston(H)
	float	H;
{
	Coord	vec[3];
	int	j, k;
	float	HH;
/*
 * do the sides
 */
	bgnqstrip();
	for (k = 0; k < NTRIANGLE; k++) {
            vec[0] = cs[k];
            vec[1] = sn[k];
            vec[2] = H;
            v3f(vec);
            vec[2] = 0;
            v3f(vec);
	}
	endqstrip();
/*
 * do the ends
 */
	color(CYAN);
	for (k = -1; k <= 1; k += 2) {
		HH = H * (k + 1) / 2;
		bgntmesh();
		vec[0] = 0;
		vec[1] = 0;
		vec[2] = HH;
		v3f(vec);
		for (j = 0; j < NTRIANGLE; j++) {
			vec[0] = cs[j];
			vec[1] = k * sn[j];
			v3f(vec);
			swaptmesh();
		}
		endtmesh();
	}
}
