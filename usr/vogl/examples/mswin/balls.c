#include <stdio.h>

#ifdef SGI
#include "gl.h"
#include "device.h"
#else
#include "vogl.h"
#include "vodevice.h"
#endif

#ifndef PI
#define PI 3.1415926535
#endif

#ifndef TC
#include <math.h>
#else
extern double	sin(), cos();
#endif

#define RADIUS 10.0
#define	SPHERE	1L

/*
 * makesphere
 *
 *	make a sphere object
 */
void
makesphere()
{
	float	r, z;
	int	i, a;

	makeobj(SPHERE);

		/*
		 * create the latitudinal rings
		 */
		for (i = 0; i < 1800; i += 200) {
			pushmatrix();
				rotate(i, 'y');
				circ(0.0, 0.0, RADIUS);
			popmatrix();
		}
		
		/*
		 * create the longitudinal rings
		 */
		pushmatrix();
			rotate(900, 'x');
			for (a = -900; a < 900; a += 200) {
				r = RADIUS * cos((double)a * PI / 180.0);
				z = RADIUS * sin((double)a * PI / 180.0);
				pushmatrix();
					translate(0.0, 0.0, -z);
					circ(0.0, 0.0, r);
				popmatrix();	
			}
		popmatrix();

	closeobj();
}

/*
 * a demonstration of objects
 */
main()
{
	short	val;

	vinit("mswin");
	winopen("balls");

	unqdevice(INPUTCHANGE);
	qdevice(KEYBD);
	/* 
	 * Wait for REDRAW event ...
	 */
	while (qread(&val) != REDRAW)
		;

	/*
	 * set up our viewing transformation
	 */
	perspective(900, 1.0, 0.001, 500.0);
	lookat(13.0, 13.0, 8.0, 0.0, 0.0, 0.0, 0);

	color(BLACK);
	clear();

	/*
	 * Call a routine to make the sphere object
	 */
	makesphere();

	/*
	 * Now draw the sphere object scaled down. We use the pushmatrix
	 * and the popmatrix to preserve the transformation matrix so
	 * that only this sphere is drawn scaled.
	 */
	color(CYAN);

	pushmatrix();
		scale(0.5, 0.5, 0.5);
		callobj(SPHERE);
	popmatrix();

	/*
	 * now we draw the same sphere translated, with a different
	 * scale and color.
	 */

	color(WHITE);

	pushmatrix();
		translate(0.0, -1.4 * RADIUS, 1.4 * RADIUS);
		scale(0.3, 0.3, 0.3);
		callobj(SPHERE);
	popmatrix();

	/*
	 * and maybe a few more times....
	 */


	color(RED);

	pushmatrix();
		translate(0.0, RADIUS, 0.7 * RADIUS);
		scale(0.2, 0.2, 0.2);
		callobj(SPHERE);
	popmatrix();

	color(GREEN);

	pushmatrix();
		translate(0.0, 1.5 * RADIUS, -RADIUS);
		scale(0.15, 0.15, 0.15);
		callobj(SPHERE);
	popmatrix();

	color(YELLOW);

	pushmatrix();
		translate(0.0, -RADIUS, -RADIUS);
		scale(0.12, 0.12, 0.12);
		callobj(SPHERE);
	popmatrix();

	color(BLUE);

	pushmatrix();
		translate(0.0, -2.0*RADIUS, -RADIUS);
		scale(0.3, 0.3, 0.3);
		callobj(SPHERE);
	popmatrix();

	hfont("times.rb");
	ortho2(0.0, 1.0, 0.0, 1.0);
	hcentertext(1);
	htextsize(0.08, 0.15);
	move2(0.8, 0.5);
	htextang(-90.0);
	hcharstr("I'm very ordinary!");

	qread(&val);

	gexit();
}
