#include <stdio.h>

#ifdef SGI
#include "gl.h"
#include "device.h"
#else
#include "vogl.h"
#include "vodevice.h"
#endif

#ifndef PI
#define PI	3.1415926535
#endif

#ifndef TC
#include <math.h>
#else
extern double	sin(), cos();
#endif

#define RADIUS	10.0
#define	SPHERE	1L

void	showroundtext(), makesphere();

/*
 * most of the things in this program have been done before but it has
 * a certain novelty value.
 */
main()
{
	int	i;
	short	val;
	float	r, z, a;

	winopen("world");
	qdevice(KEYBD);
	unqdevice(INPUTCHANGE);
	/* 
	 * Wait for REDRAW event ...
	 */
	while (qread(&val) != REDRAW)
		;

	hfont("futura.m");

	perspective(800, 1.0, 0.001, 50.0);
	lookat(13.0, 13.0, 8.0, 0.0, 0.0, 0.0, 0);

	color(BLACK);
	clear();

	makesphere();

	/*
	 * draw the main one in cyan
	 */
	color(CYAN);

	callobj(SPHERE);

	/*
	 * draw a smaller one outside the main one in white
	 */
	color(WHITE);

	pushmatrix();
		translate(0.0, -1.4 * RADIUS, 1.4 * RADIUS);
		scale(0.3, 0.3, 0.3);
		callobj(SPHERE);
	popmatrix();

	/*
	 * scale the text
	 */
	hboxfit(2.0 * PI * RADIUS, 0.25 * RADIUS, 31);

	/*
	 * now write the text in rings around the main sphere
	 */

	color(GREEN);
	showroundtext("Around the world in eighty days ");

	color(BLUE);
	/*
	 * note: that software text is rotated here as
	 * anything else would be whether you use textang
	 * or rotate depends on what you are trying to do.
	 * Experience is the best teacher here.
	 */
	rotate(900, 'x');
	showroundtext("Around the world in eighty days ");

	color(RED);
	rotate(900, 'z');
	showroundtext("Around the world in eighty days ");

	qread(&val);

	gexit();
}

/*
 * showroundtext
 *
 *	draw string str wrapped around a circle in 3d
 */
void
showroundtext(str)
	char	*str;
{
	int	i, inc;

	inc = 3600 / strlen(str);

	for (i = 0; i < 3600; i += inc) {
		pushmatrix();
			/*
			 * find the spot on the edge of the sphere
			 * by making it (0, 0, 0) in world coordinates
			 */
			rotate(i, 'y');
			translate(0.0, 0.0, RADIUS);

			move(0.0, 0.0, 0.0);

			hdrawchar(*str++);
		popmatrix();
	}
}

/*
 * makesphere
 *
 *	create the sphere object
 */
void
makesphere()
{
	float	i, r, z, a;

	makeobj(SPHERE);

		for (i = 0; i < 180; i += 20) {
			pushmatrix();
				rotate((int)i * 10, 'y');
				circ(0.0, 0.0, RADIUS);
			popmatrix();
		}
		
		pushmatrix();
			rotate(900, 'x');
			for (a = -90.0; a < 90.0; a += 20.0) {
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
