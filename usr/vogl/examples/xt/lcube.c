
#include <stdio.h>
#include "vogl.h"

#define	CUBE_SIZE	200.0
#define	TRANS		25.0
#define	SCAL		0.1
#define FACE		1
#define FILLED		2
#define OUTLINE		3

extern Object	filledthing, outlinething;
extern float	tdir;
extern float	scal;
extern int	but, nplanes;
int	i, n;

setup_lcube()
{
	window(-800.0, 800.0, -800.0, 800.0, -800.0, 800.0);
	lookat(0.0, 0.0, 1500.0, 0.0, 0.0, 0.0, 0);

	linewidth(3);

	/*
	 * Start with a very ordinary filled cube like the old demo..
	 */


	makecube(filledthing);
	makecube(outlinething);

	backbuffer(1);

}

makecube(obj)
	int	obj;
{

	void (*fun)();

	makeobj(obj);
		if (obj == outlinething)
			polymode(PYM_LINE);
		else
			polymode(PYM_FILL);

		pushmatrix();
			translate(0.0, 0.0, CUBE_SIZE);
			color(RED);
			rectf(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
		popmatrix();

		pushmatrix();
			translate(CUBE_SIZE, 0.0, 0.0);
			rotate(900, 'y');
			color(GREEN);
			rectf(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
		popmatrix();

		pushmatrix();
			translate(0.0, 0.0, -CUBE_SIZE);
			rotate(1800, 'y');
			color(BLUE);
			rectf(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
		popmatrix();

		pushmatrix();
			translate(-CUBE_SIZE, 0.0, 0.0);
			rotate(-900, 'y');
			color(CYAN);
			rectf(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
		popmatrix();

		pushmatrix();
			translate(0.0, CUBE_SIZE, 0.0);
			rotate(-900, 'x');
			color(MAGENTA);
			rectf(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
		popmatrix();

		pushmatrix();
			translate(0.0, -CUBE_SIZE, 0.0);
			rotate(900, 'x');
			color(YELLOW);
			rectf(-CUBE_SIZE, -CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
		popmatrix();

	closeobj();
}
