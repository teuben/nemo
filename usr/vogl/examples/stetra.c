/*
 * Demonstrate a rotating translating tetrahedron.
 * Modified by garym@virtual.rose.utoronto.ca (Gary Lawrence Murphy)
 * to be "stereo".
 */

#include <stdio.h>

#ifdef SGI
#include "gl.h"
#include "device.h"
#else
#include "vogl.h"
#include "vodevice.h"
#endif

#ifndef TC
#include <math.h>
#else
extern double sin(), cos();
#endif

#define	TETRAHEDRON	1L
#define	NSIDES	3
#define	NFACES	4
#define	NPNTS	4

Coord	points[NPNTS][3] = {
/* 	{-0.5, 0.866, -0.667}, */
/* 	{-0.5, -0.866, -0.667}, */
/* 	{ 1.0, 0.0, -0.667}, */
/* 	{ 0.0, 0.0, 1.334} */

  { -1, 1, -1 },
  { -1, -1, 1 },
  {  1,  1, 1 },
  {  1, -1, -1}

};

int	faces[NFACES][NSIDES] = {
	{2, 1, 0},
	{0, 1, 3},
	{1, 2, 3},
	{2, 0, 3}
};

int	colface[NFACES] = {
		GREEN,
		YELLOW,
		CYAN,
		MAGENTA
};

#define ZEYE 9.0
#define XEYE 0.1
#define LEFTEYE -1.0
#define RIGHTEYE 1.0

#define SETCAMERA( Eye )  	translate((Eye)*XEYE, 0.0, 0.0 )

int main(int argc, char	**argv)
{
	char	dev[20];
	int	i;
	int	rotval = 0, drotval = 2, irot = 0, srot;
	float	R = 3.0, tx = 0.0, tz = R;
	int	do_backface = 0;
	int	do_fill = 0;
	short	val;

	for (i = 1; i < argc; i++) {
	  if (!strcmp(argv[i], "-b"))
		do_backface = 1;
	  if (!strcmp(argv[i], "-f"))
		do_fill = 1;
	}

	prefsize(400L, 400L);

	winopen("tetra");          /* set up device */

	qdevice(ESCKEY);
	qdevice(QKEY);
	unqdevice(INPUTCHANGE);

	doublebuffer();
	gconfig();


	polymode(PYM_LINE);
	if (do_fill)
	  polymode(PYM_FILL);

	if (do_backface)
		backface(1);
	
	/*
	 * set up a perspective projection with a field of view of
	 * 36.0 degrees, aspect ratio of 1.0, near clipping plane 0.1,
	 * and the far clipping plane at 1000.0.
	 */
	
	perspective(360, 1.0, 0.001, 15.0);
	lookat(0.0, 0.0, ZEYE, 0.0, 0.0, 0.0, 0);

	/*
	 * Make a tetrahedron object
	 */

	maketetra();

/* 	scale(7, 7, 7); */
/* 	callobj(TETRAHEDRON); */
/* 	rotate(900, 'z'); */
/* 	callobj(TETRAHEDRON); */

	do {
	  for (rotval = 0; rotval < 3600; rotval += drotval) {
		color(BLACK);
		clear();

		for (i=0; i < 2; i ++) {

		  SETCAMERA( (i) ? RIGHTEYE : LEFTEYE );
		  color( (i)? RED : BLUE );
			  
		  irot = (irot + 3) % 3600;
		  srot = rotval*sin(irot/1800.0);
		  if (srot<0) srot += 3600;

		  pushmatrix();
		  {
			/*
			 * Rotate the whole scene...
			 */
			rotate(rotval, 'x');
			rotate(rotval, 'z');
			
			pushmatrix();
			{
			  rotate(900, 'x');
			  circ(0.0, 0.0, R);
			}
			popmatrix();
			
			move(0.0, 0.0, 0.0);
			
			pushmatrix();
			{
			  rotate(450,'z');
			  rotate(450,'x');
			  scale(0.4, 
					0.4+sin(rotval/1800.0)/20.0, 
					0.4+cos(irot/1800.0)/20.0);
			  
			  callobj(TETRAHEDRON);
			  rotate(900, 'x');
			  callobj(TETRAHEDRON);
			}
			popmatrix();
			
			/*
			 * Remember! The order of the transformations is
			 * the reverse of what is specified here in between
			 * the pushmatrix and the popmatrix. These ones don't
			 * accumulate because of the push and pop.
			 */
			pushmatrix();
			{
			  translate(tx+sin(irot/1800.0), 0.0, tz+cos(irot/1800.0));
			  circ(0.0,0.0,0.1);
			  rotate(450,'z');
			  rotate(450,'x');
			  circ(0.0,0.0,0.1);
			  rotate(irot, 'y');
			  scale(0.2, 0.2, 0.2);
			  callobj(TETRAHEDRON);
			  rotate(900, 'y');
			  callobj(TETRAHEDRON);
			}
			popmatrix();
		  }
		  popmatrix();
		}
		tz = R * cos((double)(rotval * 3.1415926535 / 180));
		tx = R * sin((double)(rotval * 3.1415926535 / 180));
		
		swapbuffers();
		
		if (qtest()) {
		  gexit();
		  exit(0);
		}
	  }
	  
	} while (1);
  }

/*
 * maketetra
 *
 *	draw a tetrahedron
 */
maketetra()
{
	int	i, j;

	makeobj(TETRAHEDRON);

	for (i = 0; i < NFACES; i++) {
/*		color(colface[i]); */
		bgnpolygon();
			for (j = 0; j < NSIDES; j++) 
				v3f(points[faces[i][j]]);
		endpolygon();
	}

	closeobj();
}
