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

#ifndef PI
#define PI	3.14159261358979
#endif

#define	STEP	PI / 180.0

main()
{
	float a;
	int i;
	float v[2];
	short	val;

	vinit("mswin");
	winopen("bgnline/endline test");
	unqdevice(INPUTCHANGE);
	qdevice(KEYBD);
	/* 
	 * Wait for REDRAW event ...
	 */
	while (qread(&val) != REDRAW)
		;

	color(BLACK);
	clear();
	ortho2(-0.5, 2 * PI + .5, -1.5, 1.5);

	color(GREEN);
	bgnline();
		v[0] = -0.5;
		v[1] = 0.0;
		v2f(v);

		v[0] = 2 * PI + .5;
		v[1] = 0.0;
		v2f(v);
	endline();

	color(RED);
	bgnline();
		v[0] = 0.0;
		v[1] = -1.3;
		v2f(v);

		v[0] = 0.0;
		v[1] = 1.3;
		v2f(v);
	endline();


	color(YELLOW);
	bgnline(); 
		v[0] = v[1] = 0.0;
		v2f(v);

		for (a = 0.0; a <= 2 * PI; a += STEP) {
			v[0] = a;
			v[1] = sin(a);
			v2f(v);
		}
	endline();

	qread(&val);

	gexit(); 
}

