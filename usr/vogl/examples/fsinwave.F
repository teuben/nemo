	program sine

#ifdef SGI
#include "fgl.h"
#include "fdevice.h"
#else
#include "fvogl.h"
#include "fvodevice.h"
#endif

	parameter(PI = 3.141592653589)
	integer *2 val
	real v(0:1)

	step = PI / 180.0

	call winope('bgnline/endline test', 20)
	call unqdev(INPUT)
	call qdevic(KEYBD)
c
c Read the initial REDRAW event
c
	idum = qread(val)

	call color(BLACK)
	call clear
	call ortho2(-0.5, 2 * PI + .5, -1.5, 1.5)

	call color(GREEN)
	call bgnlin
		v(0) = -0.5
		v(1) = 0.0
		call v2f(v)

		v(0) = 2 * PI + .5
		v(1) = 0.0
		call v2f(v)
	call endlin

	call color(RED)
	call bgnlin
		v(0) = 0.0
		v(1) = -1.3
		call v2f(v)

		v(0) = 0.0
		v(1) = 1.3
		call v2f(v)
	call endlin


	call color(YELLOW)
	call bgnlin
		v(0) = 0.0
		v(1) = 0.0
		call v2f(v)

		do 10 a = 0.0, 2 * PI, step
			v(0) = a
			v(1) = sin(a)
			call v2f(v)
10		continue
	call endlin

	idum = qread(val)

	call gexit 

	end
