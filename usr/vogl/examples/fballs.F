c
c makesphere
c
c	make a sphere object
c
	subroutine makesp

	integer SPHERE
	real r, z, a, RADIUS, PI
	parameter (PI = 3.1415926535, RADIUS = 10.0, SPHERE = 1)

	call makeob(SPHERE)

c
c create the latitudinal rings
c
	    do 10 i = 0, 1800, 200
		call pushma
		    call rotate(i, 'y')
		    call circ(0.0, 0.0, RADIUS)
		call popmat
10	    continue
		
c
c create the longitudinal rings
c
	    call pushma
		call rotate(900, 'x')
		do 20 a = -900, 900, 200
		    r = RADIUS * cos(a * PI / 180.0)
		    z = RADIUS * sin(a * PI / 180.0)
		    call pushma
			call transl(0.0, 0.0, -z)
			call circ(0.0, 0.0, r)
		    call popmat
20		continue
	    call popmat

	call closeo

	end

c
c a demonstration of objects
c
	program fballs

#ifdef SGI
#include "fgl.h"
#include "fdevice.h"
#else
#include "fvogl.h"
#include "fvodevice.h"
#endif

	integer *2 val
	integer SPHERE
	real RADIUS
	parameter (RADIUS = 10.0)
	parameter(SPHERE = 1)

	call winope('fballs', 6)
	call unqdev(INPUTC)
	call qdevic(KEYBD)
c
c Read the initial REDRAW event
c
	idum = qread(val)
c
c set up our viewing transformation
c
	call perspe(900, 1.0, 0.001, 500.0)
	call lookat(13.0, 13.0, 8.0, 0.0, 0.0, 0.0, 0)

	call color(BLACK)
	call clear

c
c Call a routine to make the sphere object
c
	call makesp

c
c Now draw the sphere object scaled down. We use the pushmatrix
c and the popmatrix to preserve the transformation matrix so
c that only this sphere is drawn scaled. The callobj then enables
c us to draw the sphere we generated with makeobj in makesphere.
c
	call color(CYAN)

	call pushma
	    call scale(0.5, 0.5, 0.5)
	    call callob(SPHERE)
	call popmat

c
c now we draw the same sphere translated, with a different
c scale and color.
c
	call color(WHITE)

	call pushma
	    call transl(0.0, -1.4 * RADIUS, 1.4 * RADIUS)
	    call scale(0.3, 0.3, 0.3)
	    call callob(SPHERE)
	call popmat

c
c and maybe a few more times....
c

	call color(RED)

	call pushma
	    call transl(0.0, RADIUS, 0.7 * RADIUS)
	    call scale(0.2, 0.2, 0.2)
	    call callob(SPHERE)
	call popmat

	call color(GREEN)

	call pushma
	    call transl(0.0, 1.5 * RADIUS, -RADIUS)
	    call scale(0.15, 0.15, 0.15)
	    call callob(SPHERE)
	call popmat

	call color(YELLOW)

	call pushma
	    call transl(0.0, -RADIUS, -RADIUS)
	    call scale(0.12, 0.12, 0.12)
	    call callob(SPHERE)
	call popmat

	call color(BLUE)

	call pushma
	    call transl(0.0, -2.0*RADIUS, -RADIUS)
	    call scale(0.3, 0.3, 0.3)
	    call callob(SPHERE)
	call popmat

	call hfont('times.rb', 8)
	call ortho2(0.0, 1.0, 0.0, 1.0)
	call hcente(.true.)
	call htexts(0.08, 0.15)
	call move2(0.8, 0.5)
	call htexta(-90.0)
	call hchars('I''m very ordinary!', 18)

	idum = qread(val)

	call gexit

	end
