
c
c Shows various combinations of viewing and
c projection transformations.
c
	program fviews

#ifdef SGI
#include "fgl.h"
#include "fdevice.h"
#else
#include "fvogl.h"
#include "fvodevice.h"
#endif

	integer *2 val, vminx, vmaxx, vminy, vmaxy

	call winope('fviews', 6)
c
c we are interested in keyboard events
c
	call unqdev(INPUTC)
	call qdevic(KEYBD)
c
c Read the initial REDRAW event
c
	idum = qread(val)

	call color(BLACK)
	call clear
c
c open a hershey font
c
	call hfont('times.r', 7)

c
c we want to draw just within the boundaries of the screen
c
        call getvie(vminx, vmaxx, vminy, vmaxy)
	maxx = vmaxx
	minx = vminx
	maxy = vmaxy
	miny = vminy
        call viewpo(maxx / 10, maxx / 10 * 9, maxy / 10, maxy / 10 * 9)
 
c
c set the world size
c
	call ortho2(-5.0, 5.0, -5.0, 5.0)

c
c draw a boundary frame
c
	call color(RED)
	call rect(-5.0, -5.0, 5.0, 5.0)

c
c set up a perspective projection with a field of view of
c 40.0 degrees, aspect ratio of 1.0, near clipping plane 0.1,
c and the far clipping plane at 1000.0.
c
	call perspe(400, 1.0, 0.1, 1000.0)
c
c we want the drawing to be done with our eye point at (5.0, 8.0, 5.0)
c looking towards (0.0, 0.0, 0.0). The last parameter gives a twist
c in degrees around the line of sight, in this case zero.
c
	call lookat(5.0, 8.0, 5.0, 0.0, 0.0, 0.0, 0)

	call drawte

c
c set the text size
c
	call htexts(0.6, 0.9)

	call move2(-4.5, -4.5)
	call hchars('perspective/lookat', 18)

	idum = qread(val)

c
c window can also be used to give a perspective projection. Its
c arguments are 6 clipping planes, left, right, bottom, top, near,
c and far.
c
	call window(-5.0, 5.0, -5.0, 5.0, -5.0, 5.0)
c
c as window replaces the current transformation matrix we must
c specify our viewpoint again.
c
	call lookat(5.0, 8.0, 5.0, 0.0, 0.0, 0.0, 0)

	call color(BLACK)
	call clear

	call color(GREEN)
	call rect(-5.0, -5.0, 5.0, 5.0)

	call drawte

	call htexts(0.6, 0.9)
	call move2(-4.5,-4.5)
	call hchars('window/lookat', 13)

	idum = qread(val)

c
c set up our original perspective projection again.
c
	call perspe(400, 1.0, 0.1, 1000.0)
c
c polarview also specifies our viewpoint, but, unlike lookat, in polar
c coordinates. Its arguments are the distance from the world origin, an
c azimuthal angle in the x-y plane measured from the y axis, an 
c incidence angle in the y-z plane measured from the z axis, and a
c twist around the line of sight.
c
	call polarv(15.0, 300, 300, 300)

	call color(BLACK)
	call clear

	call color(MAGENT)
	call rect(-5.0, -5.0, 5.0, 5.0)

	call drawte

	call move2(-4.5,-4.5)
	call htexts(0.6, 0.9)
	call hchars('perspective/polarview', 21)

	idum = qread(val)

c
c once more with window for comparison
c
	call window(-4.0, 4.0, -4.0, 4.0, -4.0, 4.0)
	call polarv(6.0, 2000, -3000, 7000)

	call color(BLACK)
	call clear

	call color(YELLOW)
	call rect(-5.0, -5.0, 5.0, 5.0)

	call drawte

	call move2(-4.5,-4.5)
	call htexts(0.6, 0.9)
	call hchars('window/polarview', 16)

	idum = qread(val)

	call gexit
	
	end

c
c drawtetra 
c 
c generate a tetrahedron as a series of move draws 
c
c
	subroutine drawte
	
	integer WHITE
	parameter (WHITE = 7)

	call move(-0.5,  0.866, -0.5)
	call draw(-0.5, -0.866, -0.5)
	call draw( 1.0,  0.0,   -0.5)
	call draw(-0.5,  0.866, -0.5)
	call draw( 0.0,  0.0,    1.5)
	call draw(-0.5, -0.866, -0.5)
	call move( 1.0,  0.0,   -0.5)
	call draw( 0.0,  0.0,    1.5)
	
c 
c    Label the vertices.
c
	call color(WHITE)
	call htexts(0.3, 0.5)	
	call move(-0.5,  0.866, -0.5)
	call hdrawc('a')
	call move(-0.5, -0.866, -0.5)
	call hdrawc('b')
	call move( 1.0,  0.0,   -0.5)
	call hdrawc('c')
	call move( 0.0,  0.0,    1.5)
	call hdrawc('d')

	end
