c
c Demonstrate just how much you can put in an object
c
	program fobjviews
	
#ifdef SGI
#include "fgl.h"
#include "fdevice.h"
#else
#include "fvogl.h"
#include "fvodevice.h"
#endif

	integer *2 val, vminx, vmaxx, vminy, vmaxy

	integer CUBE, TOPLEFT, TOPRIGHT, BOTTOMLEFT, BOTTOMRIGHT
	parameter(CUBE = 1, TOPLEFT = 2, TOPRIGHT = 3,
     +            BOTTOMLEFT = 4, BOTTOMRIGHT = 5)

	call winope('fobjvws', 7)

	call getvie(vminx, vmaxx, vminy, vmaxy)

	minx = vminx
	maxx = vmaxx
	miny = vminy
	maxy = vmaxy

	call unqdev(INPUTC)
	call qdevic(SPACEK)
c
c Read the initial REDRAW event
c
	idum = qread(val)
	call pushvi

	call hfont('futura.m', 8)
	call htexts(0.5, 0.9)

	call color(BLACK)
	call clear

c
c make an object that represents the cube
c
	call makecu

c
c set up an object which draws in the top left of the screen.
c
	call makeob(TOPLEFT)
		call viewpo(minx, (maxx - minx) / 2,
     +                     (maxy - miny) / 2, maxy)

		call ortho2(-5.0, 5.0, -5.0, 5.0)

		call color(RED)

		call rect(-5.0, -5.0, 5.0, 5.0)

		call perspe(400, 1.0, 0.1, 1000.0)
		call lookat(5.0, 8.0, 5.0, 0.0, 0.0, 0.0, 0)
c
c          Call an object within another object
c
		call callob(CUBE)

		call color(GREEN)

		call move2(-4.5, -4.5)
		call hchars('perspective/lookat', 18)
	call closeo

c
c now set up one which draws in the top right of the screen
c
	call makeob(TOPRIGHT)
		call viewpo((maxx - minx) / 2, maxx,
     +                      (maxy - miny) / 2, maxy)
		call ortho2(-5.0, 5.0, -5.0, 5.0)

		call color(GREEN)

		call rect(-5.0, -5.0, 5.0, 5.0)

		call window(-5.0, 5.0, -5.0, 5.0, -5.0, 5.0)
		call lookat(5.0, 8.0, 5.0, 0.0, 0.0, 0.0, 0)

		call callob(CUBE)

		call color(RED)

		call move2(-4.5, -4.5)
		call hchars('window/lookat', 13)
	call closeo

c
c try the bottom left
c
	call makeob(BOTTOMLEFT)
		call viewpo(minx, (maxx - minx) / 2,
     +                      miny, (maxy - miny) / 2)

		call ortho2(-5.0, 5.0, -5.0, 5.0)

		call color(MAGENT)

		call rect(-5.0, -5.0, 5.0, 5.0)

		call perspe(400, 1.0, 0.1, 1000.0)
		call polarv(15.0, 300, 300, 300)

		call callob(CUBE)

		call color(YELLOW)

		call move2(-4.5, -4.5)
		call hchars('perspective/polarview', 21)
	call closeo

c
c and the bottom right
c
	call makeob(BOTTOMRIGHT)
		call viewpo((maxx - minx) / 2, maxx,
     +                      miny, (maxy - miny) / 2)
		call ortho2(-5.0, 5.0, -5.0, 5.0)

		call color(CYAN)

		call rect(-5.0, -5.0, 5.0, 5.0)

		call window(-5.0, 5.0, -5.0, 5.0, -5.0, 5.0)
		call polarv(8.0, -180, -30, 180)

		call callob(CUBE)

		call color(BLUE)

		call move2(-4.5, -4.5)
		call hchars('window/polarview', 16)
	call closeo

c
c now draw them
c
	call callob(TOPLEFT)
	call callob(TOPRIGHT)
	call callob(BOTTOMLEFT)
	call callob(BOTTOMRIGHT)

	idum = qread(val)

	call gexit

	end

c
c makecube
c
c set up a cube
c
	subroutine makecu
	integer CUBE
	parameter (CUBE = 1)

	call makeob(CUBE)

c
c The border around the cube
c
	    call rect(-5.0, -5.0, 10.0, 10.0)

c
c Make the cube from 4 squares
c
	    call pushma
		call side
		call rotate(900, 'x')
		call side
		call rotate(900, 'x')
		call side
		call rotate(900, 'x')
		call side
	    call popmat

	call closeo

	end

c
c side
c
c define a face for the cube
c
	subroutine side

	call pushma
	    call transl(0.0, 0.0, 1.0)
	    call rect(-1.0, -1.0, 1.0, 1.0)
	call popmat

	end
