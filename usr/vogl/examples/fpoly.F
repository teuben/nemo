c
c drawpoly
c
c	draw some polygons
c
	subroutine drawpo

	
#ifdef SGI
#include "fgl.h"
#else
#include "fvogl.h"
#endif
	integer *2 val
c
c	An array of points for a polygon
c
	real parray(3,4)
	real vec(3)

	data parray/ -8.0, -8.0, 0.0,
     +  	    -5.0, -8.0, 0.0,
     +  	    -5.0, -5.0, 0.0,
     +  	    -8.0, -5.0, 0.0 /

	call color(YELLOW)

c
c Draw a polygon using poly, parray is our array of
c points and 4 is the number of points in it.
c
	call poly(4, parray)

	call color(GREEN)

c
c Draw a 5 sided figure by using bgnpolygon, v3f and endpolygon
c
	call polymo(PYM_LI)
	call bgnpol
		vec(1) = 0.0 
		vec(2) = 0.0 
		vec(3) = 0.0
		call v3f(vec)

		vec(1) = 3.0 
		vec(2) = 0.0 
		vec(3) = 0.0
		call v3f(vec)
		
		vec(1) = 3.0 
		vec(2) = 4.0 
		vec(3) = 0.0
		call v3f(vec)
		
		vec(1) = -1.0 
		vec(2) = 5.0 
		vec(3) = 0.0
		call v3f(vec)
		
		vec(1) = -2.0 
		vec(2) = 2.0 
		vec(3) = 0.0
		call v3f(vec)
	call endpol

	call color(MAGENT)

c
c draw a sector representing a 1/4 circle
c
	call arc(1.5, -7.0, 3.0, 0, 900)
        
	call move2(1.5, -7.0)
	call draw2(1.5, -4.0)

	call move2(1.5, -7.0)
	call draw2(4.5, -7.0)
 


	idum = qread(val)

	end

c
c drawpolyf
c
c	draw some polygons
c
	subroutine drawpf

#ifdef SGI
#include "fgl.h"
#else
#include "fvogl.h"
#endif
	
	integer *2 val
c
c	An array of points for a polygon
c
	real parray(3,4)

	data parray/ -8.0, -8.0, 0.0,
     +  	    -5.0, -8.0, 0.0,
     +  	    -5.0, -5.0, 0.0,
     +  	    -8.0, -5.0, 0.0 /

	call color(YELLOW)

	call polymo(PYM_FI)
c
c Draw a polygon using poly, parray is our array of
c points and 4 is the number of points in it.
c
	call polf(4, parray)

	call color(GREEN)

        call pmv(0.0, 0.0, 0.0)
                call pdr(3.0, 0.0, 0.0)
                call pdr(3.0, 4.0, 0.0)
                call pdr(-1.0, 5.0, 0.0)
                call pdr(-2.0, 2.0, 0.0)
        call pclos()

        call color(MAGENT)

c
c draw a filled sector representing a 1/4 circle
c
        call arcf(1.5, -7.0, 3.0, 0, 900)

	idum = qread(val)

	end
c
c Using polygons, hatching, and filling.
c
	program fpoly
#ifdef SGI
#include "fgl.h"
#include "fdevice.h"
#else
#include "fvogl.h"
#include "fvodevice.h"
#endif

	integer *2 val

	call winope('fpoly', 5)

c
c	We are interested in keyboard events
c
	call unqdev(INPUTC)
	call qdevic(KEYBD)
c
c Read the initial REDRAW event
c
	idum = qread(val)

c
c load a hershey font
c
	call hfont('futura.l', 8)
c
c clear to black
c
	call color(BLACK)
	call clear

c
c world coordinates are now in the range -10 to 10
c in x, y, and z. Note that positive z is towards us.
c
	call ortho(-10.0, 10.0, -10.0, 10.0, 10.0, -10.0)

	call color(YELLOW)

c
c write out the string "Polygon from poly()" in the
c starting at (-8.0, -4.0) and scaled to be 4.0 units long,
c 0.5 units high.
c
	call hboxte(-8.0, -4.0, 4.0, 0.5,
     +   'Polygon from poly()/ polf()', 27)

	call color(GREEN)

c
c write out a scaled string starting at (0.0, 6.0)
c
	call hboxte(0.0, 6.0, 4.5, 0.5,
     +   'Polygon from bgnpol()/ endpol()', 31)
	call hboxte(0.0, 5.0, 4.5, 0.5,
     +   '             pmv()/ pdr()/ pclos()', 34)

	call color(MAGENT)

c
c draw some polygons
c
	call drawpo

c
c  Rotate 20 degrees around x and 30 around y
c
	call rotate(200, 'x')
	call rotate(300, 'y')

c
c draw some polygons with filling
c
	call drawpf

	call gexit

	end
