c
c Demonstrate a rotating translating tetrahedron, and 
c doublebuffering.
c
	program ftetra

#ifdef SGI
#include "fgl.h"
#include "fdevice.h"
#else
#include "fvogl.h"
#include "fvodevice.h"
#endif
	integer TETRAHEDRON
	parameter (TETRAHEDRON = 1)

	real R, tx, tz, zeye
	integer rotval, drotval
	logical	dobackface, dofill
	character ans*1
	integer *2 val

	call prefsi(300, 300)

	print*,'Backfacing ON or OFF (Y/N)?'
	read(*, '(a)') ans
	dobackface = (ans .eq. 'y' .or. ans .eq. 'Y')

	print*,'Fill the polygons (Y/N)?'
	read(*, '(a)') ans
	dofill = (ans .eq. 'y' .or. ans .eq. 'Y')

	call winope('ftetra', 6)

	call double
	call gconfi

	call unqdev(INPUTC)
	call qdevic(QKEY)
	call qdevic(ESCKEY)
c
c Read the initial REDRAW event
c
	idum = qread(val)
c
c Make the tetrahedral object
c
	call makeit

	rotval = 0
	drotval = 10
	zeye = 5.0

	R = 1.6

	tx = 0.0
	tz = R

	call polymo(PYM_LI)
	if (dofill) call polymo(PYM_FI)
	if (dobackface) call backfa(.true.)

c
c set up a perspective projection with a field of view of
c 40.0 degrees, aspect ratio of 1.0, near clipping plane 0.1,
c and the far clipping plane at 1000.0.
c
	call perspe(400, 1.0, 0.001, 15.0)
	call lookat(0.0, 0.0, zeye, 0.0, 0.0, 0.0, 0)


c
c here we loop back here adnaseum until someone hits a key
c
 10	continue

	  rotval = 0

	  do 20 i = 0, int(3590 / drotval)

	    call color(BLACK)
	    call clear

c
c Rotate the whole scene...(this acumulates - hence
c drotval)
c
	    call rotate(drotval, 'x')
	    call rotate(drotval, 'z')

	    call color(RED)
	    call pushma
		call rotate(900, 'x')
		call circ(0.0, 0.0, R)
	    call popmat

	    call color(BLUE)
	    call move(0.0, 0.0, 0.0)
	    call draw(tx, 0.0, tz)
			
c
c Remember! The order of the transformations is
c the reverse of what is specified here in between
c the pushmatrix and the popmatrix. These ones don't
c accumulate because of the push and pop.
c

	    call pushma
		call transl(tx, 0.0, tz)
		call rotate(rotval, 'x')
		call rotate(rotval, 'y')
		call rotate(rotval, 'z')
		call scale(0.4, 0.4, 0.4)
 		call callob(TETRAHEDRON)
	    call popmat

	    tz = R * cos(rotval * 3.1415926535 / 180)
	    tx = R * sin(rotval * 3.1415926535 / 180)

	    call swapbu

	    if (qtest()) then
		call gexit
		stop
	    endif

	    rotval = rotval + drotval
	    if (rotval .gt. 3600) rotval = 3600

 20	  continue

	goto 10
		
	end

c
c maketheobject
c
c	generate a tetrahedron object as a series of move draws
c
	subroutine makeit

#ifdef SGI
#include "fgl.h"
#else
#include "fvogl.h"
#endif
	integer TETRAHEDRON, NSIDES, NFACES, NPNTS
	parameter (TETRAHEDRON = 1, NSIDES = 3, NFACES = 4, NPNTS = 4)

	integer colface(NFACES)

	real points(3, NPNTS), tmp(3)

	integer	faces(NSIDES, NFACES)

	integer i, j
	real x, y, z


cdata points/
c+	-0.5, 0.866, -0.667,
c+	-0.5, -0.866, -0.667,
c+	 1.0, 0.0, -0.667,
c+	 0.0, 0.0, 1.334/

      data points/
     +	-0.5, 0.866, -0.667,
     +	-0.5, -0.866, -0.667,
     +	 1.0, 0.0, -0.667,
     +	 0.0, 0.0, 1.334/


	data colface/GREEN, YELLOW, CYAN, MAGENT/

	data faces/
     +	3, 2, 1,
     +	1, 2, 4,
     +	2, 3, 4,
     +	3, 1, 4/

 	call makeob(TETRAHEDRON)

       do 20 i = 1, NFACES
                call color(colface(i))
                call bgnpol
                do 10 j = 1, NSIDES
                    call v3f(points(1, faces(j, i)))
 10             continue
                call endpol
 20     continue

 	call closeo
	end

