
c
c a program to demonstrate  double buffering and what happens
c when you hit a clipping plane. Specifying an extra argument
c turns on the filling.
c
	program cube

$include: 'fvogl.h'
$include: 'fvodevic.h'

	character ans*1
	real	t, dt
	integer nplane, r, dr
	integer *2 val
	logical ifill, s

	print*,'Fill the polygons (Y/N)?'
	read(*, '(a)') ans
	ifill = ans .eq. 'y' .or. ans .eq. 'Y'
				  

	call prefsi(300, 300)

	call winope('fcube', 5)

	call unqdev(INPUTC)
	call qdevic(KEYBD)
c
c Read the initial REDRAW event
c
	idum = qread(val)

	call polymo(PYM_LI)
	if (ifill) call polymo(PYM_FI)

	dr = 100
	dt = 0.2

	nplane = getpla()

	call color(BLACK)
	call clear

	call window(-1.5, 1.5, -1.5, 1.5, 9.0, -5.0)
	call lookat(0.0, 0.0, 12.0, 0.0, 0.0, 0.0, 0)

	call backfa(.true.)
c
c Setup drawing into the backbuffer....
c
	call double
	call gconfi

	t = 0.0

	r = 0

 10	continue
	    if (r .ge. 3600) r = 0
	    call color(BLACK)
	    call clear

	    call pushma

	    call transl(0.0, 0.0, t)
	    call rotate(r, 'y')
	    call rotate(r, 'z')
	    call rotate(r, 'x')
	    call color(WHITE)

	    call drawcu(nplane)

	    if (nplane .eq. 1 .and. ifill) then
	        call polymo(PYM_LI)
	        call color(BLACK)
	        call drawcu(nplane)
		if (ifill) call polymo(PYM_FI)
	    endif

	    call popmat

	    t = t + dt
	    if (t .gt. 3.0 .or. t .lt. -18.0) dt = -dt

	    call swapbu

	    s = qtest()
	    if (s) then
	    	call gexit
	    	stop
	    endif

	    r = r + dr

	goto 10

	end

c
c this routine draws the cube, using colours if available
c
	subroutine drawcu(nplane)
	integer nplane

$include: 'fvogl.h'

	real carray(3, 8)
	data carray/
     +     -1.0,  -1.0,   1.0,
     +      1.0,  -1.0,   1.0,
     +      1.0,   1.0,   1.0,
     +     -1.0,   1.0,   1.0,
     +     -1.0,  -1.0,  -1.0,
     +      1.0,  -1.0,  -1.0,
     +      1.0,   1.0,  -1.0,
     +     -1.0,   1.0,  -1.0/

	if (nplane.gt.1) call color(RED)

	call pmv(carray(1,1), carray(2,1), carray(3,1))
	call pdr(carray(1,2), carray(2,2), carray(3,2))
	call pdr(carray(1,3), carray(2,3), carray(3,3))
	call pdr(carray(1,4), carray(2,4), carray(3,4))
	call pclos
	
	if (nplane.gt.1) call color(GREEN)

	call pmv(carray(1,6), carray(2,6), carray(3,6))
	call pdr(carray(1,5), carray(2,5), carray(3,5))
	call pdr(carray(1,8), carray(2,8), carray(3,8))
	call pdr(carray(1,7), carray(2,7), carray(3,7))
	call pclos

	if (nplane.gt.1) call color(YELLOW)

	call pmv(carray(1,2), carray(2,2), carray(3,2))
	call pdr(carray(1,6), carray(2,6), carray(3,6))
	call pdr(carray(1,7), carray(2,7), carray(3,7))
	call pdr(carray(1,3), carray(2,3), carray(3,3))
	call pclos

	if (nplane.gt.1) call color(BLUE)

	call pmv(carray(1,1), carray(2,1), carray(3,1))
	call pdr(carray(1,4), carray(2,4), carray(3,4))
	call pdr(carray(1,8), carray(2,8), carray(3,8))
	call pdr(carray(1,5), carray(2,5), carray(3,5))
	call pclos

	if (nplane.gt.1) call color(MAGENT)

	call pmv(carray(1,3), carray(2,3), carray(3,3))
	call pdr(carray(1,7), carray(2,7), carray(3,7))
	call pdr(carray(1,8), carray(2,8), carray(3,8))
	call pdr(carray(1,4), carray(2,4), carray(3,4))
	call pclos
	
	if (nplane.gt.1) call color(CYAN)

	call pmv(carray(1,1), carray(2,1), carray(3,1))
	call pdr(carray(1,5), carray(2,5), carray(3,5))
	call pdr(carray(1,6), carray(2,6), carray(3,6))
	call pdr(carray(1,2), carray(2,2), carray(3,2))
	call pclos

	end
