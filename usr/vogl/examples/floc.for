c
c a routine to demonstrate using locator.
c
	program flocator

$include: 'fvogl.h'
$include: 'fvodevic.h'

	integer bt
	integer *2 x, y, sx, sy
	logical act, curpnt
	integer *2 val, vminx, vmaxx, vminy, vmaxy
c
c Note the declaration of the function locator below
c

	call winope('floc', 4)
c
c Read the initial REDRAW event
c
	idum = qread(val)
        call getvie(vminx, vmaxx, vminy, vmaxy)

	xmin = vminx
	xmax = vmaxx
	ymin = vminy
	ymax = vmaxy

        call ortho2(xmin, xmax, ymin, ymax)

	call color(BLACK)
	call clear

	call color(BLUE)

c
c	draw some axes
c
	y = (vmaxy - vminy) / 2
        call move2s(vminx, y)
        call draw2s(vmaxx, y)

	x = (vmaxx - vminx) / 2
        call move2s(x, vminy)
        call draw2s(x, vmaxy)

	call color(GREEN)

c
c	enable the mouse buttons
c
        call unqdev(INPUTC)
        call qdevic(LEFTMO)
        call qdevic(MIDDLE)


	act = .false.
	curpnt = .false.
c
c	qread waits for a mouse event and getval tells us
c	the valuator's value. In this case it's the X and Y
c	positions of the mouse.	Note: these come back to us in
c	screen coordinates.
c

1	continue
		bt = qread(val)
		sx = getval(MOUSEX)
		sy = getval(MOUSEY)

		if (bt .eq. MIDDLE) then
			call gexit
			stop
		else if (act) then
			act = .false.
			call move2s(sx, sy)
			call draw2s(x, y)
		else 
			act = .true.
			x = sx
			y = sy
		end if

c
c		swallow the up event...
c
		bt = qread(val)
	goto 1

	end
