
	program flcube

#ifdef SGI
#include "fgl.h"
#include "fdevice.h"
#else
#include "fvogl.h"
#include "fvodevice.h"
#endif

	integer x, y, but
	integer *2 val

	parameter(TRANS = 20.0, SC = 0.1)

	call prefsi(500, 500)

	call winope('flcube', 6)

	call unqdev(INPUTC)
	call qdevic(SKEY)
	call qdevic(XKEY)
	call qdevic(YKEY)
	call qdevic(ZKEY)
	call qdevic(EQUALK)
	call qdevic(MINUSK)
	call qdevic(ESCKEY)
	call qdevic(QKEY)
c
c Read the initial REDRAW event
c
	idum = qread(val)

	call window(-800.0, 800.0, -800.0, 800.0, -800.0, 800.0)
	call lookat(0.0, 0.0, 1500.0, 0.0, 0.0, 0.0, 0)

	tdir = TRANS
	scal = SC

	nplanes = getpla()
	if (nplanes .eq. 1) call makecu(0)

	call makecu(1)

	call backfa(.true.)
c
c Setup drawing into the backbuffer....
c
	call double
	call gconfi

1	continue
		x = 500 - getval(MOUSEX)
		y = 500 - getval(MOUSEY)
		x = x * 3
		y = y * 3
		call pushma
			call rotate(x, 'y')
			call rotate(y, 'x')
			call color(BLACK)
			call clear
			call callob(3)
			if (nplanes .eq. 1) call callob(2)
		call popmat
		call swapbu

		
		if (qtest()) then
			but = qread(val)
			if (but .eq. XKEY) then
				call transl(tdir, 0.0, 0.0)
			else if (but .eq. YKEY) then
				call transl(0.0, tdir, 0.0)
			else if (but .eq. ZKEY) then
				call transl(0.0, 0.0, tdir)
			else if (but .eq. SKEY) then
				call scale(scal, scal, scal)
			else if (but .eq. MINUSK) then
				tdir = -tdir
			
				if (scal .lt. 1.0) then
					scal = 1.0 + SC
				else
					scal = 1.0 - SC
				end if

			else if (but .eq. EQUALK) then
c
c				we are pretending it's a '+' key
c				we are supposed to see if the shift key is
c				also down - but who could be bothered!
c
				tdir = TRANS
			else if (but .eq. QKEY .or. but .eq. ESCKEY) then
				call gexit
				stop
			end if
c
c			Swallow the UP event...
c
			but = qread(val)
		end if
	goto 1
	end

	subroutine makecu(fill)
#ifdef SGI
#include "fgl.h"
#else
#include "fvogl.h"
#endif
	integer	fill

	call makeob(fill + 2)
		if (fill .ne. 0) then
			call polymo(PYM_FI)
		else
			call polymo(PYM_LI)
			call color(BLACK)
		end if

		call pushma
			call transl(0.0, 0.0, 200.0)
			if (fill .ne. 0) then 
				call color(RED)
				call rectf(-200.0, -200.0, 200.0, 200.0)
			else
				call rect(-200.0, -200.0, 200.0, 200.0)
			end if
		call popmat

		call pushma
			call transl(200.0, 0.0, 0.0)
			call rotate(900, 'y')
			if (fill .ne. 0) then
				call color(GREEN)
				call rectf(-200.0, -200.0, 200.0, 200.0)
			else
				call rect(-200.0, -200.0, 200.0, 200.0)
			end if
		call popmat

		call pushma
			call transl(0.0, 0.0, -200.0)
			call rotate(1800, 'y')
			if (fill .ne. 0) then
				call color(BLUE)
				call rectf(-200.0, -200.0, 200.0, 200.0)
			else
				call rect(-200.0, -200.0, 200.0, 200.0)
			end if
		call popmat

		call pushma
			call transl(-200.0, 0.0, 0.0)
			call rotate(-900, 'y')
			if (fill .ne. 0) then
				call color(CYAN)
				call rectf(-200.0, -200.0, 200.0, 200.0)
			else
				call rect(-200.0, -200.0, 200.0, 200.0)
			end if
		call popmat

		call pushma
			call transl(0.0, 200.0, 0.0)
			call rotate(-900, 'x')
			if (fill .ne. 0) then
				call color(MAGENT)
				call rectf(-200.0, -200.0, 200.0, 200.0)
			else
				call rect(-200.0, -200.0, 200.0, 200.0)
			end if
		call popmat

		call pushma
			call transl(0.0, -200.0, 0.0)
			call rotate(900, 'x')
			if (fill .ne. 0) then
				call color(YELLOW)
				call rectf(-200.0, -200.0, 200.0, 200.0)
			else
				call rect(-200.0, -200.0, 200.0, 200.0)
			end if
		call popmat

	call closeo

	return
	end
