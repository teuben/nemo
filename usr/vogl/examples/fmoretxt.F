c
c drawgr
c
c	draw a grid in the middle of the screen
c
	subroutine drawgr

#ifdef SGI
#include "fgl.h"
#else
#include "fvogl.h"
#endif

	real x

	call color(GREEN)

	call rect(0.1, 0.4, 0.9, 0.6)

	x = 0.2

	do 10 i = 1, 8
		call move2(x, 0.4)
		call draw2(x, 0.6)
		x = x + 0.1
10	continue

	call move2(0.1, 0.5)
	call draw2(0.9, 0.5)

	call color(YELLOW)

	end

c
c demonstrate some more features of text
c
	program	fmtxt

#ifdef SGI
#include "fgl.h"
#include "fdevice.h"
#else
#include "fvogl.h"
#include "fvodevice.h"
#endif

	character fname*80
	integer *2 val

	print*,'Enter a font name:'
	read(*,'(a)') fname

	call prefsi(512, 512)
	call winope('fmoretxt', 8)

c
c We are interested in Keyboard events...
c
	call unqdev(INPUTC)
	call qdevic(KEYBD)
c
c Read the initial REDRAW event
c
	idum = qread(val)

	call hfont(fname, nchars(fname))

	call ortho2(0.0, 1.0, 0.0, 1.0)

	call color(BLACK)
	call clear

	call drawgr

c
c show some scaled text on the grid (In the bottom part)
c
	call hboxte(0.1, 0.4, 0.8, 0.1,
     +      '{This is Some text] | $', 23)

	idum = qread(val)

	call color(BLACK)
	call clear

	call drawgr

c
c centertext causes text to be centered around the current graphics
c position this is especially usefull if you want your text to come
c out centered on a line, or a character to be centered on a point
c in a graph. A non-zero (.true.) argument turns centertext on.
c
c show a string centered on the center line
c
	call hcente(.true.)

	call hboxte(0.5, 0.5, 0.8, 0.1,
     +      '{This is Some Centered text] | $', 31)
c
c turn centertext off. We use an argument with the value zero (.false.).
c
	call hcente(.false.)

	idum = qread(val)

	call color(BLACK)
	call clear

c
c rotate the grid so that it is the same angle as the text after
c textang for text ang.
c
	call pushma
		call transl(0.5, 0.5, 0.0)
		call rotate(900, 'z')
		call transl(-0.5, -0.5, 0.0)

		call drawgr
	call popmat

c
c turn on centered text again
c
	call hcente(.true.)

c
c set the angle to 90.
c
	call htexta(90.0)

c
c draw the string
c
	call htexts(0.05, 0.05)
	call hboxte(0.5, 0.5, 0.8, 0.1,
     +      '{This is Some Rotated Centered text] | $', 39)
c
c turn off center text
c
	call hcente(.false.)

c
c set text angle back to 90
c
	call htexta(0.0)

	idum = qread(val)

	call color(BLACK)
	call clear

	call drawgr

c
c as all the software fonts are proportionally spaced we use
c the fixedwidth call to make each character take the same amount
c of horizontal space. As with centertext this is done by passing
c fixedwidth a non-zero (.true.) argument.
c
	call hfixed(.true.)

	call hboxte(0.1, 0.5, 0.8, 0.1,
     +      '{This is Some Fixedwidth text] | $', 33)

	idum = qread(val)

	call color(BLACK)
	call clear

	call drawgr

c
c now try centered and fixewidth at the same time
c
	call hcente(.true.)

	call move2(0.5, 0.5)
	call hchars('{This is Some Cent.Fixedwidth text] | $', 38)

	call hcente(.false.)
	
	idum = qread(val)
	call color(BLACK)
	call clear

	call drawgr

c
c scale the text so tha a character is the size of a box in
c the grid.
c
	call hboxfi(0.8, 0.1, 8)

c
c draw the two strings fixedwidth (it is still turned on)
c
	call move2(0.1, 0.4)
	call hchars('ABCDefgh', 8)

	call move2(0.1, 0.5)
	call hchars('IJKLmnop', 8)

	idum = qread(val)

	call gexit

	end
c
c nchars
c
c return the real length of a string padded with blanks
c
	integer function nchars(str)
	character *(*) str

	do 10 i = len(str), 1, -1
	    if (str(i:i) .ne. ' ') then
	    	nchars = i
	    	return
	    end if
10      continue

	nchars = 0

	return

	end
