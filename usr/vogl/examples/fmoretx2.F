	program moret2

#ifdef SGI
#include "fgl.h"
#include "fdevice.h"
#else
#include "fvogl.h"
#include "fvodevice.h"
#endif

	integer *2 val

	call winope('moretxt2', 8)

	call unqdev(INPUTC)
	call qdevic(KEYBD)
c
c Read the initial REDRAW event
c
	idum = qread(val)

	call hfont('futura.l', 8)
	call htexts(0.03, 0.04)

	call ortho2(0.0, 1.0, 0.0, 1.0)


	call color(RED)
	call clear()


	call drawst()

c	 Now do it all with the text rotated .... 

	call htexta(45.0)
	call drawst

	call htexta(160.0)
	call drawst()

	call htexta(270.0)
	call drawst

c	 Now with a single character 

	call htexta(0.0)
	call drwst2()

	call htextang(45.0)
	call drwst2()

	call htextang(160.0)
	call drwst2()

	call htextang(270.0)
	call drwst2()

	call gexit()
	end

	subroutine drawst
#ifdef SGI
#include "fgl.h"
#else
#include "fvogl.h"
#endif

	integer *2 val

	call color(BLACK)
	call rectf(0.1, 0.1, 0.9, 0.9)
	call color(WHITE)
	call move2(0.1, 0.5)
	call draw2(0.9, 0.5)
	call move2(0.5, 0.1)
	call draw2(0.5, 0.9)

	call color(GREEN)
	call move2(0.5, 0.5)
	call hleftj(.true.)
	call hchars('This is Left Justified text', 27)

	idum = qread(val)

	call color(BLACK)
	call rectf(0.1, 0.1, 0.9, 0.9)
	call color(WHITE)
	call move2(0.1, 0.5)
	call draw2(0.9, 0.5)
	call move2(0.5, 0.1)
	call draw2(0.5, 0.9)

	call color(YELLOW)
	call move2(0.5, 0.5)
	call hcente(.true.)
	call hchars('This is Centered text', 21)
	call hcente(.false.)

	idum = qread(val)

	call color(BLACK)
	call rectf(0.1, 0.1, 0.9, 0.9)
	call color(WHITE)
	call move2(0.1, 0.5)
	call draw2(0.9, 0.5)
	call move2(0.5, 0.1)
	call draw2(0.5, 0.9)

	call color(MAGENT)
	call move2(0.5, 0.5)
	call hright(.true.)
	call hchars('This is Right Justified text', 28)
	call hright(.false.)

	idum = qread(val)

	end

	subroutine drwst2
#ifdef SGI
#include "fgl.h"
#else
#include "fvogl.h"
#endif
	integer *2 val

	call color(BLACK)
	call rectf(0.1, 0.1, 0.9, 0.9)
	call color(WHITE)
	call move2(0.1, 0.5)
	call draw2(0.9, 0.5)
	call move2(0.5, 0.1)
	call draw2(0.5, 0.9)

	call color(GREEN)
	call move2(0.5, 0.5)
	call hleftj(.true.)
	call hdrawc('B')

	idum = qread(val)

	call color(BLACK)
	call rectf(0.1, 0.1, 0.9, 0.9)
	call color(WHITE)
	call move2(0.1, 0.5)
	call draw2(0.9, 0.5)
	call move2(0.5, 0.1)
	call draw2(0.5, 0.9)

	call color(YELLOW)
	call move2(0.5, 0.5)
	call hcente(.true.)
	call hdrawc('B')
	call hcente(.false.)

	idum = qread(val)

	call color(BLACK)
	call rectf(0.1, 0.1, 0.9, 0.9)
	call color(WHITE)
	call move2(0.1, 0.5)
	call draw2(0.9, 0.5)
	call move2(0.5, 0.1)
	call draw2(0.5, 0.9)

	call color(MAGENT)
	call move2(0.5, 0.5)
	call hright(.true.)
	call hdrawc('B')
	call hright(.false.)

	idum = qread(val)
	
	end
