	program ftrivial

#ifdef SGI
#include "fgl.h"
#include "fdevice.h"
#else
#include "fvogl.h"
#include "fvodevice.h"
#endif
c 
c  the basic test program for a driver if we can draw a line and do
c  hardware text we are almost there!
c
	integer*2 val

c
c  Open a window (or just the screen in some cases)
c
	call winope('trivial', 7)
c
c  We want to know about keyboard hits....
c
	call qdevic(INPUTC)
	call qdevic(KEYBD)
c
c Read the initial REDRAW event
c
	idum = qread(val)
c
c  we want our coordinate system to run from -1 to 1 in x and y
c
	call ortho2(-1.0, 1.0, -1.0, 1.0)
c
c  set current color to black 
c
	call color(BLACK)
c
c  clear to current color
c
	call clear
c
c  we want to draw in green 
c
	call color(GREEN)
c
c  draw a horizontal line at y = 0 
c
	call move2(-1.0, 0.0)
	call draw2(1.0, 0.0)
c
c  pause for some input 
c
	i = qread(val)
c
c  draw a line along x = 0 
c
	call move2(0.0, 0.0)
	call draw2(0.0, 1.0)
c
c  move text position to the middle of the screen 
c
	call cmov2(0.0, 0.0)
c
c  draw 'Hello' starting at the origin 
c
	call charst('Hello', 5)
c
c  pause again 
c
	idum = qread(val)
c
c  set screen back to original state 
c
	call gexit
	end
