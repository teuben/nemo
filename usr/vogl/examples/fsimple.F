	program fsimple


#ifdef SGI
#include "fgl.h"
#include "fdevice.h"
#else
#include "fvogl.h"
#include "fvodevice.h"
#endif

c  A program showing basic line drawing, text and (if applicable)
c  colour. As none of the projection routines have been called we
c  move and draw in the initial coordinate system -1.0 to 1.0.
c
	character fname*80, p*11
	logical sfont
	integer*2 val
	data p/'Hello world'/

	
	print*,'Use (h)ardware font or (H)ershey font?'
	read(*,'(a)') fname

	sfont = fname(1:1) .eq. 'H'
	
	if (sfont) then
		print*,'Enter a Hershey font name:'
		read(*,'(a)') fname
	end if
	
c
c  set up device 
c
	call winope('simple', 6)
c
c  We want to exit after a keyboard press.
c
	call unqdev(INPUTC)
	call qdevic(KEYBD)
c
c Read the initial REDRAW event
c
	idum = qread(val)
c
c  We want our coordinates to go from -1 to 1 in x and y
c
	call ortho2(-1.0, 1.0, -1.0, 1.0)
c
c  change font to the one input...and set textsize
c
	if (sfont) then
		call hfont(fname, nchars(fname))
		call htexts(0.05, 0.05)
	end if
c
c  set current color
c
	call color(BLACK)
c
c  clear screen to current color 
c
	call clear
c
c
	call color(GREEN)
c
c  2 d move to start where we want to draw a string start 
c  Uses cmov2 for 'hardware' font or move2 for software font
c  then draws string in current color 
c
	if (sfont) then
		call move2(-0.9, 0.9)
		call hchars('A Simple Example', 16)
	else
		call cmov2(-0.9, 0.9)
		call charst('A Simple Example', 16)
	end if
c
c
c  the next four lines draw the x 
c
	call move2(0.0, 0.0)
	call draw2(0.76, 0.76)
	call move2(0.0, 0.76)
	call draw2(0.76, 0.0)

c
c  Now some more text...
c
	if (sfont) then
		call move2(0.0, 0.5)
		call hchars('x done', 6)
		call hchars('next sentence', 13)
		call move2(0.0, 0.1)
		do 10 i = 1, 11
			call hdrawc(p(i:i))
10		continue
	else
		call cmov2(0.0, 0.5)
		call charst('x done', 6)
		call charst('next sentence', 13)
		call cmov2(0.0, 0.1)
		do 20 i = 1, 11
			call charst(p(i:i), 1)
20		continue
	end if

c 
c  the next five lines draw the square
c
c
c 65280 = 0xff00(hex)
c

	call deflin(34, 65280)
	call setlin(34)
	call linewi(2)
	call move2(0.0, 0.0)
	call draw2(0.76, 0.0)
	call draw2(0.76, 0.76)
	call draw2(0.0, 0.76)
	call draw2(0.0, 0.0)
c
c  wait for some input 
c
	idum = qread(val)
c
c  set the screen back to its original state 
c
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
