c
c display all the hershey fonts and demonstrate textang
c
	program fcirctxt

#ifdef SGI
#include "fgl.h"
#include "fdevice.h"
#else
#include "fvogl.h"
#include "fvodevice.h"
#endif
	character*40 str1, str2, str3, str4, fonts(22)
	character*100 buf
	integer i
	integer *2 val
	data fonts/ 'astrology', 'cursive', 'futura.l',
     +  	'futura.m', 'gothic.eng', 'gothic.ger',
     +  	'gothic.ita', 'greek', 'japanese', 'markers',
     +  	'math.low', 'math.upp', 'meteorology', 'music',
     +  	'cyrillic', 'script', 'symbolic', 'times.g',
     +  	'times.ib', 'times.i', 'times.r', 'times.rb' /

	data str1/ 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
	data str2/ 'abcdefghijklmnopqrstuvwxyz' /
	data str3/ '1234567890+-=!@#$%^&*(){}[]' /
	data str4/ '<>,./?~`\|_BONK,blark' /

	call winope("fcirctxt", 8)

c
c we are interested in Keyboard events...
c
	call unqdev(INPUTC)
	call qdevic(KEYBD)
c
c Read the initial REDRAW event
c
	idum = qread(val)

	call color(BLACK)
	call clear

c
c define the world space
c
	call ortho2(-14.0, 14.0, -14.0, 14.0)

	do 10 i = 1, 22

c
c textang is used to specify the orientation of text. As
c we want the title to come out straight we make sure it is
c zero each time we go through this loop.
c
	    call htexta(0.0)

c
c do the title
c
	    call color(YELLOW)
	    call hfont('futura.m', 8)
	    buf = ' '
	    write(buf, '(''This is Hershey font '',a)') fonts(i)
	    call hboxte(-11.0, 12.0, 20.0, 1.0, buf, 32)

c
c draw a box around the title
c
	    call rect(-11.0, 12.0, 9.0, 13.0)

	    call color(GREEN)

c
c grab a font from the table
c
	    call hfont(fonts(i), nchars(fonts(i)))

c
c show the outer ring
c
	    call htexts(1.5, 1.5)
	    call ShowCi(11.0, str1)

c
c show the second ring
c
	    call htexts(1.3, 1.3)
	    call ShowCi(8.5, str2)

c
c show the third ring
c
	    call htexts(1.1, 1.1)
	    call ShowCi(7.0, str3)

c
c show the inside ring
c
	    call htexts(0.9, 0.9)
	    call ShowCi(5.0, str4)

	    idum = qread(val)
	    if (char(val) .eq. 'q') then
		call gexit
		stop
	    end if

	    call color(BLACK)
	    call clear
10	continue

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
c
c ShowCi
c
c	show a ring of text
c
	subroutine ShowCi(r, str)
	real r
	character*(*) str

	real i, inc, x, y, a, pi
	integer j
	character*1 c
	parameter (pi = 3.1415926535)

	j = 1
	inc = 360.0 / nchars(str)

	do 10 i = 0, 360.0, inc
c
c calculate the next drawing position
c
	    c = str(j:j)
	    x = r * cos(i * pi / 180.0)
	    y = r * sin(i * pi / 180.0)
	    call move2(x, y)
c
c calculate angle for next character
c
	    a = 90.0 + i
c
c set the orientation of the next character
c
	    call htexta(a)
c
c draw the character
c
	    call hdrawc(c)
	    j = j + 1
10	continue

	end
