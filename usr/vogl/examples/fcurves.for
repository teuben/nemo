c
c using curves
c
	program fcurve

$include: 'fvogl.h'
$include: 'fvodevic.h'

	character buf*50
	real bezier(4, 4), cardin(4, 4), bsplin(4, 4)
	real geom1(3, 4), geom2(3, 6)
	integer *2 val
c
c curve basis types
c
	data bezier /
     +  	-1.0,	3.0,	-3.0,	1.0,
     +  	3.0,	-6.0,	3.0,	0.0,
     +  	-3.0,	3.0,	0.0,	0.0,
     +  	1.0,	0.0,	0.0,	0.0 
     +  /

	data cardin /
     +  	-0.5,	1.5,	-1.5,	0.5,
     +  	1.0,	-2.5,	2.0,	-0.5,
     +  	-0.5,	0.0,	0.5,	0.0,
     +  	0.0,	1.0,	0.0,	0.0
     +  /

	data bsplin /
     +          -0.166666,     0.5,     -0.5,     0.166666,
     +           0.5,         -1.0,      0.5,     0.0,
     +          -0.5,          0.0,      0.5,     0.0,
     +           0.166666,     0.666666, 0.166666, 0.0
     +  /

c
c Geometry matrix to demonstrate basic spline segments
c
	data geom1 /
     +  	 -180.0, 10.0, 0.0,
     +  	 -100.0, 110.0, 0.0,
     +  	 -100.0, -90.0, 0.0,
     +  	 0.0, 50.0, 0.0
     +  /

c
c Geometry matrix to demonstrate overlapping control points to
c produce continuous (Well, except for the bezier ones) curves
c from spline segments
c
	data geom2 /
     +  	200.0, 480.0, 0.0,
     +  	380.0, 180.0, 0.0,
     +  	250.0, 430.0, 0.0,
     +  	100.0, 130.0, 0.0,
     +  	50.0,  280.0, 0.0,
     +  	150.0, 380.0, 0.0
     +  /


	call winope('fcurves', 7)
c
c We'll use the SPACE bar to go to the next curve...
c
	call unqdev(INPUTC)
	call qdevic(SPACEK)
c
c Read the initial REDRAW event
c
	idum = qread(val)

	call ortho2(-200.0, 400.0, -100.0, 500.0)

	call color(BLACK)
	call clear()

	call color(YELLOW)

c
c label the control points in geom1
c
	do 10 i = 1, 4
	    call cmov2(geom1(1, i), geom1(2, i))
	    write(buf, '(i1)')i
	    call charst(buf, nchars(buf))
10	continue
								 
c
c label the control points in geom2
c
	do 20 i = 1, 6
	    call cmov2(geom2(1, i), geom2(2, i))
	    write(buf, '(i1)')i
	    call charst(buf, nchars(buf))
20	continue

c
c set the number of line segments appearing in each curve to 20
c
	call curvep(20)

c
c define the basis matricies
c
	call defbas(1, bezier)
	call defbas(2, cardin)
	call defbas(3, bsplin)

c
c set the current basis as a bezier basis
c
	call curveb(1)

	call color(RED)

c
c draw a curve using the current basis matrix (bezier in this case)
c and the control points in geom1
c
	call crv(geom1)

	call cmov2(70.0, 60.0)
	call charst('Bezier Curve Segment', 20)

	call cmov2(-190.0, 450.0)
	call charst('Three overlapping Bezier Curves', 31)

c
c curven draws overlapping curve segments according to geom2, the
c number of curve segments drawn is three less than the number of
c points passed, assuming there are a least four points in the
c geometry matrix (in this case geom2). This call will draw 3
c overlapping curve segments in the current basis matrix - still
c bezier.
c
	call crvn(6, geom2)

	idum = qread(val)
c
c	Eat the up event as well...
c
	idum = qread(val)

c
c load in the cardinal basis matrix
c
	call curveb(2)

	call color(MAGENT)

	call cmov2(70.0, 10.0)
	call charst('Cardinal Curve Segment', 22)

c
c plot out a curve segment using the cardinal basis matrix
c
	call crv(geom1)

	call cmov2(-190.0, 400.0)
	call charst('Three overlapping Cardinal Curves', 33)

c
c now draw a bunch of them again.
c
	call crvn(6, geom2)

	idum = qread(val)
c
c	Eat the up event as well...
c
	idum = qread(val)

c
c change the basis matrix again
c
	call curveb(3)

	call color(GREEN)

	call cmov2(70.0, -40.0)
	call charst('Bspline Curve Segment', 21)

c
c now draw our curve segment in the new basis...
c
	call crv(geom1)

	call cmov2(-190.0, 350.0)
	call charst('Three overlapping Bspline Curves', 32)

c
c ...and do some overlapping ones
c
	call crvn(6, geom2)

	idum = qread(val)
c
c	Eat the up event as well...
c
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
