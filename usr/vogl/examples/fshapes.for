	program fshapes

$include: 'fvogl.h'
$include: 'fvodevice.h'

c
c This program shows some of the simple primitives.
c
	integer edgel
	integer *2 val, vminx, vmaxx, vminy, vmaxy
c
c open a window (or the screen on some devices)
c
	call prefsi(512, 600)
	call winope('shapes', 6)
	call unqdev(INPUTC)
c
c Read the initial REDRAW event
c
	idum = qread(val)

c 
c the two lines below clear the screen to white if we have
c colours, on a monochrome device color is ignored so the
c screen will be cleared to its background color, normally black.
c
	call color(BLACK)
	call clear()

c
c set the screen to be 2.0 units wide and 2.0 units wide, with
c the drawable coordinates going from -1.0 to 1.0.
c
	call ortho2(-1.0, 1.0, -1.0, 1.0)

	call color(MAGENT)

c
c okay, so we want to draw in the range -1 to 1, but we
c only want to draw in the top lefthand corner of the
c screen. The call to viewpo allows us to do this. As
c viewpo always takes screen coordinates, we need to
c call getviewpo to found out how big our screen is
c at the moment. We use the values returned from getviewpo
c calculate the positions for our new viewpo. We note
c that on an Iris (0,0) is the bottom left pixel.
c
	call getvie(vminx, vmaxx, vminy, vmaxy)

	maxx = vmaxx
	minx = vminx
	maxy = vmaxy
	miny = vminy

	call viewpo(minx, (maxx - minx) / 2, (maxy - miny) / 2, maxy)

c
c  write out a heading 
c
	call cmov2(-0.9, -0.5)
	call charst('rect', 4)

c
c draw a rectangle around the points (-0.2, -0.2), (-0.2, 0.2),
c (0.3, 0.2), and (0.3, -0.2).
c
	call rect(-0.2, -0.2, 0.3, 0.2)

	call color(BLUE)

c
c now we want to draw in the top right corner of the screen,
c and we want to draw a circular circle so we must make sure
c our viewpo is square (if it isn't we'll get an ellipse),
c so we calculate the shortest edge and set up a viewpo which
c is in the top right region of the screen, but doesn't necessarilly
c occupy all the top right corner.
c
c find smallest edge 
c
	if (maxx - minx .gt. maxy - miny) then
		edgel = (maxy - miny) / 2
	else 
		edgel = (maxx - minx) / 2
	end if

c
c create a square viewpo - otherwise the circle will look
c like an ellipse.
c

	call viewpo((maxx - minx) / 2, (maxx - minx) / 2 + edgel,
     +           (maxy - miny) / 2, (maxy - miny) / 2 + edgel)

	call cmov2(-0.9, -0.5)
	call charst('circle', 6)

c
c draw a circle of radius 0.4 around the point (0.0, 0.0)
c
	call circ(0.0, 0.0, 0.4)

	call color(GREEN)

c
c bottom left hand corner.
c
	call viewpo(minx, (maxx - minx) / 2,
     +                miny, (maxy - miny) / 2)

	call cmov2(-0.9, -0.5)
	call charst('ellipse', 7)

c
c To draw an ellipse we change the aspect ratio so it is no longer
c 1 and call circ. In this case we use ortho2 to make the square
c viewpo appear to be higher than it is wide. Alternatively you
c could use arc to construct one.
c
c The call to pushmatrix saves the current viewing transformation.
c After the ortho2 has been done, we restore the current viewing
c transformation with a call to popmatrix. (Otherwise everything
c after the call to ortho would come out looking squashed as the
c world aspect ratio is no longer 1).
c
	call pushma
		call ortho2(-1.0, 1.0, -1.0, 2.0)
		call circ(0.0, 0.5, 0.4)
	call popmat

	call color(RED)

c
c bottom right hand corner
c
	call viewpo((maxx - minx) / 2, maxx,
     +                 miny, (maxy - miny) / 2)

	call cmov2(-0.9, -0.5)
	call charst('arc', 3)

c
c draw an arc centered at (0.0, 0.0), radius of 0.4. 0.0 is the start
c angle and 90.0 is the end angle of the arc being drawn. So this
c draws a quarter circle - unless our viewpo isn't square.
c
	call arc(0.0, 0.0, 0.4, 0, 900)

c
c we want to stop after a keyboard event
c
	call qdevic(KEYBD)

c
c wait for the event
c
	idum = qread(val)

	call gexit

	end
