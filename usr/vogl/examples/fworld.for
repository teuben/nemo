c
c most of the things in this program have been done before but it has
c a certain novelty value.
c
	program fworld

$include: 'fvogl.h'
$include: 'fvodevic.h'
	integer *2 val
	integer SPHERE
	real RADIUS, PI
	parameter (RADIUS = 10.0, PI = 3.1415926535, SPHERE = 1)

	call winope('fworld', 6)
	call hfont('futura.m', 8)

	call unqdev(INPUTC)
	call qdevic(SPACEK)
	call qdevic(QKEY)
	call qdevic(ESCKEY)
c
c Read the initial REDRAW event
c
	idum = qread(val)

	call perspe(800, 1.0, 0.001, 50.0)
	call lookat(13.0, 13.0, 8.0, 0.0, 0.0, 0.0, 0)

	call color(BLACK)
	call clear

	call makesp

c
c 	draw the main one in cyan
c
	call color(CYAN)

	call callob(SPHERE)

c
c	draw a smaller one outside the main one in white
c
	call color(WHITE)

	call pushma
		call transl(0.0, -1.4 * RADIUS, 1.4 * RADIUS)
		call scale(0.3, 0.3, 0.3)
		call callob(SPHERE)
	call popmat

c
c	scale the text
c
	call hboxfi(2.0 * PI * RADIUS, 0.25 * RADIUS, 31)

c
c	now write the text in rings around the main sphere
c

	call color(GREEN)
	call showroundtext('Around the world in eighty days ')

	call color(BLUE)
c
c	note: that software text is rotated here as
c	anything else would be whether you use textang
c	or rotate depends on what you are trying to do.
c	Experience is the best teacher here.
c
	call rotate(900, 'x')
	call showroundtext('Around the world in eighty days ')

	call color(RED)
	call rotate(900, 'z')
	call showroundtext('Around the world in eighty days ')

	idum = qread(val)

	call gexit

	end
c
c showroundtext
c
c	draw string str wrapped around a circ in 3d
c
	subroutine showroundtext(str)
	character *(*) str

	real RADIUS
	parameter (RADIUS = 10.0)
	integer j

	inc = 3600 / float(nchars(str))

	j = 1
	do 10 i = 0, 3600, inc
		call pushma
c
c 			find the spot on the edge of the sphere
c 			by making it (0, 0, 0) in world coordinates
c
			call rotate(i, 'y')
			call transl(0.0, 0.0, RADIUS)

			call move(0.0, 0.0, 0.0)

			call hdrawc(str(j:j))
			j = j + 1
		call popmat
10	continue

	end

c
c makesphere
c
c	create the sphere object
c
	subroutine makesp
	integer SPHERE
	parameter (SPHERE = 1)
	parameter(PI = 3.1415926535)
	parameter(RADIUS = 10.0)

	call makeob(SPHERE)

	do 10 i = 0, 1800, 200
		call pushma
			call rotate(i, 'y')
			call circ(0.0, 0.0, RADIUS)
		call popmat
10	continue
	
	call pushma
		call rotate(900, 'x')
		do 20 a = -90.0, 90.0, 20.0
			r = RADIUS * cos(a*PI/180.0)
			z = RADIUS * sin(a*PI/180.0)
			call pushma
				call transl(0.0, 0.0, -z)
				call circ(0.0, 0.0, r)
			call popmat	
20		continue
	call popmat

	call closeo

	return
	end
c
c nchars
c
c	find the number of characters in the string str
c
	integer function nchars(str)
	character *(*) str
	
	do 10 i = len(str), 1, -1
		if (str(i:i) .ne. ' ') then
			nchars = i
			return
		end if
10	continue
	nchars = 0
	return
	end
