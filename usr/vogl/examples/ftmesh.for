c
c Demonstrate triangular mesh
c
	program fpiston

$include: 'fvogl.h'
$include: 'fvodevic.h'

	logical	dobackface, dofill, dodouble
	character ans*1
	integer *2 val

	call prefsi(300, 300)

	print*,'Backfacing ON or OFF (Y/N)?'
	read(*, '(a)') ans
	dobackface = (ans .eq. 'y' .or. ans .eq. 'Y')

	print*,'Fill the polygons (Y/N)?'
	read(*, '(a)') ans
	dofill = (ans .eq. 'y' .or. ans .eq. 'Y')
c
	dodouble=.true.

 	print*,'double buffer (Y/N)?'
 	read(*, '(a)') ans
 	dodouble = (ans .eq. 'y' .or. ans .eq. 'Y')

	call winope('fpiston', 7)

	if( dodouble) call double
	call gconfi

	call unqdev(INPUTC)
	call qdevic(QKEY)
	call qdevic(ESCKEY)
        call qdevic(REDRAW)
c
c Read the initial REDRAW event
c
	idum = qread(val)
c
	call makecyl

	call polymo(PYM_LI)
	if (dofill) call polymo(PYM_FI)
	if (dobackface) call backfa(.true.)
c
c set up a perspective projection with a field of view of
c 40.0 degrees, aspect ratio of 1.0, near clipping plane 0.1,
c and the far clipping plane at 1000.0.
c
	call perspe(400, 1.5, 0.1, 600.0)
	call lookat(0.0, -6.0, 4., 0.0, 0.0, 0.0, 0)

c
c here we loop back here adnaseum until someone hits a key
c
 10	continue

	  do 20 i = 0, 360, 5

	    call color(BLACK)
	    call clear
            call color(RED)
            H = 1. + cos(2.*3.14159265*i/180.)
	    ixr = 500 - getval(MOUSEX)
	    iyr = 500 - getval(MOUSEY)
	    ixr = ixr * 3
	    iyr = iyr * 3
	    call pushma
	       call rotate(ixr, 'x')
	       call rotate(iyr, 'y')
               call piston(H)
	    call popmat

	    if( dodouble ) call swapbu

	    if (qtest()) then
              itest = qread(idata)
              if( itest .eq. QKEY .or. itest .eq.ESCKEY) then
                 call gexit
                 stop
              endif
	    endif

 20       continue
          goto 10
          end

          subroutine makecyl
c
          parameter ( NTRI=20)
          common /cyl/ cs(NTRI), sn(NTRI)
c
          pi = 2.*asin(1.)
          dphi = 2.*pi/(NTRI-1)
c
          do 10 k=1,ntri
            phi = (k-1)*dphi
            cs(k) = cos(phi)
            sn(k) = sin(phi)
 10       continue
          return
          end
          subroutine piston(H)

$include: 'fvogl.h'
$include: 'fvodevic.h'

          parameter ( NTRI=20)
          common /cyl/ cs(NTRI), sn(NTRI)
          real vec(3)
c
c do the sides
          call bgnqst()
          do 100 k=1,ntri
            vec(1)=cs(k)
            vec(2)=sn(k)
            vec(3)= H
            call v3f(vec)
            vec(3) = 0
            call v3f(vec)
 100      continue
          call endqst()
c do the ends
          call color(CYAN)
          do 220 i=-1,1,2
            HH = H*(I+1)/2 
            call bgntme()
            vec(1)=0
            vec(2)=0
            vec(3)=hh
            call v3f(vec)
            do 200 k=1,ntri
              vec(1) =cs(k)
              vec(2) =i*sn(k)
              call v3f(vec)
              call swaptm()
 200        continue
            call endtme()
 220      continue
         return
          end

