c-----------------------------------------------------------------------
c this curious program shows how local arrays
c (i.e. ones not with a SAVE statement, or in some common block)
c cannot exceed too much space (compiler/machine dependant)
c 
	program big
	write (*,*) "Big"
	call worker1
	end
c-----------------------------------------------------------------------
	subroutine worker1
	include 'big.inc'
c++ comment out any of the 3 lines below to get different behaviors
c	save array
c       common array
c       common /biga/array
c++
	integer i,j,k

	do k=1,NC
	 do j=1,NB
	  do i=1,NA
	   array(i,j,k) = i+j+k
	  enddo
	 enddo
	enddo
	end
c-----------------------------------------------------------------------
