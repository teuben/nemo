	program test
	
c
	call pgbegin(0,'/xw',1,1)
c	call pgenv(0.0,20.0,0.0,20.0,0,1)
	call pgswin(0.0,20.0,0.0,20.0)

	call pgsci(1)

	call pgmove(1.0,1.0)
	call pgdraw(3.0,3.0)
	call pgsci(2)
	call pgdraw(10.0,10.0)
	call pgsci(3)
	call pgdraw(15.0,15.0)
	call pgend

	end
