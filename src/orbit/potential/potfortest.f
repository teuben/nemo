	PROGRAM pottest
c
c  test potential descriptors in a fortran environment
c  if your potential(5) was compiled with fortran, this
c  is all you need to link, if it was written in C, you
c  need to link with NEMO's potentialf.c  (see also libnemo.a)
c
	DOUBLE PRECISION pos(3), acc(3), phi, time, par(10)
	INTEGER ndim, npar
	CHARACTER name*40

	ndim=3

	WRITE(6,*) 'Warning: no potential parameters allowed'
	npar=0
	CALL inipotential(npar,par,name)
	time = 0.0d0
	WRITE(6,*) 'Initpotential: par(1) returns ',par(1)

	WRITE(6,*) 'Enter test position: '
	READ(5,*) pos(1), pos(2), pos(3)
	CALL potential(ndim,pos,acc,phi,time)
	WRITE(6,*) 'Pos= ',pos(1),pos(2),pos(3)
	WRITE(6,*) 'Phi= ',phi
	WRITE(6,*) 'Acc= ',acc(1),acc(2),acc(3)

	END
	
