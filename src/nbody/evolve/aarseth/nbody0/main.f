C
C	Fortran driver for Aarseth's NBODY0 program
C	By picking the fortran MAIN, as opposed to the
C	C nemo_main(), you can run NBODY0 withouth NEMO's
C	I/O routines. You need to the fortran I/O stuff 
C	too (nbodyio.f)
C
C	 1-jul-89	created?	PJT
C	17-may-91	documented	PJT
C
      PROGRAM MAIN0
      CALL NBODY0
      END

