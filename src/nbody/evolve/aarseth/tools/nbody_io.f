C
C subroutines to read particle dump OUT3 files from NBODY1, NBODY2, NBODY5
C              read and write UNIT 4 initial conditions 
C
C used by programs:
C		u3tos, u4tos, stou4
C these are C programs, and as such these are dependant on the 
C CtoF interface.
C
C	7-apr-93  Created		Peter Teuben
C	9-apr-93  changed NB4 to REAL POS()	PJT(+SJA)
C		  changed NB3 to include NK header for A(NK)
C      29-mar-94  allow NB3 to read files without the NK header item
C      23-may-95  removed _ from name for f2c (who otherwise would create
C                 symbols with appending __
C      10-jun-95  changed meaning of nk in NB3HEADER (nk=0 are old datasets)
C      19-jun-97  fixed NB4OPEN
C
C-----------------------------------------------------------------------
      SUBROUTINE NB3OPEN(fname)
      CHARACTER fname*(*)

      OPEN (UNIT=3,STATUS='OLD',FORM='UNFORMATTED',FILE=FNAME)

      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE NB3HEADER(n,model,nrun,nk)
      INTEGER n,model,nrun,nk

      IF (nk.EQ.0) THEN
          READ (3, ERR=99, END=99)  n, model, nrun
      ELSE
          READ (3, ERR=99, END=99)  n, model, nrun, nk
      ENDIF
      RETURN
 99   n = 0
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE NB3DATA(n,nk,a,body,xs,xdot,name)
      INTEGER n,nk,name(n)
      REAL a(nk),body(n),xs(3,n),xdot(3,n)
c
      INTEGER j,k

      IF (n.LE.0 .OR. nk.LE.0) RETURN

      READ (3)  (a(k),k=1,nk), (body(j),j=1,n),
     &           ((xs(k,j),k=1,3),j=1,n), ((xdot(k,j),k=1,3),j=1,n),
     &           (name(j),j=1,n)

      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE NB3CLOSE()

      CLOSE (UNIT=3)

      RETURN
      END
C***********************************************************************
      SUBROUTINE NB4OPEN(fname)
      CHARACTER fname*(*)

      OPEN (UNIT=4,STATUS='UNKNOWN',FORM='UNFORMATTED',FILE=FNAME)

      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE NB4GET(n,body,x,xdot)
      INTEGER n
      REAL body(n), x(3,n), xdot(3,n)

      INTEGER i,k

      DO 1 i = 1,n
         READ (4,ERR=99)  body(i), (x(k,i),k=1,3), (xdot(k,i),k=1,3)
  1   CONTINUE

      RETURN
 99   n=0
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE NB4PUT(n,body,x,xdot)
      INTEGER n
      REAL body(n), x(3,n), xdot(3,n)

      INTEGER i,k

      DO 1 i = 1,n
         WRITE (4)  body(i), (x(k,i),k=1,3), (xdot(k,i),k=1,3)
  1   CONTINUE

      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE NB4CLOSE()

      CLOSE (UNIT=4)

      RETURN
      END
C-----------------------------------------------------------------------
