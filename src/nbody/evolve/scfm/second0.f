C***********************************************************************
C
C
                            FUNCTION second()
C
C
C***********************************************************************
C
C
C     Subroutine to return elapsed cpu time.
C	this is a non-working replacement for systems where second()
C	is *not* available, or it just gets you going
C
C
C=======================================================================

        REAL etime,utime,stime,x,second

        x=0
c	x=etime(utime,stime)

        second=utime+stime     

        RETURN 
        END
