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
C
C
C=======================================================================

c       REAL etime,utime,stime,x,second,ctimes(2)
c       x=etime(utime,stime)
c       second=utime+stime     

        REAL etime,ctimes(2)
	second = etime(ctimes)
        RETURN 
        END
