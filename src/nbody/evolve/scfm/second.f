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

        REAL etime,utime,stime,x,second

        x=etime(utime,stime)

        second=utime+stime     

        RETURN 
        END
