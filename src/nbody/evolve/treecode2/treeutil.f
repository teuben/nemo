C  23-oct-91 pjt  Needed to disable the locally defined SECOND routine
c		  for the Cray, called it MYSECOND for now.
C***********************************************************************
C
C
                          SUBROUTINE ranset(sd)
C
C
C***********************************************************************
C
C
C     Dummy subroutine to initialize random numbers.
C
C
C=======================================================================

        DOUBLE PRECISION sd

        RETURN 
        END

C***********************************************************************
C
C
                           FUNCTION ranf(iran)
C
C
C***********************************************************************
C
C
C     Function to return random numbers.
C
C
C=======================================================================

        DOUBLE PRECISION ranf,drand
        INTEGER iran

        ranf=DRAND(iran)

        RETURN 
        END

C***********************************************************************
C
C
                          SUBROUTINE mysecond(cpu)
C
C
C***********************************************************************
C
C
C     Subroutine to return elapsed cpu time.
C
C
C=======================================================================

        DOUBLE PRECISION cpu
        REAL etime,utime,stime,x

        x=etime(utime,stime)

        cpu=utime+stime     

        RETURN 
        END

C***********************************************************************
C
C
             SUBROUTINE timedate(daytime,date,machine,channel)
C
C
C***********************************************************************
C
C
C     Subroutine to return time of day, etc.
C
C
C=======================================================================

        CHARACTER*8 daytime,date,machine,channel
        CHARACTER*24 dstring

        machine='sun     '
        channel='a       '

        CALL fdate(dstring)

        daytime=dstring(12:19)
        date(1:3)=dstring(5:7)
        date(4:5)=dstring(9:10)
        date(6:6)=' '
        date(7:8)=dstring(23:24)

        RETURN 
        END

