C***********************************************************************
C
C
                        FUNCTION cvmgp(y1,y2,y3)
C
C
C***********************************************************************
C
C
C     Function to perform conditional vector merge on positive
C     required by CRAY version.
C
C
C=======================================================================

        DOUBLE PRECISION cvmgp,y1,y2,y3

        IF(y3.GE.0.0) THEN
           cvmgp=y1
        ELSE
           cvmgp=y2
        ENDIF

        RETURN 
        END

C***********************************************************************
C
C
                        FUNCTION cvmgt(y1,y2,y3)
C
C
C***********************************************************************
C
C
C     Function to perform conditional vector merge on true required
C     by CRAY version.
C
C
C=======================================================================

        DOUBLE PRECISION cvmgt,y1,y2
        LOGICAL y3

        IF(y3) THEN
           cvmgt=y1
        ELSE
           cvmgt=y2
        ENDIF

        RETURN 
        END

C***********************************************************************
C
C
                         FUNCTION ismax(n,x,inc)
C
C
C***********************************************************************
C
C
C     Function to locate index of maximum element of a real vector.
C
C
C=======================================================================

        DOUBLE PRECISION x(1),xmax
        INTEGER ismax,n,inc,i

        ismax=1
        xmax=x(1)

        DO 10 i=2,n,inc
           IF(x(i).GT.xmax) THEN
              ismax=i
              xmax=x(i)
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
                          FUNCTION ismin(n,x,inc)
C
C
C***********************************************************************
C
C
C     Function to locate index of minimum element of a real vector.
C
C
C=======================================================================

        DOUBLE PRECISION x(1),xmin
        INTEGER ismin,n,inc,i

        ismin=1
        xmin=x(1)

        DO 10 i=2,n,inc
           IF(x(i).LT.xmin) THEN
              ismin=i
              xmin=x(i)
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
                 FUNCTION isrchigt(n,iarray,inc,itarget)
C
C
C***********************************************************************
C
C
C     Function to return index of first element of iarray greater
C     than itarget, or n+1 if none is found.
C
C
C=======================================================================

        INTEGER isrchigt,n,iarray(1),inc,itarget,i

        isrchigt=n+1

        DO 10 i=1,n,inc
           IF(iarray(i).GT.itarget) THEN
              isrchigt=i
              RETURN
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
                         SUBROUTINE link(dummess)
C
C
C***********************************************************************
C
C
C     Dummy subroutine to handle call to CRAY linker.
C
C
C=======================================================================

        CHARACTER*(*) dummess

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE wheneq(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector equal to itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).EQ.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenfgt(n,array,inc,target,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of a real
C     vector greater than target.
C
C
C=======================================================================

        DOUBLE PRECISION array(1),target
        INTEGER index(1),n,inc,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(array(i).GT.target) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenflt(n,array,inc,target,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of a real
C     vector less than target.
C
C
C=======================================================================

        DOUBLE PRECISION array(1),target
        INTEGER index(1),n,inc,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(array(i).LT.target) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenigt(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector greater than itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).GT.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenile(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector less than or equal to itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).LE.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenilt(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector less than itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).LT.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
            SUBROUTINE whenne(n,iarray,inc,itarget,index,nval)
C
C
C***********************************************************************
C
C
C     Subroutine to return locations of all elements of an integer
C     vector not equal to itarget.
C
C
C=======================================================================

        INTEGER iarray(1),index(1),n,inc,itarget,nval,i

        nval=0

        DO 10 i=1,n,inc
           IF(iarray(i).NE.itarget) THEN
              nval=nval+1
              index(nval)=i
           ENDIF
 10     CONTINUE

        RETURN 
        END

C***********************************************************************
C
C
                          SUBROUTINE second(cpu)
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


