From lars@astro.Princeton.EDU Fri Apr  6 16:32:35 1990
Received: from Princeton.EDU ([128.112.129.117]) by astro.UMD.EDU  (4.0/5.4WLS)
	id AA11841; Fri, 6 Apr 90 16:32:32 EDT
Received: from astro.Princeton.EDU by Princeton.EDU (5.58+++/2.32/mailrelay)
	id AA15432; Fri, 6 Apr 90 16:33:44 EDT
Received: from rigel (rigel.Princeton.EDU) by astro.Princeton.EDU (4.0/1.100)
	id AA18584; Fri, 6 Apr 90 16:34:46 EDT
From: Lars Hernquist <lars@astro.Princeton.EDU>
Received: by rigel (4.0/Astro-Client)
	id AA11354; Fri, 6 Apr 90 16:34:45 EDT
Date: Fri, 6 Apr 90 16:34:45 EDT
Message-Id: <9004062034.AA11354@rigel>
To: teuben@astro.umd.edu
Status: R

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

        DOUBLE PRECISION ranf
        INTEGER iran,iseed
        REAL ran

        DATA iseed/1234/

        SAVE iseed

        ranf=RAN(iseed)

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
        INTEGER icount,icpu

        DATA icount/1/

        SAVE icount

        IF(icount.EQ.1) THEN
           CALL lib$init_timer()
           icount=0
        ENDIF

        CALL lib$stat_timer(2,icpu)

        cpu=1.e-2*icpu

        RETURN 
        END

C***********************************************************************
C
C
             SUBROUTINE timedate(daytime,vdate,machine,channel)
C
C
C***********************************************************************
C
C
C     Subroutine to return time of day, etc.
C
C
C=======================================================================

        CHARACTER*8 daytime,vdate,machine,channel
        CHARACTER*9 dstring

        machine='vax     '
        channel='a       '

        CALL TIME(daytime)
        
        CALL DATE(dstring)

        vdate(1:3)=dstring(4:6)
        vdate(4:5)=dstring(1:2)
        vdate(6:6)=' '
        vdate(7:8)=dstring(8:9)

        RETURN 
        END


