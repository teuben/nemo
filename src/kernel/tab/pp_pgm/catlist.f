      PROGRAM CATLIST
C
C     PROGRAM TO LIST CATALOG BETWEEN FREQUENCY LIMITS
C
      CHARACTER*79 LINE
      CHARACTER*64 FILX
      CHARACTER*15 MOLNAM,CFQLOW,CFQHI
      CHARACTER*1  RESP
      INTEGER*4 MOLTAG,NLINE,CATFRQ
      REAL*4 FQLOW,FQHI,Q(7)
      REAL*8 FREQ,C
      LOGICAL CNVT
      C=29979.2458D0
      WRITE(*,*) ' Enter list file'
      READ(*,999) FILX
      OPEN(50,FILE=FILX,STATUS='UNKNOWN')
      WRITE(*,*) MOLNAM,' MHz units, or  wavenumbers ? (m,w) [M] '
      READ(*,999) RESP
      IRESP=ICHAR(RESP)
      IF(IRESP.LT.ICHAR('a') ) IRESP=IRESP+32
      CNVT=IRESP.EQ.ICHAR('w')
      IF(CNVT) THEN
          WRITE(*,*) ' ENTER FREQUENCY LIMITS IN WAVENUMBERS '
          FQLOW=0.
          FQHI=10000.
          READ(*,*) FQLOW,FQHI
          FQLOW=FQLOW*C
          FQHI=FQHI*C
      ELSE
          WRITE(*,*) ' ENTER FREQUENCY LIMITS IN GHz '
          FQLOW=0.
          FQHI=10000.
          READ(*,*) FQLOW,FQHI
          FQLOW=FQLOW*1000.
          FQHI=FQHI*1000.
      ENDIF
      WRITE(CFQLOW,'(F13.4)')FQLOW
      WRITE(CFQHI ,'(F13.4)')FQHI
      MOLTAG=0
  5   IF(NXTDIR(MOLTAG).EQ.0) STOP 'LAST SPECIES READ' 
          CALL CATDIR(MOLTAG,MOLNAM,NLINE,Q,IVER)
          NLINE=CATFRQ(MOLTAG,CFQLOW,LINE)
          IF(NLINE.LE.0) GO TO 5
          IF(LINE.GT.CFQHI) GO TO 5
              WRITE(*,*) MOLNAM,' is next. OK? (y/n/a) [y] '
              READ(*,999) RESP
              IRESP=ICHAR(RESP)
              IF(IRESP.LT.ICHAR('a') ) IRESP=IRESP+32
              IF(IRESP.EQ.ICHAR('a') ) STOP 'ABORTED'
              IF(IRESP.EQ.ICHAR('n') ) GO TO 5
              WRITE(50,'(14X,I6,A)') MOLTAG,MOLNAM
 45              IF(CNVT) THEN
                     READ(LINE,'(F13.4,F8.4)') FREQ,ERR
                     ERR=ERR/C
                     FREQ=FREQ/C
                     WRITE(LINE(1:21),'(F13.5,F8.5)') FREQ,ERR
                 ENDIF
                 WRITE(50,'(A)') LINE
                 NLINE=NLINE+1
                 CALL CATRD(MOLTAG,NLINE,LINE,IERR)
                 IF(IERR.EQ.0 .AND. LINE.LE.CFQHI) GO TO 45
          GO TO 5
 999  FORMAT(A)
      END
