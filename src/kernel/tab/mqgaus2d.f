      SUBROUTINE MQGAUS2D(X,Y,Z,N,XC,YC,PEAK,FWHM,CHISQ,BADFIT)
C
C     This subroutine uses the non-linear Levenberg-Marquardt method 
C     (see "Numerical Recipes", W. H. Press et al., p. 526)
C     to fit a Gaussian to a set of measurements. It returns the offset 
C     of the Gaussian, and can optionally fit for the peak and width of
C     the Gaussian.
C
C     INPUTS:
C
C        X(*)    REAL*4   Array of offset values (arc-minutes)
C        Y(*)    REAL*4   Array of offset values (arc-minutes)
C        Z(*)    REAL*4   Array of measured amplitudes.
C        N       INTEGER  Number of points in X and Y
C        PEAK    REAL*4   For PEAK > 0, will use as initial guess of maximum of
C                         Gaussian and return best fit value. For PEAK < 0 will
C                         use ABS(PEAK) in fit and hold constant.              
C        FWHM    REAL*4   For FWHM > 0, will use as initial guess of width of
C                         Gaussian and return best fit value. For FWHM < 0 will
C                         use ABS(FWHM) in fit and hold constant.              
C     OUTPUTS:         
C                      
C        XC      REAL*4   Best fit offset position in units of X
C        YC      REAL*4   Best fit offset position in units of X
C        PEAK    REAL*4   For PEAK > 0 on input, best fit on output
C        FWHM    REAL*4   For FWHM > 0 on input, best fit on output
C        CHISQ   REAL*4   Chi-squared of fit
C        BADFIT  REAL*4   badfit = -1.0 if no convergence after 10 iterations
C

	implicit none
      integer n, i, lista(4), mfit, ier
      real    x(*), y(*), z(*), peak, xc, yc, fwhm,
     &        a(4), chisq1, chisq, alamda1, alamda,
     &        alpha(4,4), covar(4,4), badfit, fgauss2d
      
      external fgauss2d
      
      badfit = 0.0
      mfit = 2
      lista (1) = 1
      lista (2) = 2
      lista (3) = 0
      lista (4) = 0
      a(1) = 0.0
      a(2) = 0.0
      if (peak .ge. 0) then
         mfit = mfit + 1
         lista(mfit) = 3
         a(3) = abs(peak)
      else 
         a(3) = abs(peak)
      end if
      if (fwhm .ge. 0) then
         mfit = mfit + 1
         lista(mfit) = 4
         a(4) = abs(fwhm)
      else
         a(4) = abs(fwhm)
      end if
      
      alamda  = -1.0
      alamda1 = 0.01
      chisq1  = 1.0E10
      do i = 1, 10
         call mrqmin(x,y,z,n,a,4,lista,mfit,covar,
     &        alpha,4,chisq,fgauss2d,alamda,ier)
         if (ier.NE.0) goto 10
         if (abs(chisq1-chisq).lt.1.E-5) goto 20
         if (alamda .gt. alamda1) goto 10
         chisq1 = chisq
         alamda1 = alamda
      end do
 10   badfit = -1.0
      return 
      
 20   continue
      alamda=0.0
      call mrqmin(x,y,z,n,a,4,lista,mfit,covar,
     &     alpha,4,chisq,fgauss2d,alamda,ier)

      xc   = a(1)
      yc   = a(2)
      peak = a(3)
      fwhm = a(4)
      
      return
      
      end
      
      SUBROUTINE MRQMIN(X,Y,Z,NDATA,A,MA,LISTA,MFIT,
     *     COVAR,ALPHA,NCA,CHISQ,FUNCS,ALAMDA,IER)
      PARAMETER (MMAX=4)
      DIMENSION X(*),Y(*),Z(*),A(*),LISTA(*),
     *     COVAR(NCA,NCA),ALPHA(NCA,NCA),ATRY(MMAX),BETA(MMAX),DA(MMAX)
      IF(ALAMDA.LT.0.)THEN
        KK=MFIT+1
        DO 12 J=1,MA
          IHIT=0
          DO 11 K=1,MFIT
            IF(LISTA(K).EQ.J)IHIT=IHIT+1
11        CONTINUE
          IF (IHIT.EQ.0) THEN
            LISTA(KK)=J
            KK=KK+1
          ELSE IF (IHIT.GT.1) THEN
C     PAUSE 'Improper permutation in LISTA'
          ENDIF
12      CONTINUE
C     IF (KK.NE.(MA+1)) PAUSE 'Improper permutation in LISTA'
        ALAMDA=0.001
        CALL MRQCOF(X,Y,Z,NDATA,A,MA,LISTA,MFIT,ALPHA,BETA,NCA,CHISQ,F
     *UNCS)
        OCHISQ=CHISQ
        DO 13 J=1,MA
          ATRY(J)=A(J)
13      CONTINUE
      ENDIF
      DO 15 J=1,MFIT
        DO 14 K=1,MFIT
          COVAR(J,K)=ALPHA(J,K)
14      CONTINUE
        COVAR(J,J)=ALPHA(J,J)*(1.+ALAMDA)
        DA(J)=BETA(J)
15    CONTINUE
      CALL GAUSSJ(COVAR,MFIT,NCA,DA,1,1,IER)
      IF (IER.NE.0) RETURN
      IF(ALAMDA.EQ.0.)THEN
        CALL COVSRT(COVAR,NCA,MA,LISTA,MFIT)
        RETURN
      ENDIF
      DO 16 J=1,MFIT
        ATRY(LISTA(J))=A(LISTA(J))+DA(J)
16    CONTINUE
      CALL MRQCOF(X,Y,Z,NDATA,ATRY,MA,LISTA,MFIT,COVAR,DA,NCA,CHISQ,FU
     *NCS)
      IF(CHISQ.LT.OCHISQ)THEN
        ALAMDA=0.1*ALAMDA
        OCHISQ=CHISQ
        DO 18 J=1,MFIT
          DO 17 K=1,MFIT
            ALPHA(J,K)=COVAR(J,K)
17        CONTINUE
          BETA(J)=DA(J)
          A(LISTA(J))=ATRY(LISTA(J))
18      CONTINUE
      ELSE
        ALAMDA=10.*ALAMDA
        CHISQ=OCHISQ
      ENDIF
      RETURN
      END

      SUBROUTINE FGAUSS2D(X,Y,A,Z,DZDA)
      REAL*4    DX, DY, DX2, DY2, W2, B
      DIMENSION A(*),DZDA(*)
      DATA B /2.77258872/

      DX = X - A(1)
      DY = Y - A(2)
      DX2= DX**2.
      DY2= DY**2.
      W2 = A(4)**2.
      EX=EXP(-B*(DX2+DY2)/W2)
      Z=A(3)*EX
      DZDA(1) = 2.*B*Z*DX/W2
      DZDA(2) = 2.*B*Z*DY/W2
      DZDA(3) = EX
      DZDA(4) = 2.*B*Z*(DX2+DY2)/(W2*A(4))
      RETURN
      END

      SUBROUTINE COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      DIMENSION COVAR(NCVM,NCVM),LISTA(*)
      DO 12 J=1,MA-1
        DO 11 I=J+1,MA
          COVAR(I,J)=0.
11      CONTINUE
12    CONTINUE
      DO 14 I=1,MFIT-1
        DO 13 J=I+1,MFIT
          IF(LISTA(J).GT.LISTA(I)) THEN
            COVAR(LISTA(J),LISTA(I))=COVAR(I,J)
          ELSE
            COVAR(LISTA(I),LISTA(J))=COVAR(I,J)
          ENDIF
13      CONTINUE
14    CONTINUE
      SWAP=COVAR(1,1)
      DO 15 J=1,MA
        COVAR(1,J)=COVAR(J,J)
        COVAR(J,J)=0.
15    CONTINUE
      COVAR(LISTA(1),LISTA(1))=SWAP
      DO 16 J=2,MFIT
        COVAR(LISTA(J),LISTA(J))=COVAR(1,J)
16    CONTINUE
      DO 18 J=2,MA
        DO 17 I=1,J-1
          COVAR(I,J)=COVAR(J,I)
17      CONTINUE
18    CONTINUE
      RETURN
      END

      SUBROUTINE MRQCOF(X,Y,Z,NDATA,A,MA,LISTA,MFIT,ALPHA,BETA,NALP,
     *                  CHISQ,FUNCS)
      PARAMETER (MMAX=20)
      DIMENSION X(*),Y(*),Z(*),ALPHA(NALP,NALP),BETA(MA),
     *    DYDA(MMAX),LISTA(*),A(MA)
      DO 12 J=1,MFIT
        DO 11 K=1,J
          ALPHA(J,K)=0.
11      CONTINUE
        BETA(J)=0.
12    CONTINUE
      CHISQ=0.
      DO 15 I=1,NDATA
        CALL FUNCS(X(I),Y(I),A,ZMOD,DYDA,MA)
        DY=Z(I)-ZMOD
        DO 14 J=1,MFIT
          WT=DYDA(LISTA(J))
          DO 13 K=1,J
            ALPHA(J,K)=ALPHA(J,K)+WT*DYDA(LISTA(K))
13        CONTINUE
          BETA(J)=BETA(J)+DY*WT
14      CONTINUE
        CHISQ=CHISQ+DY*DY
15    CONTINUE
      DO 17 J=2,MFIT
        DO 16 K=1,J-1
          ALPHA(K,J)=ALPHA(J,K)
16      CONTINUE
17    CONTINUE
      RETURN
      END

      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP,IER)
      PARAMETER (NMAX=50)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      IER = 0
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                 IER = 1
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) IER = 1
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
