C***********************************************************************
C
C
                          SUBROUTINE hackgrav
C
C
C***********************************************************************
C
C
C     Subroutine to compute the acceleration components and potential
C     field for the body with index pskip.  This version of treegrav
C     is optimized for a three-dimensional system.
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

C   Initialize counters for number of force terms.
C   ----------------------------------------------
        nterms=1
        nqterms=1

C-----------------------------------------------------------------------
C   Determine which bodies, cells to use in force calculation.
C-----------------------------------------------------------------------

        IF(.NOT.usequad) THEN
           CALL walkmono(root)
C               --------
        ELSE
           CALL walkquad(root)
C               --------
        ENDIF

C-----------------------------------------------------------------------
C   Terminate the simulation if nterms exceeds the maximum allowed.
C-----------------------------------------------------------------------

        nterms=nterms-1
        nqterms=nqterms-1
        IF(nterms.GT.maxnterm)
     &     CALL terror(' error in hackgrav--array overflow')
C               ------

C-----------------------------------------------------------------------
C   Compute monopole contribution to potential and acceleration.
C-----------------------------------------------------------------------

        CALL monopole
C            --------

C-----------------------------------------------------------------------
C   Compute quadrupole contribution to potential and acceleration,
C   if required.
C-----------------------------------------------------------------------

        IF(usequad) CALL quadpole
C                        --------

        RETURN
        END

C***********************************************************************
C
C
                          SUBROUTINE monopole
C
C
C***********************************************************************
C
C
C     Subroutine to compute the monopole contribution to the
C     potential and acceleration components for body pskip if the
C     system is three-dimensional.
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

C=======================================================================
C   Loop over bodies, cells that contribute to potential field and
C   acceleration at pos0.
C=======================================================================

        phi(pskip)=0.
        acc(pskip,1)=0.
        acc(pskip,2)=0.
        acc(pskip,3)=0.

        DO 10 i=1,nterms
           tdrdrinv=1./(drdotdr(i)+eps2)
           phii=pmass(i)*SQRT(tdrdrinv)
           phi(pskip)=phi(pskip)-phii
           phii=phii*tdrdrinv
           acc(pskip,1)=acc(pskip,1)-dr(i,1)*phii
           acc(pskip,2)=acc(pskip,2)-dr(i,2)*phii
           acc(pskip,3)=acc(pskip,3)-dr(i,3)*phii
 10     CONTINUE

C   Remove self-force contribution.
C   -------------------------------
        phi(pskip)=phi(pskip)+mass(pskip)*epsinv

        RETURN
        END

C***********************************************************************
C
C
                          SUBROUTINE quadpole
C
C
C***********************************************************************
C
C
C     Subroutine to determine the quadrupole contribution to the
C     potential and acceleration components if the system is three-
C     dimensional.
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

C-----------------------------------------------------------------------
C   Compute quadrupole contribution to the potential and acceleration.
C-----------------------------------------------------------------------

        DO 10 i=1,nqterms
           tdrdrinv=1./(qdrdotdr(i)+eps2)
           qr5inv=tdrdrinv*tdrdrinv*SQRT(tdrdrinv)
           phiquad=(-.5*((qdr(i,1)*qdr(i,1)-qdr(i,3)*qdr(i,3))*
     &             qquad(i,1)+(qdr(i,2)*qdr(i,2)-qdr(i,3)*qdr(i,3))*
     &             qquad(i,4))-(qdr(i,1)*qdr(i,2)*qquad(i,2)+
     &             qdr(i,1)*qdr(i,3)*qquad(i,3)+qdr(i,2)*qdr(i,3)*
     &             qquad(i,5)))*qr5inv
           phi(pskip)=phi(pskip)+phiquad
           phiquad=5.*phiquad*tdrdrinv
           acc(pskip,1)=acc(pskip,1)+qdr(i,1)*phiquad+(qdr(i,1)*
     &                  qquad(i,1)+qdr(i,2)*qquad(i,2)+qdr(i,3)*
     &                  qquad(i,3))*qr5inv
           acc(pskip,2)=acc(pskip,2)+qdr(i,2)*phiquad+(qdr(i,2)*
     &                  qquad(i,4)+qdr(i,1)*qquad(i,2)+qdr(i,3)*
     &                  qquad(i,5))*qr5inv
           acc(pskip,3)=acc(pskip,3)+qdr(i,3)*phiquad+(qdr(i,3)*
     &                  (-qquad(i,1)-qquad(i,4))+qdr(i,1)*qquad(i,3)+
     &                  qdr(i,2)*qquad(i,5))*qr5inv
 10     CONTINUE

        RETURN
        END

C***********************************************************************
C
C
                         SUBROUTINE walkmono(p0)
C
C
C***********************************************************************
C
C
C     Recursive subroutine to walk through the tree and initialize
C     the arrays for the force evaluation, if the system is three-
C     dimensional and only the monopole term is used in the
C     acceleration and potential calculations.  The argument p0 is
C     a pointer to the body/cell under examination.
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

        INTEGER p0

C-----------------------------------------------------------------------
C   Compute displacement vector and distance**2 from body pskip
C   to body/cell p0.
C-----------------------------------------------------------------------

        dr(nterms,1)=pos0(1)-pos(p0,1)
        dr(nterms,2)=pos0(2)-pos(p0,2)
        dr(nterms,3)=pos0(3)-pos(p0,3)
        drdotdr(nterms)=dr(nterms,1)*dr(nterms,1)+dr(nterms,2)*
     &                  dr(nterms,2)+dr(nterms,3)*dr(nterms,3)

C-----------------------------------------------------------------------
C   If tolerance criterion is not satisfied, subdivide the cell and
C   examine all non-empty subcells; otherwise save a pointer to
C   body/cell p0 and the quantities required for the potential.
C-----------------------------------------------------------------------

        IF(drdotdr(nterms).LT.sizetol2(p0)) THEN

           IF(subp(p0,1).NE.0) CALL walkmono(subp(p0,1))
C                                   --------
           IF(subp(p0,2).NE.0) CALL walkmono(subp(p0,2))
C                                   --------
           IF(subp(p0,3).NE.0) CALL walkmono(subp(p0,3))
C                                   --------
           IF(subp(p0,4).NE.0) CALL walkmono(subp(p0,4))
C                                   --------
           IF(subp(p0,5).NE.0) CALL walkmono(subp(p0,5))
C                                   --------
           IF(subp(p0,6).NE.0) CALL walkmono(subp(p0,6))
C                                   --------
           IF(subp(p0,7).NE.0) CALL walkmono(subp(p0,7))
C                                   --------
           IF(subp(p0,8).NE.0) CALL walkmono(subp(p0,8))
C                                   --------

        ELSE

           pmass(nterms)=mass(p0)
           nterms=nterms+1

        ENDIF

        RETURN
        END

C***********************************************************************
C
C
                         SUBROUTINE walkquad(p0)
C
C
C***********************************************************************
C
C
C     Recursive subroutine to walk through the tree and initialize
C     the arrays for the force evaluation, if the system is three-
C     dimensional and both monopole and quadrupole terms are used in
C     the acceleration and potential calculations.  The argument p0
C     is a pointer to the body/cell under examination.
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

        INTEGER p0

C-----------------------------------------------------------------------
C   Compute displacement vector and distance**2 from body pskip
C   to body/cell p0.
C-----------------------------------------------------------------------

        dr(nterms,1)=pos0(1)-pos(p0,1)
        dr(nterms,2)=pos0(2)-pos(p0,2)
        dr(nterms,3)=pos0(3)-pos(p0,3)
        drdotdr(nterms)=dr(nterms,1)*dr(nterms,1)+dr(nterms,2)*
     &                  dr(nterms,2)+dr(nterms,3)*dr(nterms,3)

C-----------------------------------------------------------------------
C   If tolerance criterion is not satisfied, subdivide the cell and
C   examine all non-empty subcells; otherwise save a pointer to
C   body/cell p0 and the quantities required for the potential.
C-----------------------------------------------------------------------

        IF(drdotdr(nterms).LT.sizetol2(p0)) THEN

           IF(subp(p0,1).NE.0) CALL walkquad(subp(p0,1))
C                                   --------
           IF(subp(p0,2).NE.0) CALL walkquad(subp(p0,2))
C                                   --------
           IF(subp(p0,3).NE.0) CALL walkquad(subp(p0,3))
C                                   --------
           IF(subp(p0,4).NE.0) CALL walkquad(subp(p0,4))
C                                   --------
           IF(subp(p0,5).NE.0) CALL walkquad(subp(p0,5))
C                                   --------
           IF(subp(p0,6).NE.0) CALL walkquad(subp(p0,6))
C                                   --------
           IF(subp(p0,7).NE.0) CALL walkquad(subp(p0,7))
C                                   --------
           IF(subp(p0,8).NE.0) CALL walkquad(subp(p0,8))
C                                   --------

        ELSE

           pmass(nterms)=mass(p0)
           IF(p0.GT.nbodies) THEN
              qquad(nqterms,1)=quad(p0,1)
              qquad(nqterms,2)=quad(p0,2)
              qquad(nqterms,3)=quad(p0,3)
              qquad(nqterms,4)=quad(p0,4)
              qquad(nqterms,5)=quad(p0,5)
              qdr(nqterms,1)=dr(nterms,1)
              qdr(nqterms,2)=dr(nterms,2)
              qdr(nqterms,3)=dr(nterms,3)
              qdrdotdr(nqterms)=drdotdr(nterms)
              nqterms=nqterms+1
           ENDIF
           nterms=nterms+1
        ENDIF

        RETURN
        END

