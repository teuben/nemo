C***********************************************************************
C
C
                        SUBROUTINE hackload(p0)
C
C
C***********************************************************************
C
C
C     Higher level routine to insert a body, identified by the
C     pointer p0, into the tree.  The local variable p is set equal
C     to p0.  The arrays rp, xp, rmid, and xm are the coordinates of
C     body p0, the integerized coordinates of body p0, the midpoint
C     of the system box, and the integerized form of rmid,
C     respectively.  The local variable newt is a pointer to the new
C     root if an expansion of the system box has been performed.
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

        INTEGER p0,p,xp,xm(1:ndim),subindex
        LOGICAL intcoord

        COMMON/loadcom/p,xp(1:ndim)

        DIMENSION rmid(1:ndim),rp(1:ndim)

C   Initialize pointer p, temporary copy of body coordinates, rp.
C   -------------------------------------------------------------
        p=p0

        DO 5 k=1,ndim
           rp(k)=pos(p,k)
 5      CONTINUE

C-----------------------------------------------------------------------
C   Simulate a DO-WHILE loop--expand box until sufficiently large.
C-----------------------------------------------------------------------

 10     CONTINUE

C-----------------------------------------------------------------------
C   If a body lies outside of the system box, then expand.
C-----------------------------------------------------------------------

        IF(.NOT.intcoord(xp,rp)) THEN
C               --------

C   Compute box midpoint.
C   ---------------------
           DO 20 k=1,ndim
              rmid(k)=rmin(k)+.5*rsize
 20        CONTINUE

C-----------------------------------------------------------------------
C   Decide how to expand box--expand to left if body is left of middle.
C-----------------------------------------------------------------------

           DO 30 k=1,ndim
              IF(pos(p,k).LT.rmid(k)) rmin(k)=rmin(k)-rsize
 30        CONTINUE

C   Double box length, record the expansion.
C   ----------------------------------------
           rsize=2.*rsize
           CALL outbox
C               ------

C=======================================================================
C   If tree already exists, redefine root and make old tree a
C   subcell of the new root.
C=======================================================================

           IF(root.NE.null) THEN
              newt=makecell()
C                  --------
              IF(.NOT.intcoord(xm,rmid))
C                     --------
     &           CALL terror(' error in hackload--rmid out of range')
C                     ------
              k=subindex(xm,imax2)
C               --------
              subp(newt,k)=root
              root=newt
           ENDIF

C-----------------------------------------------------------------------
C   If system box is sufficiently large, load the body into the tree,
C   and return.
C-----------------------------------------------------------------------

        ELSE
           root=loadsub(root,imax2)
C               -------
           RETURN

        ENDIF

        GO TO 10

        END

C***********************************************************************
C
C
                     INTEGER FUNCTION loadsub(q,l)
C
C
C***********************************************************************
C
C
C     Recursive function to perform the task of inserting a body into
C     the tree.  Loadsub attempts to insert the body p into the slot
C     to which q points.  If q is null (i.e. a non-existent node),
C     the body is inserted by pointing loadsub to p.  If q points to
C     a body a new cell is created, the old body is placed at the
C     appropriate location below the new cell, and an attempt is made
C     to insert the body p at the next lower level.  Following
C     completion, loadsub is set to point to the new cell.  If q
C     points to a pre-existing cell, an attempt is made to insert
C     body p below this cell.  A pointer to the pre-existing cell
C     is returned upon completion.
C
C     The argument l identifies the current level.  The local
C     variable c is a pointer to a cell at the current level.  The
C     variable xp is the integerized coordinate vector of body p,
C     and the local variables rq and xq are the coordinates and
C     integerized coordinate vector of body q, if q points to a
C     body.  The local variable s is a pointer to a subcell of the
C     cell c.  The local variable l2 identifies the next lower
C     level (i.e. l2 = l / 2 ).
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

        INTEGER q,c,s,xq,xp,subindex,p
        LOGICAL intcoord

        COMMON/loadcom/p,xp(1:ndim)

        DIMENSION xq(1:ndim),rq(1:ndim)

C-----------------------------------------------------------------------
C   Terminate simulation if tree contains too many levels.
C-----------------------------------------------------------------------

        IF(l.EQ.0) CALL terror(' error in loadsub--out of bits')
C                       ------

C-----------------------------------------------------------------------
C   If q points to an empty slot, insert body by pointing loadsub at p.
C-----------------------------------------------------------------------

        IF(q.EQ.null) THEN
           loadsub=p
           RETURN

C-----------------------------------------------------------------------
C   If q points to a body, create a new cell.
C-----------------------------------------------------------------------

        ELSE IF(q.LE.nbodies) THEN
                c=makecell()
C                 --------

                DO 10 k=1,ndim
                   rq(k)=pos(q,k)
 10             CONTINUE

                IF(.NOT.intcoord(xq,rq))
C                       --------
     &             CALL terror(' error in loadsub--xq out of range')
C                       ------

C   Insert old body below new cell.
C   -------------------------------
                k=subindex(xq,l)
C                 --------
                subp(c,k)=q

C   If q points to a cell, insert p below it.
C   -----------------------------------------
        ELSE
           c=q

        ENDIF

C   Determine where to attempt to insert p.
C   ---------------------------------------
        k=subindex(xp,l)
C         --------
        s=subp(c,k)

C   Initialize level flag for next lower level.
C   -------------------------------------------
        l2=l/2

C   Descend the tree.
C   -----------------
        subp(c,k)=loadsub(s,l2)
C                 -------

C   Return pointer to the new or pre-existing cell.
C   -----------------------------------------------
        loadsub=c

        RETURN
        END

C***********************************************************************
C
C
                    LOGICAL FUNCTION intcoord(xp,rp)
C
C
C***********************************************************************
C
C
C     Function to compute integerized coordinates, xp, given real
C     coordinates, rp.  A value of FALSE is returned if the body
C     lies outside of the system box, otherwise intcoord is TRUE.
C     By convention, if a body coordinate is equal to the minimum
C     value of the corresponding box coordinate, the body is
C     assumed to lie in the box.  If, however, a body coordinate
C     is equal to the maximum value of the corresponding box
C     coordinate, the body is assumed to lie outside of the box.
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

        INTEGER xp

        DIMENSION xp(1:ndim),rp(1:ndim)

C   Assume body is in the box.
C   --------------------------
        intcoord=.TRUE.

C-----------------------------------------------------------------------
C   Loop over all spatial dimensions to check if body is in the box.
C-----------------------------------------------------------------------

        DO 10 k=1,ndim

C   Scale coordinates to [0,1].  Integerize if in range [0,1).
C   ----------------------------------------------------------
           xsc=(rp(k)-rmin(k))/rsize

           IF(xsc.GE.0..AND.xsc.LT.1.) THEN
              xp(k)=INT(imax*xsc)

C   If xsc = 1., body lies outside of box.
C   --------------------------------------
           ELSE
              intcoord=.FALSE.

           ENDIF

 10     CONTINUE

        RETURN
        END

C***********************************************************************
C
C
                     INTEGER FUNCTION subindex(x,l)
C
C
C***********************************************************************
C
C
C     Function to determine which subcell to select (1:nsubcell),
C     given the integerized coordinates of the body, x, and the
C     current level of the tree, l.  By convention, the value
C     returned by subindex can be interpreted as 1 + a coded binary
C     number with ndim bits.  The highest order bit determines
C     whether the x-coordinate is less than (0) or greater than (1)
C     the x-coordinate of the midpoint of the cell.  The next
C     highest order bit stores the same information for the
C     y-coordinate, etc.  Thus, for example, in two dimensions the
C     subcells would be coded as:
C
C
C                          -------------------
C                          |        |        |
C                          | 01 (2) | 11 (4) |
C                          |        |        |
C                          ---------+---------
C                          |        |        |
C                          | 00 (1) | 10 (3) |
C                          |        |        |
C                          -------------------
C
C
C     (The values returned by subindex are indicated in parentheses.)
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

        INTEGER x

        DIMENSION x(1:ndim)

C   Initialize subindex to point to lower left subcell.
C   ---------------------------------------------------
        subindex=1

C-----------------------------------------------------------------------
C   Loop over all spatial dimensions, adjusting subinedx accordingly
C   if beyond midpoint of cell.
C-----------------------------------------------------------------------

        DO 10 k=1,ndim

C=======================================================================
C   WARNING--the following line of code is specific to the CRAY CFT
C   compiler.  Under CRAY CIVIC .AND. should be replaced by .INT.,
C   while on a machine which does not allow bitwise AND operations,
C   IF(...) should be replaced by IF(MOD(x(k)/l,2).EQ.1)... .
C=======================================================================

CPJT        IF((x(k).AND.l).NE.0) subindex=subindex+nindex(k)
        IF(MOD(x(k)/l,2).EQ.1) subindex=subindex+nindex(k)
 10     CONTINUE

        RETURN
        END

C***********************************************************************
C
C
                       INTEGER FUNCTION makecell()
C
C
C***********************************************************************
C
C
C     Function to allocate cell storage.  A pointer to the new cell
C     is returned.
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

C-----------------------------------------------------------------------
C   Terminate simulation if no remaining array space for a new cell.
C-----------------------------------------------------------------------

        IF(incells.GE.ncells)
     &     CALL terror(' error in makecell--no more memory')
C               ------

C   Increment cell usage, initialize pointer to the new cell.
C   ---------------------------------------------------------
        incells=incells+1
        makecell=incells+nbodsmax

C-----------------------------------------------------------------------
C   Zero pointers to subcells of new cell.
C-----------------------------------------------------------------------

        DO 10 i=1,nsubcell
           subp(makecell,i)=null
 10     CONTINUE

        RETURN
        END

C***********************************************************************
C
C
                      SUBROUTINE hackcofm(p0,dsq)
C
C
C***********************************************************************
C
C
C     Recursive subroutine to initialize masses, center of mass
C     coordinates, and (sizes of bodies/cells / tol) ** 2 in the
C     tree.  The argument p0 points to the current body/cell and
C     the local variable r points to its subcells.  The argument
C     dsq is equal to (size of the body/cell / tol) ** 2 for the
C     body/cell p0.
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

        INTEGER p0,r

C=======================================================================
C   If current unit is a body, set sizetol2 to zero and return;
C   otherwise begin mass and center of mass coordinate accumulation,
C   initialize sizetol2 for the subcells, and descend the tree.
C=======================================================================

        IF(p0.LE.nbodies) THEN
           sizetol2(p0)=0.
           RETURN
        ENDIF

C   Zero accumulators for the cell.
C   -------------------------------
        mass(p0)=0.

        DO 10 k=1,ndim
           pos(p0,k)=0.
 10     CONTINUE

C   Assign sizetol2 for cell p0, and compute sizetol2 for subcells.
C   ---------------------------------------------------------------
        sizetol2(p0)=dsq
        dsq4=dsq/4.

C-----------------------------------------------------------------------
C   Compute cell properties as sum of properties of its subcells.
C-----------------------------------------------------------------------

        DO 20 i=1,nsubcell
           r=subp(p0,i)

C   If a subcell is non-empty, compute its properties.
C   --------------------------------------------------
           IF(r.NE.null) THEN
              CALL hackcofm(r,dsq4)
C                  --------

C-----------------------------------------------------------------------
C   Sum properties of subcells to obtain values for cell p0.
C-----------------------------------------------------------------------

              mass(p0)=mass(p0)+mass(r)

              DO 15 k=1,ndim
                 pos(p0,k)=pos(p0,k)+pos(r,k)*mass(r)
 15           CONTINUE

           ENDIF

 20     CONTINUE

C-----------------------------------------------------------------------
C   Normalize center of mass coordinates by total cell mass.
C-----------------------------------------------------------------------
	IF (mass(p0).EQ.0) THEN
	    CALL outerror('Oinks mass(p0) zero...')
            CALL outterm('p0 = ',p0)
	ENDIF
        DO 30 k=1,ndim
           pos(p0,k)=pos(p0,k)/mass(p0)
 30     CONTINUE

        RETURN
        END

C***********************************************************************
C
C
                         SUBROUTINE hackquad(p0)
C
C
C***********************************************************************
C
C
C     Recursive subroutine to initialize quadrupole moments of the
C     cells.  The argument p0 points to the current body/cell and
C     the local variable r points to its subcells.
C
C
C=======================================================================

C-----------------------------------------------------------------------
C   Standard CTSS macro for including the Cliche file.
C-----------------------------------------------------------------------

        include 'treedefs.inc'

        INTEGER p0,r

C=======================================================================
C   If current unit is a body, return; otherwise begin accumulation
C   of the quadrupole moments.
C=======================================================================

        IF(p0.LE.nbodies) RETURN

C   Zero accumulators for the cells.
C   --------------------------------
        DO 10 i=1,nquad
           quad(p0,i)=0.
 10     CONTINUE

C-----------------------------------------------------------------------
C   Compute cell properties as sum of properties of its subcells.
C-----------------------------------------------------------------------

        DO 50 i=1,nsubcell

           r=subp(p0,i)

C   If a subcell is non-empty, compute its properties.
C   --------------------------------------------------
           IF(r.NE.null) THEN

              CALL hackquad(r)
C                  --------

C-----------------------------------------------------------------------
C   Sum properties of subcells to obtain values for cell p0.
C-----------------------------------------------------------------------

              DO 40 m=1,MIN(2,ndim)
                 DO 30 n=m,ndim
                    iquad=(m-1)*(ndim-1)+n
                    quad(p0,iquad)=quad(p0,iquad)+mass(r)*
     &                             (3.*(pos(r,m)-pos(p0,m))*
     &                             (pos(r,n)-pos(p0,n)))

                    IF(m.EQ.n) THEN
                       DO 20 k=1,ndim
                          quad(p0,iquad)=quad(p0,iquad)-mass(r)*
     &                                   (pos(r,k)-pos(p0,k))*
     &                                   (pos(r,k)-pos(p0,k))
 20                    CONTINUE
                    ENDIF

                    IF(r.GT.nbodies) quad(p0,iquad)=quad(p0,iquad)+
     &                                              quad(r,iquad)
 30              CONTINUE
 40           CONTINUE

           ENDIF
 50     CONTINUE

        RETURN
        END

