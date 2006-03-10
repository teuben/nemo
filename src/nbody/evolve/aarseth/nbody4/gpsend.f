***
      SUBROUTINE GPSEND
*
*
*       Send all particles to GRAPE-6.
*       ---------------------------
*
      INCLUDE 'common4.h'
      INTEGER i,ls
*
*     do i = ifirst,ntot
*        gpaddr(i) = i - ifirst
*        gpindx(i) = i
*        gpdtj = MIN(step(i),dtk(1))
*        if(DMOD(t0(i),gpdtj).ne.0.d0)then
*           WRITE(3,'(a,2i6)')' ILLEGAL STEP IN GPSEND ',i,name(i)
*           WRITE(3,'(1p,3e16.8)')time,t0(i),t0(i)/step(i)
*           dtmax = dtk(1)
*           CALL dtchck(time,dtmax,dtk(40))
*           WRITE(3,*)'STEPS: OLD NEW ',step(i),dtmax
*           step(i) = dtmax
*           gpdtj = step(i)
*        endif
*        CALL g6_set_j_particle(gpid,gpaddr(i),gpindx(i),t0(i),gpdtj,
*    &                          body(i),d2(1,i),fdot(1,i),f(1,i),
*    &                          x0dot(1,i),x0(1,i))
*     enddo
*
* Initialize TNEXT for new sorting and form list of small steps for body > 0.
      ls = 1
      do i = ifirst,ntot
         tnext(i) = t0(i) + step(i)
         if(step(i).lt.smin.and.body(i).gt.0.0d0.and.ls.lt.40)then
            ls = ls + 1
            lshort(ls) = i
         endif
      enddo
*
* Specify membership of list for KS candidates.
      lshort(1) = ls - 1
*
* Reset GRAPE indicator (ISEND = 0 in routine INTGRT saves new send). 
      isend = 0
*
      RETURN
      END
***
