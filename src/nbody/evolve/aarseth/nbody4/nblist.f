
***
      SUBROUTINE NBLIST(I,H2,NNB)
*
*
*       Neighbour list for body #I.
*       ---------------------------
*
      INCLUDE 'common4.h'
      INTEGER i,nnb
      INTEGER iter,j,k,l,ll,nn,toolarge
      INTEGER ierr,iflag
      REAL*8 h2,rij2,ri2,resc2
*     REAL*8 fgp(3,1),f1gp(3,1),xgp(3,1),vgp(3,1),phgp(1)
*
*     nn = ntot - ifirst + 1
*     iter = 0
      toolarge = 0
*
*       Check for particle well outside the escape boundary. 
      ri2 = (x(1,i) - rdens(1))**2 + (x(2,i) - rdens(2))**2 +
     &                               (x(3,i) - rdens(3))**2
      resc2 = 4.d0*rtide**2
      if(ri2.gt.resc2.or.h2.gt.5.d0) toolarge = 1
*
*       Decide between full N loop on host or fast neighbour list on GRAPE.
 1    continue
*     if((iphase.ne.0.and.time.ne.tadj).or.kz(39).eq.-1.or.
*    &    toolarge.ne.0)then
 999     nnb = 0
*      Save index of all particles inside square distance H2 (skip j=i).
         do j = ifirst,ntot
            rij2 = (x(1,i) - x(1,j))**2
            if(rij2.gt.h2) goto 10
            rij2 = rij2 + (x(2,i) - x(2,j))**2
            if(rij2.gt.h2) goto 10
            rij2 = rij2 + (x(3,i) - x(3,j))**2
            if(rij2.lt.h2)then
               if(j.eq.i) goto 10
               nnb = nnb + 1
               if(nnb.ge.lmax-3)then
                  h2 = h2*0.9d0
                  if(nnb.gt.1000.and.toolarge.le.1)then
                     WRITE(6,*)'Far too many neighbors ',i,nnb,h2
                     toolarge = 2
                  endif
                  goto 999
               endif
               ilist(nnb+1) = j
            endif
 10         continue
         enddo
         if(toolarge.gt.1)then
            WRITE(6,*)'H2 adjusted to ',h2,ri2
            WRITE(6,*)'NNB adjusted to ',nnb,i,name(i)
         endif
*
*       Increase search distance if only one neighbour & body #I identified.
         if(nnb.le.3.and.iter.lt.2)then
*       Adopt new neighbour distance from density fitting formula or 2*H.
            rs2 = (rc**2 + ri2)/FLOAT(nc)**0.66667
            h2 = MAX(rs2,4.d0*h2)
            iter = iter + 1
            toolarge = 0
            goto 999
         endif
*        goto 30
*     endif
*
* Get neighbour list for particle I from GRAPE. 
*     ierr = 0
*     iflag = 0
*20   continue
*     if(iflag.ne.0)then
*        if(iflag.lt.0)then
*           CALL gpwipe(gpid,time)
*           CALL gpsend         
*        endif
*        iflag = 0
*        ierr = ierr + 1
*        if(ierr.ge.2)then
*           WRITE(6,*)' USE HOST FOR NBLIST ',i,name(i),step(i),h2
*           toolarge = 1
*           goto 1
*        endif
*     endif
*       Copy variables for particle i only.
*     do k = 1,1
*        do kk = 1,3
*           xgp(kk,k) = x(kk,i)
*           vgp(kk,k) = xdot(kk,i)
*           fgp(kk,k) = 2.d0*f(kk,i)
*           f1gp(kk,k) = 6.d0*fdot(kk,i)
*        enddo
*        gph2(k) = 0.d0
*        phgp(k) = phi(i)
*        index(k) = i
*     enddo
*     gph2(1) = h2 + 1.0d-10
*     gph2(2) = gph2(1)
*     eps2 = 0.d0
*
*     CALL g6_set_ti(gpid,time)
*     CALL g6calc_firsthalf(gpid,nn,1,index,xgp,vgp,
*    &                      fgp,f1gp,phgp,eps2,gph2)
*       Choose between full neighbour list and closest body.
*     if (KZ(39).EQ.0) then
*        gperr = g6calc_lasthalf(gpid,nn,1,index,xgp,vgp,
*    &                           eps2,gph2,gpacc,gpjerk,gppot)
*     else
*        gperr = g6calc_lasthalf2(gpid,nn,1,index,xgp,vgp,
*    &                            eps2,gph2,gpacc,gpjerk,gppot,nnbindex)
*        ilist(2) = nnbindex(1)
*        nnb = 1
*        goto 50
*     endif
*     if(gperr.ne.0)then
*        WRITE(6,*)'GRAPE HARDWARE ERROR: RESET REQUIRED',time
*        WRITE(6,*)'g6calc_lasthalf returned ',gperr,' in nblist'
*        WRITE(6,*)' IPH N I NMI ',iphase,n,i,name(i)
*        iflag = -1
*        ierr = 2
*        goto 20
*     endif
*
*     gperr = g6_read_neighbour_list(gpid)
*     if(gperr.ne.0)then
*        iflag = gperr
*        if(gperr.lt.0)then
*           WRITE(6,*)'GRAPE HARDWARE ERROR: RESET REQUIRED',time
*           WRITE(6,*)'g6_read_neighbour_list returned ',gperr,
*    &                ' in nblist'
*           ierr = 2
*        else
*           h2 = 0.5d0*h2
*        endif
*        goto 20
*     endif
*     gperr = g6_get_neighbour_list(gpid,0,1024,nnb,ilist(2))
*     if(gperr.ne.0)then
*        iflag = gperr
*        if(gperr.lt.0)then
*           WRITE(6,*)'GRAPE HARDWARE ERROR: RESET REQUIRED',time
*           WRITE(6,*)'g6_read_neighbour_list returned ',gperr,
*    &                ' in nblist'
*           WRITE(6,*)' I NNB H2 ',i,nnb,h2
*           ierr = 2
*        else
*           h2 = 0.5d0*h2
*        endif
*        goto 20
*     endif
*
      if(nnb.gt.128.or.nnb.ge.lmax-3)then
         nwarn = nwarn + 1
         if(nwarn.lt.1000)then
            WRITE(7,22)nnb,iphase,i,time,SQRT(h2),SQRT(ri2)
            CALL flush(7)
         endif
***
*
* Notes from Jun on GRAPE-6 neighbour list and possible limitations 
* (16/02/2001).
*
* In GRAPE-6, unlike GRAPE-4, each chip has its own memory unit. To be
* more precise, each chip has 48 (logical) pipelines, and 16 pipelines
* share the same memory unit. So there are three memory units in a
* chip. The depth of the memory unit is 256 words.
* 
* The 16 chips on one boards calculate the forces on the same set of 48
* particles. Thus, for one chip, the maximum possible number of
* neighbours which can be stored in a GRAPE-6 board is 16x256=4096.
* 
* However, this is of course the upper limit. There are two factors
* which reduce the actual number of neighbours stored. The first one is
* that one memory unit is shared by 16 particles. If the neighbour lists
* have no overlap, average length is reduced to 256. The second one is
* that the list for one set of 16 particles is actually distributed to
* 16 chips, each of which calculates the force (and therefore the
* neighbour list) from its own share of particles. If the neighbours of
* particles somehow concentrate to one chip, the neighbour list memory
* of that chip would overflow, even though others are empty.
* 
* Jun thinks the latter is unlikely to occur, since with current
* implementation of the library the particles are distributed according 
* to their indices (particles 0, 16, 32 ... go to chip 0, 1, 17 ... go to
* chip 1 ...).
* 
* Overall, it's hard to predict in what situation the overflow can occur. 
* If the average length of the neighbour list is 64 or less, Jun believes 
* overflow to be very unlikely. He does not know what would happen if the 
* average length is 128 or larger... 
* 
***
*        toolarge = 1
*        goto 1
      endif
*
 22   FORMAT(' WARNING!    NBLIST   NNB IPHASE I T H R ',
     &                     2i5,i6,f12.5,f8.4,f7.3)
*
*       Increase search distance if only one neighbour & body #I identified.
      if(nnb.le.1.and.iter.lt.2)then
         ri2 = (x(1,i) - rdens(1))**2 + (x(2,i) - rdens(2))**2 +
     &                                  (x(3,i) - rdens(3))**2
*
*       Adopt new neighbour distance from density fitting formula or 2*H.
         rs2 = (rc**2 + ri2)/FLOAT(nc+10)**0.66667
         h2 = MAX(rs2,4.d0*h2)
         nwarn = nwarn + 1
         if(nwarn.lt.1000.and.n.ge.1000)then
            rks = 0.d0
            if(i.gt.n)then
               rks = r(i-n)
            endif
            WRITE(14,25)i,nnb,time,SQRT(h2),SQRT(ri2),rks
         endif
         iter = iter + 1
*        if(h2.gt.5.d0)then
            toolarge = 1
            goto 1
*        else
*           goto 20
*        endif
      endif
*
 25   FORMAT(' WARNING!    NBLIST    I NNB T H RI RKS ',
     &                               i6,i3,f9.2,1p,3e10.2)
*
*       Include warning for exceeding maximum value.
*30   if(nnb.ge.lmax)then
*        nwarn = nwarn + 1
*        rks = 0.d0
*        if(i.gt.n)then
*           rks = r(i-n)
*        endif
*        if(nwarn.lt.1000.and.n.ge.1000)then
*           WRITE(14,35)iphase,nnb,i,time,SQRT(h2),rks
*           CALL flush(14)
*        endif
*        WRITE(6,35)iphase,nnb,i,time,SQRT(h2),rks
*        CALL gpfree
*        STOP
*     endif
*
*35   FORMAT(' WARNING!    IPH NNB > LMAX     NNB I T H RKS ',
*    &                                        i3,2i5,f9.2,1p,2e10.2)
*
*       Increase event counter.
   50 nbcall = nbcall + 1
*
*     WRITE (6,60)  I, NNB, (ILIST(L),L=2,LL)
*  60 FORMAT (' NBLIST     I NNB IL ',I5,I4,20I4,2(/,5X,20I4))
*
      RETURN
      END
***
 
